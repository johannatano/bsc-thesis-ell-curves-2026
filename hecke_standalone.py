from sage.all import *
from utils.common import Colors
from utils.class_nr import getClassNumber
import argparse
from typing import Optional, List, Dict, Tuple, Set, Any
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from sage.all import CuspForms, Gamma1

from fractions import Fraction
from dataclasses import dataclass


@dataclass
class HeckeResult:
    p: int              # prime p
    n: int              # extension degree n
    q_equiv_N: int       # q % N
    p_equiv_N: int       # p % N
    diff: int            # -T - sage_T
    phi_N: int           # euler_phi(N)
    error_1: int        # cusp_term(q, N, k)
    error_2: int   # higher n,
    legendre_N_q: int         # legendre symbol (q/N)
    legendre_q_N: int         # legendre symbol (N/q)
    legendre_p_N: int  # legendre symbol (p/N)
    legendre_N_p: int  # legendre symbol (N/p),
    max_r: int          # maximum r value encountered (for debugging)
    NSS: int            # number of supersingular curves encountered (for debugging)
    j0_SS: bool         # whether j0 is supersingular (for debugging)
    j1728_SS: bool      # whether j1728 is supersingular (for debugging)
    cusp_dim: int       # dimension of the cusp space (for debugging)


def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument(
        "-p", "--p", type=int, required=False, default=-1, help="Field char p"
    )
    p.add_argument(
        "-n",
        "--n",
        type=int,
        required=False,
        default=1,
        help="Field extension degree n",
    )
    p.add_argument("-N", "--N", type=int, required=False, default=-1, help="Level N")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")

    p.add_argument(
        "--pmax", type=int, default=100, help="Upper prime bound for interactive mode"
    )
    
    p.add_argument(
        "--compare", type=bool, default=False, help="Compare to Sage trace (for debugging)"
    )
    
    p.add_argument(
        "--plist",
        type=int,
        nargs="*",
        default=None,
        help="List of specific primes to process (for debugging)",
    )
    return p.parse_args()


from math import gcd, isqrt


def divisors(n):
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)


def A3_cusp_Gamma1(n, N, k):
    """Cusp/hyperbolic term in tr(T_n) on S_k(Gamma_1(N))."""
    total = 0
    sqrt_n = isqrt(n)
    phiN = euler_phi(N)
    for d in divisors(n):
        if d > sqrt_n:
            break
        e = n // d
        diff = e - d
        inner = 0
        for c in divisors(N):
            g = gcd(c, N // c)
            if (diff == 0 or diff % g == 0) and gcd(d, c) == 1 and gcd(e, N // c) == 1:
                inner += 1
        contribution = d ** (k - 1) * phiN * inner
        if d * d == n:
            contribution = contribution // 2
        total += contribution
    return total


def num_supersingular_curves_q_square(p: int) -> int:
    return (
        (p + 6 - 4 * kronecker(-3, p) - 3 * kronecker(-4, p)) // 6
        + 2
        - 2 * kronecker(-3, p)
        + 1
        - kronecker(-4, p)
    )


def quaternion_class_number(p):
    chi3 = kronecker(-3, p)
    chi4 = kronecker(-4, p)
    H_p = (p + 6 - 4 * chi3 - 3 * chi4) // 12
    n_j0 = (1 - chi3) // 2
    n_j1728 = (1 - chi4) // 2
    n_generic = H_p - n_j0 - n_j1728
    
    
    #print(f"Quaternion class number for p={p}: H_p={H_p}, n_generic={n_generic}, n_j0={n_j0}, n_j1728={n_j1728}")
    # return (p-1) / 12
    test = Fraction(int(n_generic)) + Fraction(int(n_j0), 3) + Fraction(int(n_j1728), 2)
    
    #print(((p-1) / 12), n_generic + n_j0 / 3 + n_j1728 / 2, test)
    
    #if test != Fraction(int(p - 1), 12):
    #    print(f"{Colors.FAIL}Warning: quaternion class number formula mismatch for p={p}: computed {test}, expected {(p-1)/12}{Colors.ENDC}")

    return Fraction(int(p - 1), 12)
    #return n_generic + n_j0 / 3 + n_j1728 / 2
    return Fraction(int(n_generic)) + Fraction(int(n_j0), 3) + Fraction(int(n_j1728), 2)


def get_j_invariants_from_order(D: int, f: int, F) -> List:
    """Return the j-invariants attached to an order via its Hilbert class polynomial.

    Args:
        D: Discriminant of the order D_K f^2.
        f: Conductor of the order. Included for API clarity.
        fq: Finite field data.

    Returns:
        Roots of the Hilbert class polynomial in F_q.
    """
    j_invs = []
    try:
        H = hilbert_class_polynomial(D)
        H_fq = H.change_ring(F)
        # Find roots of H(x) in F_q
        for j, m in H_fq.roots(multiplicities=True):
            for _ in range(m):
                j_invs.append(j)
    except Exception as e:
        print(f"Warning: Could not compute HCP for D={D}: {e}")
    return j_invs


def fundamental_discr(D):
    """Return the fundamental discriminant Δ and conductor f such that D = f² · Δ."""
    if D == 0:
        return 0, 0

    sign = -1 if D < 0 else 1
    abs_D = abs(D)

    # Remove all square factors to find squarefree part
    sf = 1
    n = abs_D
    d = 2
    while d * d <= n:
        exp = 0
        while n % d == 0:
            n //= d
            exp += 1
        if exp % 2 == 1:
            sf *= d
        d += 1
    if n > 1:
        sf *= n
    sf *= sign  # restore sign

    # Fundamental discriminant: sf if sf ≡ 1 mod 4, else 4·sf
    delta = sf if sf % 4 == 1 else 4 * sf
    f = math.isqrt(abs(D // delta))
    return delta, f


def num_P(level, f_pi, f, N_pts, q):
    if level == 1:
        return 1
    result = 1
    for l, a in factor(ZZ(level)):

        h = max(0, ZZ(f_pi).valuation(l) - ZZ(f).valuation(l)) if f_pi * f > 0 else None
        v_q1 = ZZ(q - 1).valuation(l)
        v_N = ZZ(N_pts).valuation(l)

        e1 = min(h, v_q1, v_N // 2) if h is not None else min(v_q1, v_N // 2)
        e2 = v_N - e1
        s1 = min(a, e1)
        s2 = min(a, e2)
        # exact-order-l^a count in Z/l^s1 x Z/l^s2
        result *= l ** (s1 + s2) - l ** (min(a - 1, s1) + min(a - 1, s2))
    return int(result)


def num_C(order_sage, t):
    return int(order_sage.class_number() if t != 0 else 2*order_sage.class_number())

def hk_eval(t, k):
    return int(t) if k == 1 else 1


class Hk:
    """Helpers for the symmetric polynomials appearing in Hecke computations."""

    @staticmethod
    def construct(k: int):
        """Build the polynomial h_k(X,Y)=sum_{i=0}^k X^{k-i}Y^i."""
        R = PolynomialRing(ZZ, ["X", "Y"])
        X, Y = R.gens()
        terms = [X ** (k - i) * Y**i for i in range(k + 1)]
        return R(sum(t for t in terms))

    @staticmethod
    def dickson_recursive(k: int, t: int, q: int) -> int:
        """Evaluate the Dickson-style recurrence for `h_k` without building the full polynomial."""
        # Base cases for h_k
        if k == 0:
            return 1
        if k == 1:
            return t
        # Recurrence: h_k = t * h_{k-1} - n * h_{k-2}
        # This avoids building the huge polynomial object
        hk_prev, hk_curr = 1, t
        for _ in range(k - 1):
            hk_prev, hk_curr = hk_curr, t * hk_curr - q * hk_prev
        return hk_curr

    def eval(fx, t:int, k: int):
        hk_symbolic = Hk.construct(k)
        R = fx.parent()
        x = R.gen()
        frob = x
        frob_dual = t - x
        result_multi = hk_symbolic.subs(X=frob, Y=frob_dual)
        result_uni = R(result_multi)
        final_value = result_uni.quo_rem(fx)[1]
        return final_value


def h_order(f, h_OK, D_K):
    """h(O_f) for order of conductor f."""
    
    if f == 1:
        return h_OK
    
    w = 3 if D_K == -3 else (2 if D_K == -4 else 1)
    
    result = f * h_OK
    n = f
    l = 2
    while l * l <= n:
        if n % l == 0:
            result = result * (l - kronecker_symbol(D_K, l)) // l
            while n % l == 0:
                n //= l
        l += 1
    if n > 1:
        result = result * (n - kronecker_symbol(D_K, n)) // n
    
    #print(f"h(O_{f}) = {result} (h_OK={h_OK}, D_K={D_K}, f={f})")
    return result // w

def process_t(p, n, t, N, k):

    q = p**n

    N_pts = q + 1 - t

    if N_pts % N != 0:
        return 0,0,0

    j0_SS = p % 3 == 2
    j1728_SS = p % 4 == 3

    D_pi = t**2 - 4*q
    D_K, f_pi = fundamental_discr(D_pi)

    is_quaternion = D_K == 0

    # temp sage double checks
    R = PolynomialRing(ZZ, 'x')
    x = R.gen()
    fx_pi = x**2 - t*x + q

    val = Fraction(0)

    F = GF(q)

    NSS = 0

    # Trace = [-5694.0]

    max_r = 1

    if not is_quaternion:
        K = NumberField(fx_pi, "x")

        h_OK = getClassNumber(D_K, 1)

        for f in ZZ(f_pi).divisors():

            # hk = hk_eval(t, k)

            hk = Hk.eval(fx_pi, t, k)
            # NC = h_order(f, h_OK, D_K)
            # NC = num_C(order_sage=K.order_of_conductor(f), t=t)

            ord_sage = K.order_of_conductor(f)
            h0 = ord_sage.class_number()
            NC = h0

            js = get_j_invariants_from_order(D_K * f**2, f, F)

            NP = num_P(N, f_pi, f, N_pts, q)

            r = 2 if NP == N**2-1 else 1
            if r > max_r:
                max_r = r
                
            if r == 2:
                print(f"{Colors.WARNING}COMPUTING RANK {r} at D_K={D_K}, t={t}, f={f}, NP={NP}, NC={NC}, f_pi={f_pi}{Colors.ENDC}")

            aut_size = 2


            if D_K == -4 and f == 1:
                aut_size = 4
                #print(f"{Colors.BLUE}(4, q-1)={gcd(4, q-1)}{Colors.ENDC}")
            elif D_K == -3 and f == 1:
                aut_size = 6
                #print(f"{Colors.GREEN}(6, q-1)={gcd(6, q-1)}{Colors.ENDC}")

            # if D_K == -4 or D_K == -3:

            '''if len(js) != NC:
                print(
                    f"{Colors.GREEN}p={p}, SPECIAL D_K={D_K}, t={t}, f={f}, f_pi={f_pi}, f_pi mod N={f_pi % N}, NC={NC}, NC_hcp ={len(js)}, NP={NP}{Colors.ENDC}"
                )'''

            if t % p == 0:
                # SUPERSINGULAR BUT QUARATIC ORDERS
                if f % p == 0:
                    print(f"{Colors.WARNING}Skipping trace t={t} and conductor f={f} for F_{q} since they are divisible by char {p}{Colors.ENDC}")
                    continue
                '''if t == 0 and not j1728_SS:
                    continue
                if abs(t) == math.isqrt(q) and not j0_SS:
                    continue'''
                # for SS curves we need to weight the class number count by 2 if inert and 1 if split, since we have 2 curves in the split case and 1 in the inert case, but the class number count does not see this splitting/inertness behavior and just counts the number of ideals in the order
                # if f > 1:
                # print(f"{Colors.FAIL}SKIPPING SS t={t}, f={f}, NC={NC}{Colors.ENDC}")
                #   continue

                # we force back, the class nr is already weighted by the extra auts
                #aut_size = 2

                # print(f"p={p},---------------WE HAVE SS t={t}, f={f}, NC={NC}, NP={NP}, hk={hk}")
                if f == 1:
                    inert_factor = 2 if kronecker(D_K, p) == -1 else 1
                    NC *= inert_factor
                else:
                    inert_factor = 1
                # print(f"inert factoring at t={t}: D_K={D_K}, p={p}, kronecker={kronecker_symbol(D_K, p)}, inert_factor={inert_factor}")
                # NC *= inert_factor

                j_spec = f""
                if j0_SS:
                    j_spec += "j0_SS "
                if j1728_SS:
                    j_spec += f"j1728_SS = {1728 % p}"

                print(
                    f"(6, q-1)={gcd(6, q-1)}, (4, q-1)={gcd(4, q-1)}, h_OK ={K.order_of_conductor(1).class_number()}, D_K={D_K}, -4p={-4*p}"
                )

                NSS += NC

                print(
                    f"{Colors.FAIL}p={p}, f={f}, {j_spec}, Special SS t={t}, h0={h0}, D_K={D_K}, NC={NC}, NSS={NC}, js ={js}, aut_size={aut_size}, inert_factor={inert_factor}{Colors.ENDC}"
                )

            # TODO: if D_K = 0, then we have full queternion and aut size 4 or 6 but solution is to weight the class nr count
            val += Fraction(int(hk * NC * NP), aut_size)
    else:
        # hk = hk_eval(t, k)
        hk = Hk.eval(fx_pi, t, k)
        NC = quaternion_class_number(p)
        NP = num_P(N, f_pi, -1, N_pts, q) # set to -1 so we do NOT look at conductors, only weil pairing and valuation of N_pts
        max_r = 2 if NP == N**2-1 else 1
        val = Fraction(Fraction(int(hk * NP)) * NC, 2) # we have already weighted the class nr count
        NSS += NC

    # print(f"t={t}: hk={hk}, NC={NC}, NP={NP}, contribution to count: {val}")
    # print(f"p={p}, t={t}, D_K={D_K}, adding val={val}")
    return val, NSS, max_r


from math import gcd, lcm


def cusp_term(p, n, N, k):
    total = 0
    q = p**n
    
    for c in divisors(N):
        g = gcd(c, N // c)
        if (q - 1) % g != 0:
            continue
        L = lcm(c, N // c)  # = N // g
        if q % L == 1:
            sgn = 1
        elif q % L == L - 1:
            sgn = -1
        else:
            sgn = 0
        total += euler_phi(c) * euler_phi(N // c) * (1 + sgn**k)
    return total // 4


def run(p, n, N, k, compare=False):

    '''rndm_prime = random_prime(10**4)
    if p == -1:
        p = rndm_prime
        n = 1
        print(f"Random prime : {rndm_prime}")
        q = p**n
    else:'''
    q = p**n

    # SS_poly = supersingular_j_polynomial(p)
    # SS_poly_Fq = SS_poly.change_ring(GF(q))

    # ss_js = []
    # for r, m in SS_poly_Fq.roots(multiplicities=True):
    #    for _ in range(m):
    #        ss_js.append(r)

    # print(f"Supersingular j-invariants (except special js) in F_{q}: {ss_js}")

    HB = math.isqrt(4 * q)
    T = 0

    j0_SS = p % 3 == 2
    j1728_SS = p % 4 == 3

    q_squarefree = n % 2 == 1

    '''print(
        f"{Colors.HEADER}Processing F_{q}, j0_SS={j0_SS}, j1728_SS={j1728_SS}, 1728 mod q = {1728 % q}{Colors.ENDC}, q mod N = {q % N}, p mod N = {p % N}, legendre={kronecker(q, N)}{Colors.ENDC}"
    )'''

    NSS = 0
    max_r = 0

    for t in range(-HB, HB+1):
        # we skip all SS here and add them back below
        if t % p == 0:
            continue
        T_incr, NSS_incr, r = process_t(p, n, t, N, k)
        NSS += NSS_incr
        T += T_incr
        if r > max_r:
            max_r = r

    # ONLY t = 0 occurs as SS trace, all such curves will have |Aut| = 2, ie D_K notin (-3, -4)
    if q_squarefree:
        T_incr, NSS_incr, r = process_t(p, n, 0, N, k)
        NSS += NSS_incr
        T += T_incr
        if r > max_r:
            max_r = r
    else:
        # j0 SS has t in pm sqrt(q)
        if j0_SS:
            T_incr, NSS_incr, r = process_t(p, n, HB // 2, N, k)
            NSS += NSS_incr
            T += T_incr
            if r > max_r:
                max_r = r
            T_incr, NSS_incr, r = process_t(p, n, -HB // 2, N, k)
            NSS += NSS_incr
            T += T_incr
            if r > max_r:
                max_r = r

        if j1728_SS:
            T_incr, NSS_incr, r = process_t(p, n, 0, N, k)
            NSS += NSS_incr
            T += T_incr
            if r > max_r:
                max_r = r

        # Quaternion algebras
        T_incr, NSS_incr, r = process_t(p, n, HB, N, k)
        NSS += NSS_incr
        T += T_incr
        if r > max_r:
            max_r = r
        T_incr, NSS_incr, r = process_t(p, n, -HB, N, k)
        NSS += NSS_incr
        T += T_incr
        if r > max_r:
            max_r = r

    sum_value = int(T)

    # Final count for F_841: -2176, cusp_gamma1=6, FINAL_VALUE=-2182

    cusp_dim = 0

    if N == 17 and p == 11 and n == 4 and k == 1:
        sage_T = -14200
    elif N == 17 and p == 11 and n == 3 and k == 1:
        sage_T = 1936
    # elif N == 11 and p == 19 and n == 3 and k == 1:
    #    sage_T = -5275
    # elif N == 13 and p == 29 and n == 3 and k == 1:
    #    sage_T = -27852
    else:
        cform = CuspForms(Gamma1(N), k + 2)
        sage_T = cform.hecke_operator(q).trace() if compare else 0
        cusp_dim = cform.dimension()
        # print(cform)

    # Final count for F_24389: -57072, sage trace: -27852, difference: -29220
    # Final count for F_24389: -58116, cusp_gamma1=12, FINAL_VALUE=-58128

    diff = -sum_value - sage_T

    phi_N = euler_phi(N)

    sgn = 1 if q % N == 1 else (-1 if q % N == N - 1 else 0)
    sgn_k = sgn**(k+2) #avoid k == 0 sgn**k = 1 issue
    # W = sum(euler_phi(c) * euler_phi(N // c) for c in divisors(N))
    # comp_diff2 = (W / 4) * (1 + sgn_k)
    # comp_diff = (phi_N / 2) * (1 + sgn_k)
    # print(f"sgn_k = {sgn_k}, sgn={sgn}")
    cusp_gamma1_1 = cusp_term(p, n, N, k)
    # cusp_gamma1_1 = (phi_N / 2) * (1 + sgn_k)

    # CUSP EST NEEDS TO CHANGE FOR HIGHER k
    # for n = 2, the below works for p

    c = cusp_dim
    if k % 2 == 0 or sgn == 0:
        c = gcd(cusp_dim, phi_N)
        # GCD IS 0!!!!!!!!!!!!!!!!!!!, MAKE SURE NOT CUSP EST IS 0

    # check q if k odd, check p if k even??

    recursion_q = q#**(n-1)
    # CHECK -3/p anbd -4/p, there are somee xtra cusps then? if SS?

    sgn = 1 if recursion_q % N == 1 else (-1 if recursion_q % N == N - 1 else 0)
    # RESULTS for n=2, N=13, k=3, phi_N=12:

    sgn_k = sgn ** (k + 2)  # avoid k == 0 sgn**k = 1 issue

    # sometimes gcd sometimes not
    cusp_gamma1_2 = -1 * sgn_k * c * (p) ** (k + 1)
    # k = odd works
    print(p, sgn_k, sgn)
    # cusp_gamma1 = (phi_N / 2) * (1 + sgn_k)
    # Processing F_2113, j0_SS=False, j1728_SS=False, 1728 mod q = 1728, q mod N = 13, p mod N = 13
    # comp_diff = phi_N
    # comp_diff = int(phi_N / 2) * ((phi_N) + sgn**k)

    '''comp_str = (
        f"diff={diff}, sage={sage_T}, comp_diff={comp_diff}, comp_diff2={comp_diff2}, cusp_gamma1={cusp_gamma1}"
        if compare
        else f"cusp_gamma1={cusp_gamma1}"
    )
    print(
        f"Final count for F_{q}: {-sum_value}, {comp_str}, FINAL_VALUE={(-sum_value-cusp_gamma1)}, NSS={NSS}"
    )'''

    '''if compare:
        phi_N = euler_phi(N)
        # floor_diff_q = math.floor((diff) / q)
        # ceil_diff_q = math.ceil((diff) / q)
        print(f"phi(N) = {phi_N}")
        if diff != 0:
            print(f"DIFF = {ZZ(diff).factor()}")
        left = diff - comp_diff
        print(f"left = (phi_N / 2) * (1 + sgn_k)={comp_diff}")
        # cusp_estimate = cusp_est(p, n, N, k)
        if left != 0:
            print(f"diff - (phi_N / 2) * (1 + sgn_k) = {ZZ(int(left)).factor()}")'''
    # print(f"LEFT = {ZZ(left / n).factor()}")

    # print(f"T1 : = {ZZ((phi_N / 2) * p**(n-1)).factor()}")
    # print(f"T2 : = {ZZ((phi_N) * p**(n-1)).factor()}")
    # print(float((left) / p), (p - 1) * (phi_N))

    # works for k = 1, N = 11, p = 19, n = 3, q mod N = 6, p mod N = 8
    # print((left) - (p*(1 + (p - 1) * (phi_N))))

    # works for k = 1, N = 13, p = 19, n = 3, q mod N = 8, p mod N = 6
    # print(left - 4*(p-1)*p)

    # works for k = 1, N = 17, p = 19, n = 3, q mod N = 8, p mod N = 2
    # print(left - 4 * (p - 1) * p)

    # print(ZZ(int(cusp_eisenstein_correction(q, N, k+2))).factor())

    # print(A1_identity(q, N, k+2))
    # print(A3_cusp(q, N, k+2))
    # print(A4_eisenstein(q, N, k+2))

    # print(f"phi_N + ceil_diff_q * q={phi_N + ceil_diff_q * q}")
    # print(f"phi_N // 2 + ceil_diff_q * q={phi_N // 2 + ceil_diff_q * q}")

    # print(f"phi_N - ceil_diff_q * q={phi_N - ceil_diff_q * q}")
    # print(f"phi_N // 2 - ceil_diff_q * q={phi_N // 2 - ceil_diff_q * q}")

    # works for non q equiv 1

    return HeckeResult(
        p=p,
        n=n,
        q_equiv_N=-1 if q % N == N - 1 else int(q % N),
        p_equiv_N=-1 if p % N == N - 1 else int(p % N),
        diff=int(diff),
        phi_N=int(phi_N),
        error_1=int(cusp_gamma1_1),
        error_2=int(cusp_gamma1_2),
        legendre_N_p=kronecker(N, p),
        legendre_p_N=kronecker(p, N),
        legendre_N_q=kronecker(N, q),
        legendre_q_N=kronecker(q, N),
        max_r=max_r,
        NSS=NSS,
        j0_SS=j0_SS,
        j1728_SS=j1728_SS,
        cusp_dim=cform.dimension() if cform else None,
    )


if __name__ == "__main__":
    args = parse_args()

    results = []
    if args.p == -1:
        for _p in prime_range(5, args.pmax):
            # print(f"\n{'='*40}\nRunning for p={_p}, n={args.n}, N={args.N}, k={args.k}\n{'='*40}")
            result = run(_p, args.n, args.N, args.k, args.compare)
            results.append(result)
    else:
        run(args.p, args.n, args.N, args.k, args.compare)

    print(f"RESULTS for n={args.n}, N={args.N}, k={args.k}, phi_N={euler_phi(args.N)}:\n")
    M = ModularForms(Gamma1(args.N), args.k + 2)
    dim = M.dimension()
    dim_eis = M.eisenstein_subspace().dimension()
    dim_cusp = M.cuspidal_subspace().dimension()
    print(f"Dimension of M_{args.k+2}(Gamma1({args.N})): {dim} (Eisenstein: {dim_eis}, Cusp: {dim_cusp})\n")

    if results:

        header = f" {'p':>6} {'n':>6} {'q≡N':>6} {'p≡N':>6} {'max_r':>6} {'NSS':>6} {'j0_SS':>6} {'j1728_SS':>6} {'q':>10} {'diff':>10} {'q-diff':>10} {'cusp_dim':>10}  {'cusp_est':>10}  {'adjusted_val_1':>20} {'adjusted_val_2':>20}"
        print(header) 
        print("-" * len(header))
        for r in results:
            clr = Colors.BOLD
            if r.p_equiv_N == 1:
                clr = Colors.GREEN
            elif r.p_equiv_N == -1:
                clr = Colors.WARNING
            elif r.p_equiv_N == 0:
                clr = Colors.FAIL

            adjusted_val_1 = r.diff - r.error_1
            adjusted_val_2 = r.diff - r.error_1 - r.error_2

            print(
                f" {clr}{r.p:>6} {r.n:>6} {r.q_equiv_N:>6} {r.p_equiv_N:>6} {r.max_r:>6} {float(r.NSS):>6.1f} {r.j0_SS:>6} {r.j1728_SS:>6} {(r.p**r.n):>10} {r.diff:>10} {(r.p**r.n - abs(r.diff)):>10} {r.cusp_dim:>10}  {r.error_1:>10} {str(ZZ(int(adjusted_val_1)).factor() if adjusted_val_1 != 0 else 0):>20} {str(ZZ(int(adjusted_val_2)).factor() if adjusted_val_2 != 0 else 0):>20} {Colors.ENDC}"
            )
