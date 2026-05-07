# from sage.all import *
from sage.all import *
import random
from utils.common import Colors
import argparse
from typing import Optional, List, Dict, Tuple, Set, Any
# from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
# from sage.all import CuspForms, Gamma1

from fractions import Fraction
from dataclasses import dataclass, field

from sage.libs.pari import pari
from sympy import primerange, factorint
from tqdm import tqdm

import math
import time

PRECOMPUTE_INVARIANTS = False
# CACHED_DK = {}  # D_K -> first t that produced it
LEVEL_FACTORS = []
AUT_SIZE = {-4: 4, -3: 6}
FILTER_Q_ON_LEVEL = True

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
        "--pmax", type=int, default=100, help="Upper prime bound"
    )

    p.add_argument(
        "--pmin", type=int, default=5, help="Upper prime bound"
    )

    p.add_argument(
        "--compare",
        action="store_true",
        default=False,
        help="Use Sage's built-in Hecke operator for comparison (only for small p and N)",
    )
    
    p.add_argument(
        "--random",
        action="store_true",
        default=False,
        help="Render one random prime in primerange ( n = 1 )",
    )

    p.add_argument(
        "--filter",
        action="store_true",
        default=False,
        help="When enumerating traces in HB, only process those that can have level structure",
    )
        
    p.add_argument(
        "--plist",
        type=int,
        nargs="*",
        default=None,
        help="List of specific primes to process (for debugging)",
    )
    return p.parse_args()


# GEENRIC HELPERS
def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result

def legendre(a, p):
    if p == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1

def valuation(n, l):
    if n == 0:
        return float("inf")
    v = 0
    while n % l == 0:
        n //= l
        v += 1
    return v

def fmt_factored(n):
    if n == 0:
        return "0"
    
    if n > 10**6:
        return str("")
    factors = factorize(abs(n))
    if not factors:
        return str(n)
    s = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
    return f"-{s}" if n < 0 else s

def factorize(n):
    #return list(factorint(n).items())  # → [(2, 2), (3, 1)]
    """Return list of (prime, exponent) pairs for n > 1."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

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

def quaternion_class_number(p, rescale_weights=False):

    """print(f"p={p}, chi3={chi3}, chi4={chi4}")
    H_p = (p + 6 - 4 * chi3 - 3 * chi4) / 12
    print((p + 6 - 4 * chi3 - 3 * chi4))
    return Fraction(int(p + 6 - 4 * chi3 - 3 * chi4), 12)
    n_j0 = (1 - chi3) // 2
    n_j1728 = (1 - chi4) // 2
    n_generic = H_p - n_j0 - n_j1728"""

    # just return N(t)
    if not rescale_weights:
        chi3 = legendre(-3, p)
        chi4 = legendre(-4, p)
        return Fraction(int(p + 6 - 4 * chi3 - 3 * chi4), 12)

    # we weigh the j=0 curves by 1/3 and the j=1728 curves by 1/2, then we apply 1/2 for all auts, this end up in this simplified form
    return Fraction(p - 1, 12)
    # return Fraction(int(n_generic)) + Fraction(int(n_j0), 3) + Fraction(int(n_j1728), 2)

def cusp_term(p, n, N, k):
    q = p**n

    split = 0
    non_split = 0

    for d in divisors(N):
        Nd = N // d
        # split cusps: need (N/d) | (q-1)
        if (q - 1) % Nd == 0:
            split += euler_phi(d) * euler_phi(Nd) // 2
        # non-split cusps: need d | 2 AND (N/d) | (q+1)
        if d in (1, 2) and (q + 1) % Nd == 0:
            non_split += euler_phi(d) * euler_phi(Nd) // 2

    # Each cusp contributes a_1^k where a_1 = +1 (split) or -1 (non-split)
    return split * 1 + non_split * ((-1) ** k)

# HELPER CLASS TO EVALUATE Hk
class Hk:
    @staticmethod
    def eval(q:int, t:int, k: int) -> int:
        return int(sum(
            math.comb(k - j, j) * (-q) ** j * t ** (k - 2 * j)
            for j in range(k // 2 + 1)
        ))

# HELPER METHODS RELATED TO QUADRATIC FIELDS
@dataclass
class QFOrder:
    f:int
    D:int
    h_ord: int #class nr for this order
    inv: Tuple[int, int] # all curves in this order have same invariants

class JTQuadraticField:
    def __init__(self, D: int) -> None:
        self.D_K = self._fundamental_discr(D)
        '''if self.D_K in CACHED_DK:
            clr = Colors.WARNING# if CACHED_DK[self.D_K][2] == p else Colors.FAIL
            print(
                f"{clr}Warning: D_K={self.D_K} already processed, skipping redundant class number computation{Colors.ENDC}"
            )
            self.h_OK = CACHED_DK[self.D_K]
        else:
            # TODO: replace with db lookup, using pari for now
            self.h_OK = -1 if self.D_K == 0 else int(pari(self.D_K).qfbclassno())
            CACHED_DK[self.D_K] = self.h_OK'''

        # Unit index [O_K^x : O_f^x]
        if self.D_K == -3:
            self.w = 3
        elif self.D_K == -4:
            self.w = 2
        else:
            self.w = 1

        # SAGE DEBUG
        if self.D_K != 0:
            #self.K = QuadraticField(self.D_K)
            #int(self.K.class_number())  #
            self.h_OK = int(pari(self.D_K).qfbclassno())
        else:
            self.h_OK = 0

        self.is_gaussian = self.D_K == -4
        self.is_eisentein = self.D_K == -3

    def _fundamental_discr(self, D:int) -> int:
        """Return the fundamental discriminant Δ """
        if D == 0:
            return 0
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
        # f = math.isqrt(abs(D // delta))
        return int(delta)#, f

    def h(self, f:int, q:int) -> int:
        """
        Cox primes, Thm 7.24
        
        Class number of the order of conductor f in the imaginary quadratic field
        with fundamental discriminant D_K.
        Matches Sage's K.order_of_conductor(f).class_number().
        Parameters:
            f    : conductor (positive integer)
        """
        if f == 1:
            return self.h_OK
        # Product over distinct primes dividing f
        result = Fraction(f * self.h_OK, self.w)

        for p, _ in factorize(f):
            result *= Fraction(p - legendre(self.D_K, p), p)

        if result.denominator != 1:
            print(
                f"{Colors.FAIL}Warning: non-integer class number for D_K={self.D_K}, f={f}, h_OK={self.h_OK}, w={self.w}, intermediate result={result}{Colors.ENDC}"
            )
        '''if self.K is not None:
            order_sage = self.K.order_of_conductor(f)
            h_O = order_sage.class_number()
            js = self.j_invariants(f, GF(q))
            if int(result) != int(h_O):
                print(
                    f"{Colors.FAIL}Discrepancy in class number for D_K={self.D_K}, f={f}: computed {legendre(self.D_K, p)} vs Sage's {h_O}, legendre(self.D_K, p) / p={(p - (legendre(self.D_K, p))) / p}{Colors.ENDC}"
                )'''
        # return int(h_O)
        return int(result)

    def j_invariants(self, f:int, F) -> List:
        j_invs = []
        try:
            H = hilbert_class_polynomial(self.D_K*f**2)
            H_fq = H.change_ring(F)
            for j, m in H_fq.roots(multiplicities=True):
                for _ in range(m):
                    j_invs.append(j)
        except Exception as e:
            print(f"Warning: Could not compute HCP for D={self.D_K*f**2}, probably because sage is not loaded: {e}")
        return j_invs

@dataclass
class LevelStructureResult:
    t: int = 0
    has_full: int = 0
    val: Fraction = field(default_factory=Fraction)
    NSS: Fraction = field(default_factory=Fraction)
    NC: Fraction = field(default_factory=Fraction)
    NP: int = 0
    valid: bool = False

@dataclass
class HeckeResult:
    p: int  # prime p
    n: int  # extension degree n
    q_equiv_N: int  # q % N
    p_equiv_N: int  # p % N
    eis: int  # -T - sage_T
    phi_N: int  # euler_phi(N)
    computed_sum: int  # trace computed by our implementation
    Tr: int  # trace computed by Sage (for comparison)
    error_1: int  # cusp_term(q, N, k)
    error_2: int  # higher n,
    legendre_N_q: int  # legendre symbol (q/N)
    legendre_q_N: int  # legendre symbol (N/q)
    legendre_p_N: int  # legendre symbol (p/N)
    legendre_N_p: int  # legendre symbol (N/p),
    max_r: int  # maximum r value encountered (for debugging)
    NSS: int  # number of supersingular curves encountered (for debugging)
    NC: int
    NP: int
    j0_SS: bool  # whether j0 is supersingular (for debugging)
    j1728_SS: bool  # whether j1728 is supersingular (for debugging)

def num_P(level, f_pi, f, N_pts, q):
    if level == 1:
        return 1
    result = 1
    for l, a in LEVEL_FACTORS:
        h = max(0, valuation(f_pi, l) - valuation(f, l)) if f_pi * f > 0 else None
        v_q1 = valuation(q - 1, l)
        v_N = valuation(N_pts, l)
        e1 = min(h, v_q1, v_N // 2) if h is not None else min(v_q1, v_N // 2)
        e2 = v_N - e1
        s1 = min(a, e1)
        s2 = min(a, e2)
        # exact-order-l^a count in Z/l^s1 x Z/l^s2
        result *= l ** (s1 + s2) - l ** (min(a - 1, s1) + min(a - 1, s2))
    return result

def process_t(p, n, t, N, k):
    q = p**n
    N_pts = q + 1 - t
    # safe to early exit, note: for N not prime, we might still return NP = 0 below, so passing this is NOT equivalent to having at least one point of order N

    if N_pts % N != 0:
        return LevelStructureResult(t=t)

    D_pi = t**2 - 4*q
    qf = JTQuadraticField(D_pi)

    D_K = qf.D_K
    f_pi = math.isqrt(abs(D_pi // D_K)) if D_K != 0 else 0
    curves_contrib = Fraction(0)
    NSS = Fraction(0)
    NC = Fraction(0)
    NP = 0

    if D_K != 0: # Imaginary Quadratic Case
        is_SS = t % p == 0
        for f in divisors(f_pi):
            # NON ALLOWED CONDUCTORS
            if is_SS and f % p == 0:
                continue

            # js = qf.j_invariants(f, GF(q))
            # print(f"t={t}, f={f}, D_K={D_K}, f_pi={f_pi}, N_pts={N_pts}, j_invariants={js}")
            NP_ord = num_P(N, f_pi, f, N_pts, q)
            if NP_ord > 0:
                NC_ord = qf.h(f, q)
                NSS += NC_ord if is_SS else 0

                """
                (2 if is_SS else 1)
                Schoof Thm 4.6
                
                we do NOT weigh the class numbers by mass (ie Huruwitz) we instead compensate for the aut weight here, we also compensat for the extra inert factor f*h(O) as described in Schoof, which is precicely 2 iff (-3/p)=-1 or (-4/p)=-1 which occurs iff j0 or j1728 is SS, hence we have inert factor 2 iff is_gaussian or is_eisen and is_SS
                see: Schoof Thm 4.5 and 4.6
                
                """
                aut_size = (AUT_SIZE[D_K] // (2 if is_SS else 1)) if ((qf.is_gaussian or qf.is_eisentein) and f == 1) else 2

                curves_contrib += (
                    Hk.eval(q, t, k)
                    * NC_ord
                    * NP_ord
                    * Fraction(1, aut_size)
                )
                NP += NP_ord
                NC += NC_ord 

    else: # Quaterion Case
        NP = num_P(N, f_pi, -1, N_pts, q) # set to -1 so we do NOT look at conductors, only weil pairing and valuation of N_pts
        NC = Fraction(p - 1, 12) #quaternion_class_number(p, rescale_weights=True)
        NSS = NC
        # we have already weighted the special js with 1/3 and 1/2 in this count, since it might be a mix of generic and special
        curves_contrib = (
            Hk.eval(q, t, k)
            * NC
            * NP
            * Fraction(1, 2)
        )

    return LevelStructureResult(
        t=t,
        has_full=False,
        val=curves_contrib,
        NSS=NSS,
        NC=NC,
        NP=NP,
        valid=True,
    )
# q=31: diff=4, T=-104, sage_trace=-108, NC=13, NSS=13,  [q≡1 mod ell] [q≡1 mod 5]

def run(p, n, N, k, compare=False, filter_q_level=False):
    q = p**n
    '''SS_poly = supersingular_j_polynomial(p)
    SS_poly_Fq = SS_poly.change_ring(GF(q))

    ss_js = []
    for r, m in SS_poly_Fq.roots(multiplicities=True):
        for _ in range(m):
            ss_js.append(r)

    print(f"Supersingular j-invariants (except special js) in F_{q}: {ss_js}")'''

    SQRT_Q = math.isqrt(q)
    HB = math.isqrt(4*q)
    T = 0

    j0_SS = p % 3 == 2
    j1728_SS = p % 4 == 3

    partial_results = []
    partial_results_2 = []

    # NOTE: for q squarefree, |int(sqrt(4q)=int(HB)| is NOT SS trace, hence we need the full range

    i_min = (q + 1 - HB + N - 1) // N
    i_max = (q + 1 + HB) // N

    #ts_ = []

    if False:
        for i in range(i_min, i_max + 1):
            t = q + 1 - i*N
            # print(f"Processing trace t={t} for F_{q} with Hasse bound {HB}")
            if t % p == 0:
                continue
            # ts_.append(t)
            partial_results.append(process_t(p, n, t, N, k))
    else:
        for t in range(1, HB + 1):
            # print(f"Processing trace t={t} for F_{q} with Hasse bound {HB}")
            if t % p == 0:
                continue
            # N(t) = H(t^2-4q), t^2 < HB, p nmid t, Schoof
            partial_results.append(process_t(p, n, t, N, k))
            partial_results.append(process_t(p, n, -t, N, k))

    #ts_ = [r.t for r in partial_results if r.valid]
    '''
    0 7
    6 7
    25 7
    3 2
    5 2
    9 2
    14 2
    8 -3
    27 -3
    11 -8
    17 -8
    13 -8
    19 -8
    '''

    '''ts_2 = [r.t for r in partial_results if r.valid]
    if len(ts_) != len(ts_2) or set(ts_) != set(ts_2):
        print(
            f"\n{Colors.FAIL}Warning: mismatch in expected vs processed traces, expected {(sorted(ts_))}, got {(sorted(ts_2))}{Colors.ENDC}"
        )
    else:
        print(
            f"\n{Colors.GREEN}COORECT PROCESSED TS: expected vs processed traces, expected {(sorted(ts_))}, got {(sorted(ts_2))}{Colors.ENDC}"
        )'''

    '''for t in range(1, HB + 1):
        if t % p == 0:
            continue
        # N(t) = H(t^2-4q), t^2 < HB, p nmid t, Schoof
        partial_results.append(process_t(p, n, t, N, k))
        partial_results.append(process_t(p, n, -t, N, k))'''

    '''for t in range(-HB, HB+1):
        if t % p == 0:
            continue
        # N(t) = H(t^2-4q), t^2 < HB, p nmid t, Schoof
        partial_results.append(process_t(p, n, t, N, k))'''

    # SUPERSINGULAR TRACES
    if n % 2 == 1: # CASE 1: q squarefree
        # N(0) = H(-4q), ie t = 0, Schoof
        # print(f"N(t), t = 0 | H(-4q) | {pari(4*p).qfbhclassno()}")
        partial_results.append(process_t(p, n, 0, N, k))
    else:
        # Schoof special case for j depending wheter the are SS in p
        # TODO: add p = 2,3 also

        # print(f"N(t), t = 0  | {1-legendre(-4,p)}, j1728_SS={j1728_SS}")
        if j1728_SS:
            partial_results.append(process_t(p, n, 0, N, k))

        # print(f"N(t), |t| = SQRT_Q | {2*(1-legendre(-3, p))}, j0_SS={j0_SS}")
        if j0_SS:
            partial_results.append(process_t(p, n, SQRT_Q, N, k))
            partial_results.append(process_t(p, n, -SQRT_Q, N, k))

        # Quaternion algebras |t| = HB, DK = 0, t^2 = 4q
        # print(f"Quat, |t| = HB | {2*quaternion_class_number(p)}") # remember one of each pm HB
        partial_results.append(process_t(p, n, HB, N, k))
        partial_results.append(process_t(p, n, -HB, N, k))

    # SUM TOTAL FROM PARTIALS
    NSS = int(0)
    for r in partial_results:
        NSS += r.NSS
    NSS = int(NSS)

    S_Gamma_Y = Fraction(0)
    for r in partial_results:
        S_Gamma_Y += r.val
    if S_Gamma_Y.denominator != 1:
        print(f"{Colors.FAIL}Total S_Gamma_Y is not an integer: {S_Gamma_Y}{Colors.ENDC}")

    S_Gamma_Y = int(S_Gamma_Y)

    NC = 0
    for r in partial_results:
        NC += float(r.NC)
    NC = int(NC)
    NP = 0
    for r in partial_results:
        NP += r.NP
    NP = int(NP)

    cform = CuspForms(Gamma1(N), k + 2) if compare else None
    Tr_Tq = cform.hecke_operator(q).trace() if compare else 0

    # From Theorem 2.1
    # Tr(F_q | S([H,k+2])) + e_k = - eis_H,k(q) - SUM a_1(E)^(k-2)/#AutE
    # eis_H,k(q) = - SUM a_1(E)^(k-2)/#AutE - Tr(F_q | S([H,k+2]))
    # eis_H,k(q) = - S_Gamma_Y - Tr_Tq
    # this is the expected eis_H_k, we compute it now to be able to find a formula to compute this, and therefore compute Tr
    
    #S_Gamma_Y = 501120
    
    eis_H_k = - S_Gamma_Y - Tr_Tq

    phi_N = euler_phi(N)
    sgn = 1 if q % N == 1 else (-1 if q % N == N - 1 else 0)
    sgn_k = sgn**(k+2) #avoid k == 0 sgn**k = 1 issue

    # TODO: make it work for N < 5
    cusp_gamma1_1 = cusp_term(p, n, N, k+2)
    # if N == 2:
    # cusp_gamma1_1 = (phi_N // 2) * (1 + sgn_k)
    # if N == 4:
    #    cusp_gamma1_1 = (3 / 2) * (1 + sgn_k)

    
    return HeckeResult(
        p=p,
        n=n,
        q_equiv_N=-1 if q % N == N - 1 and N != 2 else int(q % N),
        p_equiv_N=-1 if p % N == N - 1 and N != 2 else int(p % N),
        eis=int(eis_H_k),
        phi_N=int(phi_N),
        computed_sum=int(S_Gamma_Y),
        Tr=int(Tr_Tq),
        error_1=int(cusp_gamma1_1),
        error_2=int(0),
        legendre_N_p=legendre(N, p),
        legendre_p_N=legendre(p, N),
        legendre_N_q=legendre(N, q),
        legendre_q_N=legendre(q, N),
        max_r=0,
        NSS=NSS,
        NC=NC,
        NP=NP,
        j0_SS=j0_SS,
        j1728_SS=j1728_SS,
    )


def print_result(result: HeckeResult):
    if result:
        header = f" {'p':>6} {'n':>6} {'q≡N':>6} {'p≡N':>6} {'max_r':>6} {'NC':>6} {'NSS':>6} {'j0_SS':>10} {'j1728_SS':>10} {'q':>10} {'Computed':>20} {'Tr(T_q)':>20} {'eis':>10} {'cusp_est':>10} {'Error':>30}"
        print(header)
        print("-" * len(header))
        for r in results:
            adjusted_val_1 = r.eis - r.error_1
            clr = Colors.BOLD
            if adjusted_val_1 == 0:
                clr = Colors.GREEN
            # else:
            #    clr = Colors.WARNING
            if r.p_equiv_N == 0:
                clr = Colors.FAIL
            print(
                f" {clr}{r.p:>6} {r.n:>6} {r.q_equiv_N:>6} {r.p_equiv_N:>6} {r.max_r:>6} {int(r.NC):>6} {float(r.NSS):>6.1f} {r.j0_SS:>10} {r.j1728_SS:>10} {(r.p**r.n):>10} {r.computed_sum:>20} {r.Tr:>20} {r.eis:>10} {r.error_1:>10} {fmt_factored(int(adjusted_val_1)):>30}{Colors.ENDC}"
            )

if __name__ == "__main__":
    args = parse_args()

    '''if args.compare:
        try:
            from sage.all import *
        except ImportError:
            print(f"{Colors.FAIL}Error: Sage is required for --compare option{Colors.ENDC}")
            exit(1)'''
    start_t = time.time()
    
    for l, a in factorize(args.N):
        LEVEL_FACTORS.append((l, a))

    results = []

    # no fixed prime, we loop over a range or select random in range
    if args.p == -1:
        primes = [p for p in primerange(args.pmin, args.pmax+1)]
        if args.random:
            rnd_idx = random.randint(0, len(primes)-1)
            results.append( run(primes[rnd_idx], args.n, args.N, args.k, args.compare, args.filter))
        else:
            for _p in tqdm(primes, desc="processing primes", unit="ic", ncols=80, ascii=True):
                # for _p in primes:
                # print(f"\n{'='*40}\nRunning for p={_p}, n={args.n}, N={args.N}, k={args.k}\n{'='*40}")
                result = run(_p, args.n, args.N, args.k, args.compare, args.filter)
                results.append(result)
    else:
        # set p overrides --random flag
        results.append(run(args.p, args.n, args.N, args.k, args.compare, args.filter))

    end_t = time.time()

    q_info = (
        f"q={args.p**args.n}" if args.p > 0 else f"prange={args.pmin}-{args.pmax}"
    )
    print(
        f"RESULTS for {q_info} N={args.N}, k={args.k}, phi_N={euler_phi(args.N)}, computed in {end_t - start_t:.5f} seconds:\n"
    )

    print_result(results)
    # 295175736
    # T=-295175736
    '''M = ModularForms(Gamma1(args.N), args.k + 2)
    E = M.eisenstein_subspace()
    # EISEN IS ALWAYS #CUSPS
    S = M.cuspidal_subspace()
    dim_M = M.dimension()
    dim_E = E.dimension()
    dim_S = S.dimension()
    S_old = S.old_submodule()
    S_new = S.new_submodule()
    dim_S_old = S_old.dimension()
    dim_S_new = S_new.dimension()
    print(f"Dimension of M_{args.k+2}(Gamma1({args.N})): {dim_M} (Eisenstein: {dim_E}, Cusp: {dim_S}, Old: {dim_S_old}, New: {dim_S_new})\n")'''

    # from brute force
    # p=11, q=1331, q equiv ell = 5, T=195936, sage_T=227874, diff=-31938, NC=221, NSS=221, full_r=False
    # Trace = [195936.0]

    # q=6859: diff=1376711046, T=3958253518, sage_trace=2581542472, NC=5146, NSS=5146,  [q≡1 mod ell] [q≡1 mod 3]
    # q=13
    # q=13: diff=2, T=640, sage_trace=638, NC=12, NSS=12,  [q≡1 mod ell] [q≡1 mod 3]
    '''
          5      3     -1     -1      0    112    2.0          1          0        125                56032               -37284     -18748          2               -1 * 2 * 3 * 5^5               -1 * 2 * 3 * 5^5 
      7      3      1      1      0    190    0.0          0          1        343             -1952842              1280560     672282          2                  2^3 * 5 * 7^5                  2^3 * 5 * 7^5 
     11      3     -1     -1      0   1113    4.0          1          1       1331            -93092150              2259384   90832766          2            2^2 * 3 * 11^5 * 47            2^2 * 3 * 11^5 * 47 
     13      3      1      1      0   1224    0.0          0          0       2197            450960728           -214075796 -236884932          2        -1 * 2 * 11 * 13^5 * 29        -1 * 2 * 11 * 13^5 * 29 
     17      3     -1     -1      0   4148    4.0          1          0       4913           3070812652          -1818498780 -1252313872          2      -1 * 2 * 3^2 * 7^2 * 17^5      -1 * 2 * 3^2 * 7^2 * 17^5 
     19      3      1      1      0   3750    0.0          0          1       6859          -3958253518           2581542472 1376711046          2               2^2 * 19^5 * 139               2^2 * 19^5 * 139 
     '''
