# from sage.all import *
from sage.all import *
import random
from utils.common import Colors
import argparse
from typing import Optional, List, Dict, Tuple, Set, Any
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
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
A1q: dict = {}

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
# ---------------------------------------------------------------------------
# count_A1q(p, r)
# For q = p^r, iterates over all possible Frobenius traces a with |a| <= 2*sqrt(q).
# For each isogeny class (determined by trace a), computes the weighted count of
# elliptic curves over F_q in that class, where each curve is weighted by 1/|Aut_k(E)|.
# Stores the result table in A1q[q] = {a: weighted_count}.
# The grand total sum_{a} (weighted count) should equal q.
#
# The 9 cases arise from Honda-Tate theory for abelian varieties over finite fields:
#   Case 1: gcd(a, p) = 1         — ordinary curves, generic case
#   Case 2: a=0, r odd            — supersingular, all curves have the same endomorphism algebra
#   Case 3: a^2=2q, p=2, r odd   — special supersingular for p=2
#   Case 4: a^2=3q, p=3, r odd   — special supersingular for p=3
#   Case 5: a=2√q, p=2, r even   — purely inseparable Frobenius
#   Case 6: a=2√q, p=3, r even   — purely inseparable Frobenius
#   Case 7: a=2√q, general r even — Frobenius = scalar; requires correction for j=0 and j=1728
#   Case 8: a^2=q, r even         — Frobenius has order 2 over F_{p^{r/2}}
#   Case 9: a=0, r even           — purely imaginary Frobenius
# ---------------------------------------------------------------------------
def count_A1q(p: int, r: int):
    q = p**r
    if q >= 10**7 / 4:
        print("q too large")
        return -1

    Res: dict = {}
    amax = floor(2 * sqrt(q))
    total = Fraction(0)

    for a in range(0, amax + 1):

        if a % p != 0:
            # Case 1: ordinary — use Hurwitz-Kronecker class number for discriminant a^2-4q
            dscr = a**2 - 4 * q
            res = HKclass(dscr)

        elif a == 0 and r % 2 == 1:
            # Case 2: supersingular at a=0, odd extension degree
            dscr = -4 * p
            res = HKclass(dscr)

        elif a**2 == 2 * q and p == 2 and r % 2 == 1:
            # Case 3: p=2, r odd, a=sqrt(2q) — 4 automorphisms
            res = Fraction(1, 4)

        elif a**2 == 3 * q and p == 3 and r % 2 == 1:
            # Case 4: p=3, r odd, a=sqrt(3q) — 6 automorphisms
            res = Fraction(1, 6)

        elif a**2 == 4 * q and r % 2 == 0 and p == 2:
            # Case 5: p=2, r even, a=2sqrt(q)
            res = Fraction(1, 24)

        elif a**2 == 4 * q and r % 2 == 0 and p == 3:
            # Case 6: p=3, r even, a=2sqrt(q)
            res = Fraction(1, 12)

        elif a**2 == 4 * q and r % 2 == 0:
            # Case 7: general p, r even, a=2sqrt(q)
            # Base count from a formula involving p and Legendre symbols at -3 and -4
            res = Fraction(
                p + 6 - 4 * legendre(-3, p) - 3 * legendre(-4, p), 24
            )
            # Corrections for the j=0 (CM by Z[ω], extra 6 auts) and j=1728 (CM by Z[i], extra 4 auts) curves
            if p % 3 != 1:  # p is inert or ramified at 3 → j=0 curve appears here
                res += Fraction(-1, 2) + Fraction(1, 6)
            if p % 4 != 1:  # p is inert or ramified at 2 → j=1728 curve appears here
                res += Fraction(-1, 2) + Fraction(1, 4)

        elif a**2 == q and r % 2 == 0:
            # Case 8: a=sqrt(q), r even — involves only -3 Legendre symbol
            res = Fraction(1 - legendre(-3, p), 6)

        elif a == 0 and r % 2 == 0:
            # Case 9: a=0, r even — involves only -4 Legendre symbol
            res = Fraction(1 - legendre(-4, p), 4)

        else:
            res = Fraction(0)  # empty isogeny class

        Res[a] = res
        Res[-a] = (
            res  # trace a and -a give isomorphic (dual) isogeny classes with same count
        )

        # Count a=0 once, all other ±a pairs count twice
        total += res if a == 0 else 2 * res

    if total != q:
        print(f"\n{Colors.FAIL}Mistake! Total does not equal q.{Colors.ENDC}")
    '''else:
        print(f"\n{Colors.GREEN}Success! Total A1q equals q as expected.{Colors.ENDC}")'''

    A1q[q] = Res

    return total


# ---------------------------------------------------------------------------
# HKclass(dscr)
# Computes the Hurwitz-Kronecker class number H(dscr).
# Uses the convention that weights each class by 1/|Aut|, so H(-3)=1/3, H(-4)=1/2.
# Requires dscr < 0 and dscr ≡ 0 or 1 (mod 4).
# ---------------------------------------------------------------------------

LstClnmb: dict  # e.g. {-3: Fraction(1,1), -4: Fraction(1,1), -7: Fraction(1,1), ...}

def HKclass(dscr) -> Fraction:
    
    
    #if dscr >= 0 or (dscr % 4 != 0 and dscr % 4 != 1):
    #    print(f"{Colors.FAIL}Invalid discriminant {dscr} for HKclass, must be negative and ≡ 0 or 1 mod 4{Colors.ENDC}")
    #    return 0
    
    # Factor |dscr| into a list of (prime, exponent) pairs
    fac_dscr = list(factorint(abs(dscr)).items())  # [(p1, e1), (p2, e2), ...]

    # Write dscr = cnd^2 * fund_dscr where fund_dscr is a fundamental discriminant.
    # For each prime factor p^e of dscr:
    #   - if e is even:  p^(e/2) goes entirely into cnd
    #   - if e is odd:   p^((e-1)/2) goes into cnd, and p goes into fund_dscr
    cnd = 1
    fund_dscr = -1  # will accumulate the odd-exponent primes
    for prime, exp in fac_dscr:
        if exp % 2 == 0:
            cnd *= prime ** (exp // 2)
        else:
            cnd *= prime ** ((exp - 1) // 2)
            fund_dscr *= prime

    # Adjust so that fund_dscr is a true fundamental discriminant (≡ 0 or 1 mod 4).
    # If fund_dscr ≡ 2 or 3 (mod 4), we need to absorb an extra factor of 4.
    if fund_dscr % 4 == 2 or fund_dscr % 4 == 3:
        fund_dscr *= 4
        cnd //= 2  # compensate: (2*cnd')^2 * fund_dscr' = cnd^2 * (4*fund_dscr'), so cnd' = cnd/2

    # H(dscr) = sum_{d | cnd} h(fund_dscr) * d * prod_{p | d} (1 - (fund_dscr/p)/p)
    # where (fund_dscr/p) is the Legendre symbol and h is the class number of fund_dscr.
    res = Fraction(0)

    for d in divisors(cnd):

        clnmb = Fraction(int(pari(fund_dscr).qfbclassno()) * d)
        fac_d = list(factorint(d).items())  # prime factors of this divisor d

        start_idx = 0  # index into fac_d; may skip p=2 (handled specially below)

        if len(fac_d) > 0:
            if fac_d[0][0] == 2:
                # Special Euler factor at p=2 depends on fund_dscr mod 8
                start_idx = 1
                if fund_dscr % 8 == 1:
                    clnmb *= Fraction(1, 2)  # (1 - 1/2)
                elif fund_dscr % 8 == 5:
                    clnmb *= Fraction(3, 2)  # (1 + 1/2)
                # if fund_dscr % 8 == 4 or 0: no factor (ramified at 2, contributes 1)

            # Euler factors at odd primes: multiply by (1 - (fund_dscr/p)/p)
            for prime, _ in fac_d[start_idx:]:
                clnmb *= 1 - Fraction(legendre_symbol(fund_dscr, prime), prime)

        res += clnmb

    # Divide by the number of automorphisms of the fundamental order:
    # |Aut(O_{fund_dscr})| = 6 if D=-3, 4 if D=-4, 2 otherwise.
    if fund_dscr == -3:
        res /= 6
    elif fund_dscr == -4:
        res /= 4
    else:
        res /= 2

    return res


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
    if n > 10**20:
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
            split += euler_phi(d) * euler_phi(Nd)
        # non-split cusps: need d | 2 AND (N/d) | (q+1)
        if d in (1, 2) and (q + 1) % Nd == 0:
            non_split += euler_phi(d) * euler_phi(Nd)

    #print(f"cusp_term: p={p}, n={n}, N={N}, k={k}, split={split}, non_split={non_split}")
    # Each cusp contributes a_1^k where a_1 = +1 (split) or -1 (non-split)
    return (split * 1 + non_split * ((-1) ** k)) // 2


def dim_sk(k: int) -> int:
    if k < 0 or k % 2 == 1:
        return 0
    elif k == 2:
        return -1  # S_2 has dimension 0, but the formula gives -1 (boundary case)
    elif k % 12 == 2:
        return floor(k / 12) - 1
    else:
        return floor(k / 12)

# HELPER CLASS TO EVALUATE Hk
class Hk:
    @staticmethod
    def hk_poly(q: int, a: int, k: int) -> int:
        return sum(
            (-1)**i * q**i * math.comb(k - i, i) * a**(k - 2*i)
            for i in range(k // 2 + 1)
        )

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
    @staticmethod
    def H(D: int) -> int:
        """
        Kronecker class number H(Δ) as defined in Schoof (1987), Prop 2.2.
        Returns the cardinality of SL2(Z)-orbits of positive definite 
        binary quadratic forms of discriminant D.
        """
        if D >= 0:
            return 0

        abs_D = abs(D)
        total_H = 0

        # Iterating through d where d^2 divides D (Schoof Prop 2.2)
        # Schoof counts orbits, not weighted values.
        for d in range(1, int(math.isqrt(abs_D)) + 1):
            if abs_D % (d**2) == 0:
                delta_prime = D // (d**2)

                # The discriminant of the order must be ≡ 0 or 1 (mod 4)
                if delta_prime % 4 == 0 or delta_prime % 4 == 1:
                    # qfbclassno(Δ) returns the primitive class number h(Δ)

                    total_H += int(pari(delta_prime).qfbclassno())

        return total_H

    def D0(D: int) -> int:
        """Return the fundamental discriminant Δ"""
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
        return int(delta)  # , f

    def __init__(self, D: int) -> None:
        self.D_K = JTQuadraticField.D0(D)
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
            # self.K = QuadraticField(self.D_K)
            # int(self.K.class_number())  #
            self.h_OK = int(pari(self.D_K).qfbclassno())
        else:
            self.h_OK = 0

        self.is_gaussian = self.D_K == -4
        self.is_eisentein = self.D_K == -3

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
    N_gauss: int = 0
    N_eisen: int = 0
    valid: bool = False

@dataclass
class HeckeResult:
    p: int  # prime p
    n: int  # extension degree n
    q_equiv_N: int  # q % N
    N_equiv_p: int  # p % N
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
    N_gauss: int
    N_eisen: int
    j0_SS: bool  # whether j0 is supersingular (for debugging)
    j1728_SS: bool  # whether j1728 is supersingular (for debugging)
    DIM:int

def num_P(level, f_pi, f, N_pts, q):
    if level == 1:
        return 1
    result = 1
    for l, a in LEVEL_FACTORS:
        h = max(0, valuation(f_pi, l) - valuation(f, l)) if f_pi * f > 0 else None
        v_q1 = valuation(q - 1, l)
        v_N = valuation(N_pts, l)
        # h is really the cap, this decides the SPREAD, ie the invariant ranges from s // 2 to s cyclic, but we also only care about vN // 2. the q-1 is the EXTRA contition for optimixzation
        e1 = min(h, v_N // 2)#v_N // 2 #min(h, v_q1, v_N // 2) if h is not None else min(v_q1, v_N // 2)
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

    # if N_pts % N != 0:
    #    return LevelStructureResult(t=t)

    D_pi = t**2 - 4*q

    qf = JTQuadraticField(D_pi)

    D_K = qf.D_K
    f_pi = math.isqrt(abs(D_pi // D_K)) if D_K != 0 else 0

    '''if t == 0:
        print(
            f"{Colors.GREEN}n={n}, SS trace 0 t={t} yields D_K ={ZZ(D_K).factor() if D_K != 0 else 'N/A'}, D_pi={ZZ(D_pi).factor() if D_pi != 0 else 'N/A'}, D_K//D_pi={ZZ(D_pi // D_K).factor() if D_K != 0 else 'N/A'} and f_pi={ZZ(f_pi).factor() if f_pi != 0 else 'N/A'}{Colors.ENDC}"
        )

    if abs(t) == sqrt(q):
        print(
            f"{Colors.BRIGHT_YELLOW}n={n}, SS trace sqrt(q) t={t} yields D_K ={ZZ(D_K).factor() if D_K != 0 else 'N/A'}, D_pi={ZZ(D_pi).factor() if D_pi != 0 else 'N/A'}, D_K//D_pi={ZZ(D_pi // D_K).factor() if D_K != 0 else 'N/A'} and f_pi={ZZ(f_pi).factor() if f_pi != 0 else 'N/A'}{Colors.ENDC}"
        )'''

    curves_contrib = Fraction(0)
    NSS = Fraction(0)
    NC = Fraction(0)
    NP = 0
    N_gauss = 0
    N_eisen = 0

    valid = False

    # this is the Hurwitz-Kronecker class number H(t^2-4q) that counts the number of isogeny classes with trace t, weighted by aut size, see Schoof 1987, Prop 2.2

    j0_SS = p % 3 == 2
    j1728_SS = p % 4 == 3

    if D_K != 0: # Imaginary Quadratic Case

        is_SS = t % p == 0

        # Schoof: if q non square and ss t, then only maximal orders occur, the other conductors will be divisible by p. we can make early exit here
        H_t = HKclass(D_pi)
        H_t_accum = 0

        if is_SS:
            inert_scaling = 2 if legendre(D_K, p) == -1 else 1
            H_t = HKclass(D_K) * inert_scaling
            # HKclass(D_K) if ( n % 2 == 1 and is_SS) else HKclass(D_pi)

        K = QuadraticField(D_K)

        for f in divisors(f_pi):
            # NON ALLOWED CONDUCTORS
            D = D_K * f**2

            if is_SS and f % p == 0:
                # print(f"{Colors.RED}Skipping conductor f={f}, D_K={D_K}, D={D} for supersingular trace t={t} divisible by p={p}{Colors.ENDC}")
                continue

            if is_SS:

                j0_curve = EllipticCurve(GF(q), j=0)
                for et in j0_curve.twists():

                    print(
                        f"twist t = {et.trace_of_frobenius()} has {len(et.automorphisms())} automorphisms"
                    )

                js = qf.j_invariants(f, GF(q))
                print(
                    f"{Colors.RED}p={p}, runing SS t = {t}, f_pi ={f_pi}, D_K={D_K}, f={f}, j0_SS={j0_SS}, j1728_SS={j1728_SS}, 1728={1728 % p}, js={js}{Colors.ENDC}"
                )

            # discr of order

            is_maximal = f == 1
            is_maximal_at_N = f % N != 0

            js = []#qf.j_invariants(f, GF(q)) if (is_SS) else [] 

            # SCHOOF LEMMA 4.8
            '''if is_SS:
                js = qf.j_invariants(f, GF(q)) if (is_SS) else []
                for j in js:
                    sage_curve = EllipticCurve(GF(q), j=j)
                    if sage_curve:
                        twists = sage_curve.twists()
                        for et in twists:
                            # if abs(et.trace_of_frobenius()) == abs(t):
                            # cyclic if q equiv 4 is not 1
                            print(
                                f"\n p={p}, t={t}, D_K={D_K}, twist t={et.trace_of_frobenius()}, invariants={et.abelian_group().invariants()}, q equiv 4={q % 4}"
                            )'''

            NP_per_curve = num_P(N, f_pi, f, N_pts, q)
            if NP_per_curve > 0:

                valid = NP_per_curve > 0 # if we enter here at least one has torsion

                """
                (2 if is_SS else 1)
                Schoof Thm 4.6
                
                we do NOT weigh the class numbers by mass (ie Huruwitz) we instead compensate for the aut weight here, we also compensat for the extra inert factor f*h(O) as described in Schoof, which is precicely 2 iff (-3/p)=-1 or (-4/p)=-1 which occurs iff j0 or j1728 is SS, hence we have inert factor 2 iff is_gaussian or is_eisen and is_SS
                see: Schoof Thm 4.5 and 4.6
                
                """
                aut_size = (AUT_SIZE[D_K] // (2 if is_SS else 1)) if ((qf.is_gaussian or qf.is_eisentein) and f == 1) else 2

                if qf.is_gaussian and f == 1:
                    N_gauss = 1
                elif qf.is_eisentein and f == 1:
                    N_eisen = 1
                    # if p == 19:
                    #    print(f"Debug: Found Eisenstein prime p={p}, t={t}, D_K={D_K}, f_pi={f_pi}, N_pts={N_pts}, NC_ord={NC_ord}, NP_per_curve={NP_per_curve}, aut_size={aut_size}, isSS={is_SS}")

                # Schoof
                # note, this only happens for q square and SS curves?? - DOUBLE CHEKC
                '''inert_scaling = 2 if legendre(D_K*f**2, p) == -1 else 1
                clr = Colors.GREEN if is_SS else Colors.BOLD
                if inert_scaling == 2:
                    clr = Colors.WARNING
                print(
                    f"{clr}p={p}, t={t}, f={f}, D_K={D_K}, f_pi={f_pi}, N_pts={N_pts}, j_invariants={js}, j1728={1728 % p}, is_SS={is_SS}, aut_size={aut_size}, (-3/p)={legendre(-3, p)}, (-4/p)={legendre(-4, p)}, NC_ord={NC_ord}, NP_ord={NP_ord}, inert_scaling={inert_scaling}{Colors.ENDC}"
                )
                
                sage_curve = EllipticCurve(GF(q), j=0) if js else None
                
                if sage_curve:
                    
                    twists = sage_curve.twists()
                    for et in twists:
                        print(f"twist t={et.trace_of_frobenius()}, invariants={et.abelian_group().invariants()}")
                        pts = [p for p in sage_curve.points() if p.order() == N]
                        print(f"Points of order {N} on Sage's curve: {pts}, count={len(pts)}")'''

                # NON weighted actual N(t) count from schoof

                aut_size = AUT_SIZE[D_K]  if ( (qf.is_gaussian or qf.is_eisentein) and is_maximal ) else 2

                # if p > 3, and q even, we can ONLY end up with is_SS and D_K = -3 or -4 in case of special js:
                # ( (qf.is_gaussian or qf.is_eisentein) and is_maximal ) should be inertscaling == 2

                inert_scaling = 2 if legendre(D, p) == -1 else 1
                inert_scaling_2 = (
                    (1 - legendre(D_K, p))
                    if ((qf.is_gaussian or qf.is_eisentein) and is_maximal and is_SS)
                    else 1
                )

                NE_ord_true = qf.h(f, q) * inert_scaling

                clr = (Colors.GREEN if inert_scaling == 1 else Colors.WARNING) if is_SS else Colors.BOLD
                if inert_scaling != inert_scaling_2:
                    clr = Colors.FAIL
                    print(
                        f"{clr}{'GAUSS' if qf.is_gaussian else 'EISEN' if qf.is_eisentein else ''} || p={p}, t={t}, D_K={D_K}, f={f}, D={D}, f_pi={f_pi}, N_pts={N_pts}, is_SS={is_SS}, aut_size={aut_size}, inert_scaling={inert_scaling}, legendre({D}, p)={legendre(D, p)}, inert_scaling_2={inert_scaling_2}, NE_ord_true={NE_ord_true}, js={js}{Colors.ENDC}"
                    )

                # rescale NC by the aut size
                NE_ord_weight = NE_ord_true * Fraction(1, aut_size)
                # first count ACTUAL TOTAL points in this class order ( no aut weights)

                curves_contrib += Hk.hk_poly(q, t, k) * NP_per_curve * NE_ord_weight #later, only save the NP weight and avry hk eval over k

                # for global isogeny class count
                NC += NE_ord_true
                NP += NP_per_curve * NE_ord_true # we accum the ACTUAL nr of points
                NSS += NE_ord_true if is_SS else 0

                H_t_accum += NE_ord_weight

                # now HKclass should equal NC_ord =

                js = []#qf.j_invariants(f, GF(q))
                # order_sage = K.order_of_conductor(f)
                # h_O = order_sage.class_number()

                '''if is_SS:
                    js = qf.j_invariants(f, GF(q))
                    print(
                        f"{Colors.BOLD}p={p}, D_K={D_K}, D={D}, f={f}, t={t}, for supersingular trace t={t}, q_square={n % 2 == 0}, NE_ord_true={NE_ord_true}, NE_ord_weight={NE_ord_weight}, js={len(js)}, inert_scaling={inert_scaling}{Colors.ENDC}"
                    )'''

                '''print(
                    f"SS={is_SS}, p={p}, D_K={D_K}, D={D}, f={f}, t={t}, Huruwitz weighted={HKclass(D)}, NC={NC}, NE_ord_weight={NE_ord_weight}, NE_ord_true={NE_ord_true }, aut_size={aut_size}, js={len(js)}, h_O={h_O}, inert_scaling={inert_scaling}"
                )'''

        '''if H_t != H_t_accum:
            print(
                f"{Colors.FAIL}Discrepancy in class number for t={t}, D_K={D_K}: expected={H_t}, got NC={NC}, H(-4p)={HKclass(-4*p)}, H_t_accum={H_t_accum}{Colors.ENDC}"
            )'''
        """else:
            print(f"{Colors.GREEN} Finished processing I({t}) : Expecting H({D_pi})={H_t}, got NC={NC}{Colors.ENDC}")"""

    else: # Quaterion Case
        # TODO: we know the invariants of the curve here.
        # the other SS always has cyclic
        # set to -1 so we do NOT look at conductors, only weil pairing and valuation of N_pts

        NP_per_curve = num_P(N, f_pi, -1, N_pts, q)

        '''j0_curve = EllipticCurve(GF(q), j=0)
        # _t = j0_curve.trace_of_frobenius()
        # if _t != t:
        for et in j0_curve.twists():
            _t = et.trace_of_frobenius()
            if _t == t:
                print(
                    f"{Colors.WARNING}t={t} Found twist of j=0 with trace t={_t}, invariants={et.abelian_group().invariants()}, NP_per_curve={NP_per_curve}{Colors.ENDC}"
                )
                pts = [pt for pt in et.points() if pt.order() == N]
                print(f"Points of order {N} count={len(pts)}")
                curves_contrib += Hk.eval(q, t, k) * len(pts) * Fraction(1, 6) # we weigh the j=0 curves by 1/3 since they have aut size 6 instead of 2
                NSS += 1 if len(pts) > 0 else 0
                break

        j1728_curve = EllipticCurve(GF(q), j=1728)
        # _t = j0_curve.trace_of_frobenius()
        # if _t != t:
        for et in j1728_curve.twists():
            _t = et.trace_of_frobenius()
            if _t == t:
                print(
                    f"{Colors.WARNING}t={t} Found twist of j=1728 with trace t={_t}, invariants={et.abelian_group().invariants()}, NP_per_curve={NP_per_curve}{Colors.ENDC}"
                )

                pts = [pt for pt in et.points() if pt.order() == N]
                print(f"Points of order {N} count={len(pts)}")
                curves_contrib += Hk.eval(q, t, k) * len(pts) * Fraction(1, 4)
                NSS += 1 if len(pts) > 0 else 0
                break'''

        # aut_size = 2 if p > 3 else (24 if p == 2 else 12)
        NE_ord_weight = Fraction(p-1, 24) # for p = 2, we have 1/24, for p = 3, we have 1/12, for p > 3, we have mixed aut weights

        NE_ord_true = quaternion_class_number(p, rescale_weights=False) # true number, integer, no weights

        NC += NE_ord_true
        NP += NP_per_curve * NE_ord_true # we accum the ACTUAL nr of points
        NSS += NE_ord_true if NP_per_curve > 0 else 0

        valid = NP_per_curve > 0
        print(
            f"{Colors.HEADER}Quaternion case: p={p}, t={t}, D_K={D_K}, f_pi={f_pi}, Fraction(p-1, 12)={Fraction(p-1, 12)}, NE_ord_weight={NE_ord_weight}, NSS={NSS}, j0SS={p % 3 == 2}, j1728SS={p % 4 == 3}{Colors.ENDC}"
        )
        # for p = 2,3 there is only one curve, but we have to adjust aut size, ie applies to this single curve, for p > 3 we might have different aut sizes, so we apply an offset weighting to the sum in quaternon H
        aut_size = 2 if p > 3 else (24 if p == 2 else 12) # for p=2 we have 24 auts, for p=3 we have 12 auts, for p>3 we have 2 auts, this is the reason for the rescaling of the class number above, since we are not weighing by mass, we need to compensate here by dividing by the aut size, which end up in this simplified form
        # we have already weighted the special js with 1/3 and 1/2 in this count, since it might be a mix of generic and special
        curves_contrib = Hk.eval(q, t, k) * NE_ord_weight * NP_per_curve
        """print(
            f"{Colors.HEADER}Quaternion case: p={p}, t={t}, D_K={D_K}, f_pi={f_pi}, N_pts={N_pts}, NC={NC}, NP={NP_per_curve}, aut_size={Fraction(1, aut_size)}, Hk.eval(q, t, k)={Hk.eval(q, t, k)}{Colors.ENDC}"
        )
        sage_curve = EllipticCurve(GF(q), j=1728)
        if sage_curve:
            twists = sage_curve.twists()
            for et in twists:
                # if abs(et.trace_of_frobenius()) == abs(t):
                print(f"twist t={et.trace_of_frobenius()}, invariants={et.abelian_group().invariants()}")"""

    '''if NSS > 0 and curves_contrib != 0:
        print(
            f"{Colors.GREEN}\n p={p}, SS CONTRIB = {ZZ(int(curves_contrib)).factor()}, frac={curves_contrib}, Hk.eval(q, t, k)={ZZ(int(Hk.eval(q, t, k))).factor() if Hk.eval(q, t, k) != 0 else 0}{Colors.ENDC}"
        )'''
    return LevelStructureResult(
        t=t,
        has_full=False,
        val=curves_contrib,
        NSS=NSS,
        NC=NC,
        NP=NP,
        valid=valid,
        N_gauss=int(N_gauss),
        N_eisen=int(N_eisen),
    )
# q=31: diff=4, T=-104, sage_trace=-108, NC=13, NSS=13,  [q≡1 mod ell] [q≡1 mod 5]

def run(p, n, N, k, compare=False, filter_q_level=False):
    q = p**n

    SS_poly = supersingular_j_polynomial(p)
    SS_poly_Fq = SS_poly.change_ring(GF(q))

    ss_js = []
    for r in SS_poly_Fq.roots(multiplicities=False):
        ss_js.append(r)

    #
    if len(ss_js) > 0:
        print(f"\n{Colors.FAIL}FOUND GENERIC supersingular curve with js={ss_js}{Colors.ENDC}")

    SQRT_Q = math.isqrt(q)
    HB = math.isqrt(4*q)
    T = 0

    j0_SS = p % 3 == 2
    j1728_SS = p % 4 == 3

    partial_results = []

    # NOTE: for q squarefree, |int(sqrt(4q)=int(HB)| is NOT SS trace, hence we need the full range

    i_min = (q + 1 - HB + N - 1) // N
    i_max = (q + 1 + HB) // N

    '''print(
        f"\n{Colors.CYAN}Running for p={p}, n={n}, q={q}, HB={HB}, j0_SS={j0_SS}, j1728_SS={j1728_SS}, expected_NC={count_A1q(p, n)}-------------------{Colors.ENDC}\n"
    )'''
    if filter_q_level:
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

    # SUPERSINGULAR TRACES
    if n % 2 == 1: # CASE 1: q squarefree
        # N(0) = H(-4q), ie t = 0, Schoof
        #print(f"{Colors.HEADER}N(0) = H(-4p) = {JTQuadraticField.H(-4*p)}{Colors.ENDC}")
        partial_results.append(process_t(p, n, 0, N, k))

        # Schoof Thm 4.2
        if p == 2 or p == 3:
            #print(f"{Colors.HEADER}N(pm sqrt(p^(n+1))={p**((n+1)//2)}) = 1 (p = 2,3) {Colors.ENDC}")
            t = p**((n+1) // 2) # n is odd, this is integer pm sqrt(2q) and pm sqrt(3q) are SS traces
            partial_results.append(process_t(p, n, t, N, k))
            partial_results.append(process_t(p, n, -t, N, k))

    else:
        # Schoof special case for j depending wheter the are SS in p
        # TODO: add p = 2,3 also
        N_gauss = 1 - legendre(-4, p)
        N_eisen = 1 - legendre(-3, p)

        #print(f"N(0) = {1-legendre(-4,p)}")

        if N_gauss > 0:
            partial_results.append(process_t(p, n, 0, N, k))

        #print(f"N(pm SQRT_Q={SQRT_Q}) = {(1-legendre(-3, p))}")

        if N_eisen > 0:
            partial_results.append(process_t(p, n, SQRT_Q, N, k))
            partial_results.append(process_t(p, n, -SQRT_Q, N, k))

        # THIS ALWAYS HAPPENS, SINCE q EVEN
        # Quaternion algebras |t| = HB, DK = 0, t^2 = 4q
        # print(f"Quat, |t| = HB | {2*quaternion_class_number(p)}") # remember one of each pm HB
        partial_results.append(process_t(p, n, HB, N, k))
        partial_results.append(process_t(p, n, -HB, N, k))

    # SUM TOTAL FROM PARTIALS
    NC = sum(r.NC for r in partial_results)

    NSS = sum(r.NSS for r in partial_results)
    NP = sum(r.NP for r in partial_results)
    N_eis = sum(r.N_eisen for r in partial_results)
    N_gauss = sum(r.N_gauss for r in partial_results)

    NUM_TRACES = sum(1 for r in partial_results if r.valid)

    '''for r in partial_results:
        if r.valid:
            print(
                f"{Colors.GREEN}t={r.t} had {r.NP} {N}-torsion points total, across {r.NC} curves, contributing={r.val} to accumed sum {Colors.ENDC}"
            )
        else:
            print(
                f"{Colors.FAIL}t={r.t} had zero {N}-torsion points total, across {r.NC} curves, contributing={r.val} to accumed sum {Colors.ENDC}"
            )

    S_Gamma_Y = sum(r.val for r in partial_results)'''

    S_Gamma_Y = sum(r.val for r in partial_results)
    # print(f"\n{Colors.CYAN}--------------q={q} yieldsd Total {S_Gamma_Y} sum, NUM_TRACES={NUM_TRACES}{Colors.ENDC}\n")

    if S_Gamma_Y.denominator != 1:
        print(f"{Colors.FAIL}Total S_Gamma_Y is not an integer: {S_Gamma_Y}{Colors.ENDC}")

    S_Gamma_Y = int(S_Gamma_Y)

    cform = CuspForms(Gamma1(N), k + 2) if compare else None
    Tr_Tq = cform.hecke_operator(q).trace() if compare else 0

    # From Theorem 2.1
    # Tr(F_q | S([H,k+2])) + e_k = - eis_H,k(q) - SUM a_1(E)^(k-2)/#AutE
    # eis_H,k(q) = - SUM a_1(E)^(k-2)/#AutE - Tr(F_q | S([H,k+2]))
    # eis_H,k(q) = - S_Gamma_Y - Tr_Tq
    # this is the expected eis_H_k, we compute it now to be able to find a formula to compute this, and therefore compute Tr

    # S_Gamma_Y = 501120

    eis_H_k = - S_Gamma_Y - Tr_Tq

    phi_N = euler_phi(N)
    sgn = 1 if q % N == 1 else (-1 if q % N == N - 1 else 0)
    sgn_k = sgn**(k+2) #avoid k == 0 sgn**k = 1 issue

    # TODO: make it work for N < 5

    if k == 0 and N == 1:
        cusp_gamma1_1 = -q
    else:
        cusp_gamma1_1 = cusp_term(p, n, N, k+2)
    # else:
    #    cusp_gamma1_1 = -(q if N == 1 else q-1)
    # if N == 2:
    # cusp_gamma1_1 = (phi_N // 2) * (1 + sgn_k)
    # if N == 4:
    #    cusp_gamma1_1 = (3 / 2) * (1 + sgn_k)

    return HeckeResult(
        p=p,
        n=n,
        q_equiv_N=-1 if q % N == N - 1 and N != 2 else int(q % N),
        N_equiv_p=-1 if N % p == p - 1 and p != 2 else int(N % p),
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
        N_eisen=N_eis,
        N_gauss=N_gauss,
        j0_SS=j0_SS,
        j1728_SS=j1728_SS,
        DIM=NUM_TRACES
    )


def print_result(result: HeckeResult):
    if result:
        header = f" {'p^n':>6} {'n':>6} {'q≡N':>6} {'NC':>6} {'NSS':>6} {'j0_SS':>10} {'j1728_SS':>10} {'NEis':>6} {'NGauss':>6} {'q':>10} {'Computed':>20} {'Tr(T_q)':>20} {'eis':>10} {'cusp_est':>10} {'DIM':>6} {'Error':>30}"
        print(header)
        print("-" * len(header))
        for r in results:
            adjusted_val_1 = r.eis - r.error_1
            clr = Colors.BOLD
            if adjusted_val_1 == 0:
                clr = Colors.GREEN
            # else:
            #    clr = Colors.WARNING
            if r.N_equiv_p == 0:
                clr = Colors.FAIL
            print(
                f" {clr}{r.p:>6} {r.n:>6} {r.q_equiv_N:>6} {int(r.NC):>6} {float(r.NSS):>6.1f} {r.j0_SS:>10} {r.j1728_SS:>10} {r.N_eisen:>6} {r.N_gauss:>6} {(r.p**r.n):>10} {r.computed_sum:>20} {r.Tr:>20} {r.eis:>10} {r.error_1:>10} {r.DIM:>6} {fmt_factored(int(adjusted_val_1)):>30} {Colors.ENDC}"
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
    
    if args.compare and args.N > 1:
        M = ModularForms(Gamma1(args.N), args.k + 2)
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
        print(f"Dimension of M_{args.k+2}(Gamma1({args.N})): {dim_M} (Eisenstein: {dim_E}, Cusp: {dim_S}, Old: {dim_S_old}, New: {dim_S_new}) sk={dim_sk(args.k+2)}\n")

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
