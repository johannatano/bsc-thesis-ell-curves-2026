# from sage.all import *
from sage.all import *
import random
from utils.common import Colors
import argparse
from typing import Optional, List, Dict, Tuple, Set, Any
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
# from sage.all import CuspForms, Gamma1

from nt.common import legendre

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


class QF:
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


def __D0(D: int) -> int:
    fund_dscr = -1  # will accumulate the odd-exponent primes
    fac_dscr = list(factorint(dscr).items())
    print(f"factorization of D={D} is {fac_dscr}")
    for (prime, exp) in fac_dscr:
        if exp % 2 == 0:
            cnd *= prime ** (exp // 2)
        else:
            cnd *= prime ** ((exp - 1) // 2)
            fund_dscr *= prime

    # Adjust so that fund_dscr is a true fundamental discriminant (≡ 0 or 1 mod 4).
    # If fund_dscr ≡ 2 or 3 (mod 4), we need to absorb an extra factor of 4.
    if fund_dscr % 4 == 2 or fund_dscr % 4 == 3:
        fund_dscr *= 4
        cnd //= 2  # compensate: (2*cnd')^2 * fund_
    return fund_dscr

def run(p, n, N, k):
    q = p**n

    SS_poly = supersingular_j_polynomial(p)
    SS_poly_Fq = SS_poly.change_ring(GF(q))

    ss_js = []
    for r in SS_poly_Fq.roots(multiplicities=False):
        ss_js.append(r)

    print(f"\n{Colors.FAIL}FOUND supersingular curve with js={ss_js}{Colors.ENDC}")

    num_generic_quat = len(ss_js)
    num_quat_j0 = Fraction(1 - legendre(-3, p), 2)
    num_quat_j1728 = Fraction(1 - legendre(-4, p), 2)
    
    print(f"legendre num_j0_quat={num_quat_j0}, num_j1728_quat={num_quat_j1728}")

    num_j0_quat = 1 if p % 3 == 2 else 0
    num_j1728_quat = 1 if p % 4 == 3 else 0

    num_quat = Fraction(1, 12)*(p + 6 - 4 * legendre(-3 , p) - 3 * legendre(-4 , p))

    w_num_generic = Fraction(num_quat - num_j0_quat - num_j1728_quat, 2)
    w_j0 = Fraction(num_j0_quat, 6)
    w_j1728 = Fraction(num_j1728_quat, 4)

    w_total = w_num_generic + w_j0 + w_j1728

    simplified = Fraction(p-1,24)

    print(
        f"{Colors.BRIGHT_GREEN}|| NUM QUAT = {num_quat},  num_generic_ss={num_generic_quat}, num_j0_quat={num_j0_quat}, num_j1728_quat={num_j1728_quat}{Colors.ENDC} \n"
    )

    print(
        f"{Colors.BRIGHT_CYAN}|| NUM QUAT W = {w_total},  w_num_generic={w_num_generic}, w_j0={w_j0}, w_j1728={w_j1728}, simplified={simplified}{Colors.ENDC} \n"
    )

    j0 = EllipticCurve(GF(q), j=0)
    j1728 = EllipticCurve(GF(q), j=1728)

    '''print(
        f"\n j0={j0.abelian_group().invariants()}"
    )
    print(
        f"\n j1728={j1728.abelian_group().invariants()}"
    )

    print(
        f"\n iso ={j0.is_isomorphic(j1728)}, aut_size j0={len(j0.automorphisms())}, aut_size j1728={len(j1728.automorphisms())}"
    )'''

    is_j0_SS = j0.is_supersingular()
    is_j1728_SS = j1728.is_supersingular()

    '''print(
        f"{Colors.HEADER}|| gcd(4, p-1)={gcd(4, q-1)}, p equiv 3 mod 4={p % 4 == 3}, j1728_SS={is_j1728_SS} ||  gcd(6, p-1)={gcd(6, q-1)}, p equiv 2 mod 3={p % 3 == 2}, j0_SS={is_j0_SS}, || H(-4p)={pari(4*p).qfbclassno()} ie nr of curves with t=0, HKclass={QF.H(-4*p)}{Colors.ENDC}"
    )'''

    for E in j0.twists():
        t = E.trace_of_frobenius()
        D_pi = t**2 - 4*q
        if D_pi == 0:
            D_K = 0
            f_pi = 0
        else:
            K = QuadraticField(D_pi)
            D_K = K.discriminant()
            f_pi = math.isqrt(abs(D_pi // D_K)) if D_K != 0 else 0

        clr = Colors.GREEN if E.is_isomorphic(j0) else Colors.BOLD
        if D_K == 0:
            clr = Colors.WARNING
        if t == 0:
            clr = Colors.BLUE
        '''print(
            f"{clr}j0 | t={t}, , D_K={D_K}, f_pi={f_pi}, invariants={E.abelian_group().invariants()}, f_pi={f_pi}, iso to j0={E.is_isomorphic(j0)}, aut_size={len(E.automorphisms())}{Colors.ENDC}"
        )'''
        if f_pi != 0:
            for f in ZZ(f_pi).divisors():
                if f % p == 0:
                    continue
                j_invs = []

                D = D_K * f**2
                H = hilbert_class_polynomial(D)
                H_fq = H.change_ring(GF(q))
                for j, m in H_fq.roots(multiplicities=True):
                    for _ in range(m):
                        j_invs.append(j)
                inert = 2 if kronecker(D, p) == -1 else 1
                print(
                    f"j0 | t={t}, D_K={D_K}, f={f}, inert={inert}, N(t)={1 - kronecker(D, p)}, kronecker(D, p)={kronecker(D, p)}, D={D}"
                )

    if p >= 5:
        for E in j1728.twists():
            t = E.trace_of_frobenius()
            D_pi = t**2 - 4*q
            if D_pi == 0:
                D_K = 0
                f_pi = 0
            else:
                K = QuadraticField(D_pi)
                D_K = K.discriminant()
                f_pi = math.isqrt(abs(D_pi // D_K)) if D_K != 0 else 0

            clr = Colors.GREEN if E.is_isomorphic(j1728) else Colors.BOLD
            if D_K == 0:
                clr = Colors.WARNING

            if t == 0:
                clr = Colors.BLUE
            print(
                f"{clr}j1728 | t={t}, , D_K={D_K}, f_pi={f_pi}, invariants={E.abelian_group().invariants()}, f_pi={f_pi}, iso to j1728={E.is_isomorphic(j1728)}, aut_size={len(E.automorphisms())}{Colors.ENDC}"
            )
            if f_pi != 0:
                for f in ZZ(f_pi).divisors():
                    if f % p == 0:
                        continue
                    j_invs = []

                    D = D_K * f**2
                    H = hilbert_class_polynomial(D)
                    H_fq = H.change_ring(GF(q))
                    for j, m in H_fq.roots(multiplicities=True):
                        for _ in range(m):
                            j_invs.append(j)
                    inert = 2 if kronecker(D, p) == -1 else 1
                    print(
                        f"j1728 | t={t}, D_K={D_K}, f={f}, inert={inert}, N(t)={1 - kronecker(D, p)}, kronecker(D, p)={kronecker(D, p)}, D={D}"
                    )


if __name__ == "__main__":
    args = parse_args()
    run(args.p, args.n, args.N, args.k)
