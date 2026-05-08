# Translation of reference_maple.py from Maple to Python (for readability only — not runnable as-is).
#
# Maple setup notes:
#   - `read "classnmb"` loads a precomputed table LstClnmb of class numbers of fundamental discriminants.
#   - `with(numtheory)` gives: ifactors (integer factorization), divisors, legendre (Legendre symbol).
#   - `with(SF)` loads the Symmetric Functions package (used only in charpol_Sk).
#   - ifactors(n)[2] returns a list of [prime, exponent] pairs (1-indexed in Maple, 0-indexed here).
#   - nops(L) = len(L)
#   - foldl(`+`, 0, seq(...)) = sum(...)
#   - foldl(`*`, 1, seq(...)) = math.prod(...)
#   - modp(a, m) = a % m
#   - a &^ b in Maple = pow(a, b, m) (modular exponentiation, context-dependent)
#   - convert(n, base, b) = digits of n in base b, least significant first

from fractions import Fraction
from math import floor, sqrt, comb
from sympy import factorint, divisors, legendre_symbol, isprime

# LstClnmb is a preloaded dict: LstClnmb[D] = class number h(D) for fundamental discriminant D.
# Loaded from the external "classnmb" file — must be provided separately.
LstClnmb: dict  # e.g. {-3: Fraction(1,1), -4: Fraction(1,1), -7: Fraction(1,1), ...}

# A1q[q] stores the result of count_A1q(p, r) as a dict {a: weighted_count}.
A1q: dict = {}


# ---------------------------------------------------------------------------
# HKclass(dscr)
# Computes the Hurwitz-Kronecker class number H(dscr).
# Uses the convention that weights each class by 1/|Aut|, so H(-3)=1/3, H(-4)=1/2.
# Requires dscr < 0 and dscr ≡ 0 or 1 (mod 4).
# ---------------------------------------------------------------------------
def HKclass(dscr) -> Fraction:
    # Factor |dscr| into a list of (prime, exponent) pairs
    fac_dscr = list(factorint(dscr).items())  # [(p1, e1), (p2, e2), ...]

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
        clnmb = Fraction(LstClnmb[fund_dscr] * d)
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
        return

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
                p + 6 - 4 * legendre_symbol(-3, p) - 3 * legendre_symbol(-4, p), 24
            )
            # Corrections for the j=0 (CM by Z[ω], extra 6 auts) and j=1728 (CM by Z[i], extra 4 auts) curves
            if p % 3 != 1:  # p is inert or ramified at 3 → j=0 curve appears here
                res += Fraction(-1, 2) + Fraction(1, 6)
            if p % 4 != 1:  # p is inert or ramified at 2 → j=1728 curve appears here
                res += Fraction(-1, 2) + Fraction(1, 4)

        elif a**2 == q and r % 2 == 0:
            # Case 8: a=sqrt(q), r even — involves only -3 Legendre symbol
            res = Fraction(1 - legendre_symbol(-3, p), 6)

        elif a == 0 and r % 2 == 0:
            # Case 9: a=0, r even — involves only -4 Legendre symbol
            res = Fraction(1 - legendre_symbol(-4, p), 4)

        else:
            res = Fraction(0)  # empty isogeny class

        Res[a] = res
        Res[-a] = (
            res  # trace a and -a give isomorphic (dual) isogeny classes with same count
        )

        # Count a=0 once, all other ±a pairs count twice
        total += res if a == 0 else 2 * res

    if total != q:
        print("Mistake! Total does not equal q.")

    A1q[q] = Res


# ---------------------------------------------------------------------------
# count_A1q_all(strt, fin)
# Runs count_A1q for every prime power q with strt <= q <= fin.
# ---------------------------------------------------------------------------
def count_A1q_all(strt: int, fin: int):
    for n in range(strt, fin + 1):
        fac = factorint(n)
        if len(fac) == 1:  # n is a prime power
            [(p, r)] = fac.items()
            print(f"q = {p**r}")
            count_A1q(p, r)


# ---------------------------------------------------------------------------
# count_A1q_Q(Q)
# Runs count_A1q for every prime power in the list Q.
# ---------------------------------------------------------------------------
def count_A1q_Q(Q: list):
    for n in Q:
        fac = factorint(n)
        if len(fac) == 1:
            [(p, r)] = fac.items()
            print(f"q = {p**r}")
            A1q[n] = count_A1q(p, r)


# ---------------------------------------------------------------------------
# sk(k)
# Dimension of the space S_k of cuspidal elliptic modular forms of weight k for SL_2(Z).
# Formula from Riemann-Roch / Gauss-Bonnet on the modular curve.
# ---------------------------------------------------------------------------
def sk(k: int) -> int:
    if k < 0 or k % 2 == 1:
        return 0
    elif k == 2:
        return -1  # S_2 has dimension 0, but the formula gives -1 (boundary case)
    elif k % 12 == 2:
        return floor(k / 12) - 1
    else:
        return floor(k / 12)


# ---------------------------------------------------------------------------
# hk_poly(a, k, q) — helper
# Evaluates the "Hecke polynomial" (characteristic polynomial of Frobenius π on an
# elliptic curve with trace a over F_q, at the power k) at e1=a.
#
#   hk_poly(a, k, q) = sum_{i=0}^{k//2} (-1)^i * q^i * C(k-i, i) * a^{k-2i}
#
# This is the degree-k Newton polynomial in α,β where α+β=a, αβ=q.
# It equals α^k + β^k where α,β are the Frobenius eigenvalues.
# (Same recurrence as models.py hk(k) method.)
# ---------------------------------------------------------------------------
def hk_poly(a: int, k: int, q: int) -> int:
    return sum(
        (-1) ** i * q**i * comb(k - i, i) * a ** (k - 2 * i) for i in range(k // 2 + 1)
    )


# ---------------------------------------------------------------------------
# count_Sk(p, r, kk)
# Computes Tr(Frob_{p^r} | S_{kk}) using the Eichler-Selberg trace formula.
# Same 9-case structure as count_A1q, but each isogeny class contributes
# weighted by hk_poly(a, kk-2, q) instead of just 1.
#
# The hk_poly factor is the trace of Frobenius on the kk-weight motive of an
# elliptic curve with trace a — i.e., it sums α^{kk-2} + β^{kk-2} over the
# two eigenvalues.
# ---------------------------------------------------------------------------
def count_Sk(p: int, r: int, kk: int) -> Fraction:
    q = p**r
    k = kk - 2  # weight shift: S_{kk} corresponds to degree-k Hecke polynomial

    if q >= 10**7 / 4:
        print("q too large")
        return Fraction(0)

    if kk % 2 == 1:
        return Fraction(0)  # odd weight: trace is always 0 by symmetry

    Res = Fraction(0)
    amax = floor(2 * sqrt(q))

    for a in range(0, amax + 1):
        w = hk_poly(a, k, q)  # weight factor for this isogeny class

        if a % p != 0:
            # Case 1: factor of 2 because we sum over a and -a together
            res = HKclass(a**2 - 4 * q) * 2 * w

        elif a == 0 and r % 2 == 1:
            # Case 2: no factor of 2 since a=0 only appears once
            res = HKclass(-4 * p) * w

        elif a**2 == 2 * q and p == 2 and r % 2 == 1:
            res = Fraction(1, 4) * 2 * w  # case 3

        elif a**2 == 3 * q and p == 3 and r % 2 == 1:
            res = Fraction(1, 6) * 2 * w  # case 4

        elif a**2 == 4 * q and r % 2 == 0 and p == 2:
            res = Fraction(1, 24) * 2 * w  # case 5

        elif a**2 == 4 * q and r % 2 == 0 and p == 3:
            res = Fraction(1, 12) * 2 * w  # case 6

        elif a**2 == 4 * q and r % 2 == 0:
            # Case 7: same correction formula as count_A1q case 7
            weight = Fraction(
                p + 6 - 4 * legendre_symbol(-3, p) - 3 * legendre_symbol(-4, p), 24
            )
            if p % 3 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 6)
            if p % 4 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 4)
            res = weight * 2 * w

        elif a**2 == q and r % 2 == 0:
            res = Fraction(1 - legendre_symbol(-3, p), 6) * 2 * w  # case 8

        elif a == 0 and r % 2 == 0:
            res = Fraction(1 - legendre_symbol(-4, p), 4) * w  # case 9

        else:
            res = Fraction(0)

        Res += res

    # The -1 at the end is the contribution of the identity class (the trivial representation)
    return -(Res + 1)


# ---------------------------------------------------------------------------
# count_Sk_all(strt, fin, k)
# Returns list of Tr(Frob_{p^r} | S_k) for all prime powers q=p^r in [strt, fin].
# ---------------------------------------------------------------------------
def count_Sk_all(strt: int, fin: int, k: int) -> list:
    traces = []
    prime_powers = []
    for n in range(strt, fin + 1):
        fac = factorint(n)
        if len(fac) == 1:
            [(p, r)] = fac.items()
            prime_powers.append(p**r)
            traces.append(count_Sk(p, r, k))
    print("Q", prime_powers)
    return traces


# ---------------------------------------------------------------------------
# charpol_Sk(pp, k)
# Computes the characteristic polynomial of Frob_p acting on S_k.
#
# Strategy (uses the SF symmetric functions package in Maple):
#   dim = sk(k)                          — dimension of S_k
#   qq  = pp^{k-1}                       — functional equation scaling factor
#   L   = [Tr(Frob_pp^i | S_k) for i=1..dim]   — power sums p_1,...,p_dim
#
#   Convert power sums → elementary symmetric polynomials e_1,...,e_dim
#   via Newton's identities:
#       i * e_i = sum_{j=1}^{i} (-1)^{j-1} e_{i-j} p_j
#
#   The characteristic polynomial is then (using the Weil palindrome):
#       f(x) = x^{2d} + sum_{i=1}^{d}   (-1)^i e_i  x^{2d-i}
#                     + sum_{i=1}^{d-1} (-1)^i e_i qq^{d-i}  x^i
#                     + qq^d
#
# In Maple: `top(cat(e,j))` retrieves the j-th elementary symmetric polynomial
#           after substituting the power sums via the SF package.
# ---------------------------------------------------------------------------
def charpol_Sk(pp: int, k: int):
    dim = sk(k)
    qq = pp ** (k - 1)

    # Power sums: p_i = Tr(Frob_{pp}^i | S_k) for i = 1 ... dim
    power_sums = [count_Sk(pp, i, k) for i in range(1, dim + 1)]

    # Newton's identities: compute elementary symmetric polynomials e_1,...,e_dim
    # from power sums p_1,...,p_dim.
    # e[0] = 1 by convention; e[i] corresponds to e_i.
    e = [Fraction(0)] * (dim + 1)
    e[0] = Fraction(1)
    for i in range(1, dim + 1):
        e[i] = Fraction(1, i) * sum(
            (-1) ** (j - 1) * e[i - j] * power_sums[j - 1] for j in range(1, i + 1)
        )

    # Build the palindromic characteristic polynomial (as a list of coefficients):
    #   f(x) = x^{2d} + sum_{i=1}^{d} (-1)^i e_i x^{2d-i}
    #                 + sum_{i=1}^{d-1} (-1)^i e_i qq^{d-i} x^i
    #                 + qq^d
    # Coefficients indexed by degree (degree 0 = constant term):
    coeffs = [Fraction(0)] * (2 * dim + 1)
    coeffs[2 * dim] = Fraction(1)  # x^{2d}
    coeffs[0] = Fraction(qq**dim)  # constant term qq^d
    for i in range(1, dim + 1):
        coeffs[2 * dim - i] += (-1) ** i * e[i]
    for i in range(1, dim):
        coeffs[i] += (-1) ** i * e[i] * qq ** (dim - i)

    return coeffs  # coeffs[j] is the coefficient of x^j


# ---------------------------------------------------------------------------
# count_ak(p, r, kk)
# Variant of count_Sk — computes a simpler trace using a^kk as the weight
# instead of the full Hecke polynomial hk_poly(a, kk-2, q).
# (Used for computing a_p(f)^k directly rather than the Frobenius trace on the
#  symmetric power representation.)
# ---------------------------------------------------------------------------
def count_ak(p: int, r: int, kk: int) -> Fraction:
    q = p**r

    if q >= 10**7 / 4:
        print("q too large")
        return Fraction(0)

    if kk % 2 == 1:
        return Fraction(0)

    Res = Fraction(0)
    amax = floor(2 * sqrt(q))

    for a in range(0, amax + 1):
        w = a**kk  # simpler weight: just a^kk

        if a % p != 0:
            res = HKclass(a**2 - 4 * q) * 2 * w  # case 1

        elif a == 0 and r % 2 == 1:
            res = HKclass(-4 * p) * w  # case 2 (w=0 since a=0, so res=0 always here)

        elif a**2 == 2 * q and p == 2 and r % 2 == 1:
            res = Fraction(1, 4) * 2 * w  # case 3

        elif a**2 == 3 * q and p == 3 and r % 2 == 1:
            res = Fraction(1, 6) * 2 * w  # case 4

        elif a**2 == 4 * q and r % 2 == 0 and p == 2:
            res = Fraction(1, 24) * 2 * w  # case 5

        elif a**2 == 4 * q and r % 2 == 0 and p == 3:
            res = Fraction(1, 12) * 2 * w  # case 6

        elif a**2 == 4 * q and r % 2 == 0:
            weight = Fraction(
                p + 6 - 4 * legendre_symbol(-3, p) - 3 * legendre_symbol(-4, p), 24
            )
            if p % 3 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 6)
            if p % 4 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 4)
            res = weight * 2 * w  # case 7

        elif a**2 == q and r % 2 == 0:
            res = Fraction(1 - legendre_symbol(-3, p), 6) * 2 * w  # case 8

        elif a == 0 and r % 2 == 0:
            res = Fraction(1 - legendre_symbol(-4, p), 4) * w  # case 9

        else:
            res = Fraction(0)

        Res += res

    return Res


# ---------------------------------------------------------------------------
# bincoeff_lucas(m, n, ell)
# Computes C(m, n) mod ell via Lucas's theorem (valid when ell is prime).
# Lucas: C(m,n) ≡ prod_i C(m_i, n_i) (mod ell)
# where m = sum m_i * ell^i and n = sum n_i * ell^i are base-ell expansions.
# ---------------------------------------------------------------------------
def bincoeff_lucas(m: int, n: int, ell: int) -> int:
    if n > m:
        return 0

    def to_base(x: int, b: int) -> list:
        # Returns digits of x in base b, least significant first
        digits = []
        while x > 0:
            digits.append(x % b)
            x //= b
        return digits

    M = to_base(m, ell)
    N = to_base(n, ell)
    N += [0] * (len(M) - len(N))  # pad with leading zeros to match length

    # Multiply the digit-wise binomial coefficients mod ell
    res = 1
    for mi, ni in zip(M, N):
        res = (res * comb(mi, ni)) % ell
    return res


# ---------------------------------------------------------------------------
# count_Sk_mod(p, r, kk, m)
# Computes Tr(Frob_{p^r} | S_{kk}) mod m.
# Same as count_Sk but works mod m throughout to avoid large integers.
# If m is prime, uses Lucas's theorem for binomial coefficients in the Hecke polynomial.
# ---------------------------------------------------------------------------
def count_Sk_mod(p: int, r: int, kk: int, m: int) -> int:
    q = p**r
    k = kk - 2

    # Precompute the coefficients of the Hecke polynomial mod m:
    #   Pol[i] = (-1)^i * q^i * C(k-i, i)  mod m
    if isprime(m):
        # Lucas's theorem for C(k-i, i) mod prime m
        Pol = [
            pow(-1, i) * pow(q, i, m) * bincoeff_lucas(k - i, i, m) % m
            for i in range(k // 2 + 1)
        ]
    else:
        Pol = [
            pow(-1, i) * pow(q, i, m) * comb(k - i, i) % m for i in range(k // 2 + 1)
        ]

    if q >= 10**7 / 4:
        print("q too large")
        return 0

    if kk % 2 == 1:
        return 0

    Res = 0
    amax = floor(2 * sqrt(q))

    for a in range(0, amax + 1):
        # Evaluate Hecke polynomial at a, mod m
        # TrE = sum_{i=0}^{k//2} Pol[i] * a^{k-2i} mod m
        TrE = sum(Pol[i] * pow(a, k - 2 * i, m) for i in range(k // 2 + 1)) % m

        if a % p != 0:
            res = int(HKclass(a**2 - 4 * q) * 2 * TrE) % m  # case 1

        elif a == 0 and r % 2 == 1:
            res = int(HKclass(-4 * p) * TrE) % m  # case 2

        elif a**2 == 2 * q and p == 2 and r % 2 == 1:
            res = int(Fraction(1, 2) * TrE) % m  # case 3: (1/4)*2 = 1/2

        elif a**2 == 3 * q and p == 3 and r % 2 == 1:
            res = int(Fraction(1, 3) * TrE) % m  # case 4: (1/6)*2 = 1/3

        elif a**2 == 4 * q and r % 2 == 0 and p == 2:
            res = int(Fraction(1, 12) * TrE) % m  # case 5: (1/24)*2 = 1/12

        elif a**2 == 4 * q and r % 2 == 0 and p == 3:
            res = int(Fraction(1, 6) * TrE) % m  # case 6: (1/12)*2 = 1/6

        elif a**2 == 4 * q and r % 2 == 0:
            weight = Fraction(
                p + 6 - 4 * legendre_symbol(-3, p) - 3 * legendre_symbol(-4, p), 24
            )
            if p % 3 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 6)
            if p % 4 != 1:
                weight += Fraction(-1, 2) + Fraction(1, 4)
            res = int(weight * 2 * TrE) % m  # case 7

        elif a**2 == q and r % 2 == 0:
            res = int(Fraction(1 - legendre_symbol(-3, p), 6) * 2 * TrE) % m  # case 8

        elif a == 0 and r % 2 == 0:
            res = int(Fraction(1 - legendre_symbol(-4, p), 4) * TrE) % m  # case 9

        else:
            res = 0

        Res = (Res + res) % m

    return -(Res + 1) % m
