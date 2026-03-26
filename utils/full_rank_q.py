"""
Print q = p^n (p prime, n >= 1) where some prime ell < 100 divides q - 1,
AND the Hasse bound interval [q+1-HB, q+1+HB] contains a multiple of ell
(i.e., there exist actual curves over F_q with full ell-torsion).
"""
from sympy import primerange
from math import isqrt, log
import argparse

from common import Colors

parser = argparse.ArgumentParser(description="Find q=p^n with full ell-torsion potential.")
parser.add_argument("--p-max",   type=int,   default=10**6,  help="Primes p up to this (default: 10^6)")
parser.add_argument("--n-max",   type=int,   default=5,      help="Extension degrees n to consider (default: 5)")
parser.add_argument("--ell-max", type=int,   default=101,    help="Torsion primes ell < this (default: 101)")
parser.add_argument("--ell-min", type=int,   default=90,     help="Only consider ell > this for interesting cases (default: 90)")
parser.add_argument("--q-max",   type=int,   default=10**16, help="Skip q above this (default: 10^16)")
args = parser.parse_args()

P_MAX   = args.p_max
N_MAX   = args.n_max
ELL_MAX = args.ell_max
ELL_MIN = args.ell_min
Q_MAX   = args.q_max

primes   = [p for p in primerange(5, P_MAX)]  # exclude char 2 and 3
ell_list = list(primerange(2, ELL_MAX))

def hasse_valid(ell: int, q: int) -> bool:
    """True if there exists an *ordinary* t in the Hasse interval with ell^2 | q+1-t.
    (Rank-2 ell-torsion requires ell^2 | #E = q+1-t and disc = t^2-4q != 0.)"""
    HB = isqrt(4 * q)
    ell2 = ell * ell
    # t ≡ q+1 (mod ell^2)  →  i = (q+1-t)/ell^2 is an integer
    i_min = (q + 1 - HB + ell2 - 1) // ell2  # ceil
    i_max = (q + 1 + HB) // ell2              # floor
    if i_min > i_max:
        return False
    for i in range(i_min, i_max + 1):
        t = q + 1 - i * ell2
        if t * t - 4 * q != 0:   # exclude disc=0 (degenerate/boundary case)
            return True
    return False

print(f"{'q':>12}  {'p':>4}  {'n':>2}  ell dividing q-1 (Hasse-valid)")
print("-" * 65)
results = []
for p in primes:
    for n in range(1, N_MAX + 1):
        q = p ** n
        if q > Q_MAX:
            break
        #print(f"Checking F_{q}, p={p}, n={n}... ", end="")
        dividing = [ell for ell in ell_list if ell > ELL_MIN and (q - 1) % ell == 0 and hasse_valid(ell, q)]
        if dividing:
            results.append((q, p, n, dividing))

results.sort(key=lambda r: max(r[3]), reverse=True)


def complexity_note(q: int, ell: int) -> str:
    """Return side-by-side asymptotic notes for two rank-test costs."""
    l2 = ell * ell
    logq = max(log(q), 1.0)
    logl = max(log(ell), 1.0)
    c1 = l2 + ell * logq
    c2 = l2 * logl * logq
    ratio = c2 / c1 if c1 else float("inf")
    ratio_text = f"×{ratio:.2f}"
    if ratio < 1:
        ratio_text = f"{Colors.FAIL}{ratio_text}{Colors.ENDC}"
    elif ratio > 100:
        ratio_text = f"{Colors.GREEN}{ratio_text}{Colors.ENDC}"
    return (
        f"ℓ={ell}: {c1:.1f} vs {c2:.1f} ({ratio_text})"
    )


for q, p, n, dividing in results:
    notes = "; ".join(complexity_note(q, ell) for ell in dividing)
    print(f"{q:>12}  {p:>4}  {n:>2}  {dividing}  |  {notes}")
