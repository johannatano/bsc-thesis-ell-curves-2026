"""Enumeration and classification of elliptic curves over finite fields.

This module implements three complementary workflows:
- direct enumeration of j-invariants and their twist families,
- CM/Hilbert-class-polynomial enumeration of curves by endomorphism order,
- class-number-based counting when only aggregate arithmetic data is needed.

It organizes the resulting curves into trace-indexed isogeny classes and
supports downstream computations such as torsion data and isogeny volcanoes.
"""

from sage.all import *
from math import gcd
import math
import hashlib
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Set, Any

from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from sympy import primerange

from utils.common import Logger, Colors
from lib.ell_torsion_subgroup import *
from lib.curves import *
from lib.nr_fields import *
from lib.nr_fields_classifier import NumberFieldCatalogue

# from utils.mod_poly import _classical_modular_polynomial

highest_ell = 0

# =============================================================================
# Curve enumeration and classification
# =============================================================================

'''
def num_order_N(N, n1, n2):
    result = 1
    for l, a in factor(ZZ(N)):
        e1 = ZZ(n1).valuation(l)
        e2 = ZZ(n2).valuation(l)
        result *= l ** (e1 + e2) - l ** (min(a - 1, e1) + min(a - 1, e2))
    return result'''

'''
def num_order_N(N, n1, n2):
    if n2 % N != 0:
        return 0  # max order is n2; if N doesn't divide it, nothing
    # else: count elements of order exactly N in Z/n1 x Z/n2
    return phi(N) * gcd(N, n1)'''


def phi(n):
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

from math import gcd


def num_order_N(N, n1, n2):
    if n2 % N != 0:
        return 0
    # distinct primes of N
    primes = []
    n = N
    p = 2
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes.append(n)

    # inclusion-exclusion
    total = 0
    r = len(primes)
    for mask in range(1 << r):
        d = N
        bits = 0
        for i in range(r):
            if mask & (1 << i):
                d //= primes[i]
                bits += 1
        sign = -1 if bits % 2 else 1
        total += sign * gcd(d, n1) * gcd(d, n2)
    return total


'''
def l_sylow_structure(l, f_pi, f, N_pts):
    """Returns (s1, s2) with s1 <= s2, s1 + s2 = v_l(#E).
    Full l-Sylow of E(Fq), no capping at N."""

    s1 = 0
    vl_pi = ZZ(f_pi).valuation(l)
    vl_pts = ZZ(N_pts).valuation(l)

    for k in range(1, vl_pts // 2 + 1):
        if ZZ(f).valuation(l) + k <= vl_pi:  # E[l^k] subset E(Fq)
            s1 = k
        else:
            break
    s2 = vl_pts - s1
    return s1, s2'''


def l_sylow_structure(l, a, f_pi, f, N_pts, q):
    """Returns (s1, s2) capped at a."""
    vl_pts = ZZ(N_pts).valuation(l)
    vl_pi = ZZ(f_pi).valuation(l)
    vl_f = ZZ(f).valuation(l)
    # Hard upper bound on s1
    max_s1 = min(a, vl_pts // 2, ZZ(q - 1).valuation(l))
    # Embedding height
    h = vl_pi - vl_f
    s1 = min(max_s1, h)  # embedding gives at most h, capped at max_s1
    s2 = min(a, vl_pts - s1)
    return s1, s2


def invariants(N, f_pi, f, N_pts, q):
    """Returns (n1, n2) with n1 | n2, n1*n2 = #E[N](Fq)
    such that E[N](Fq) cap E(Fq) = Z/n1 x Z/n2."""
    n1, n2 = 1, 1
    for l in ZZ(N).prime_factors():
        a = ZZ(N).valuation(l)
        s1, s2 = l_sylow_structure(l, a, f_pi, f, N_pts, q)
        n1 *= l**s1
        n2 *= l**s2
    return n1, n2


def num_P(level, f_pi, f, N_pts, q):
    if level == 1:
        return 1
    result = 1
    for l, a in factor(ZZ(level)):
        
        h = max(0, ZZ(f_pi).valuation(l) - ZZ(f).valuation(l)) if ZZ(f_pi)*ZZ(f) > 0 else None
        v_q1 = ZZ(q - 1).valuation(l)
        v_N = ZZ(N_pts).valuation(l)
        
        e1 = min(h, v_q1, v_N // 2) if h is not None else min(v_q1, v_N // 2)
        e2 = v_N - e1
        s1 = min(a, e1)
        s2 = min(a, e2)
        # exact-order-l^a count in Z/l^s1 x Z/l^s2
        result *= l ** (s1 + s2) - l ** (min(a - 1, s1) + min(a - 1, s2))
    return result


def num_P_SS(level, N_pts, q):
    if level == 1:
        return 1
    result = 1
    for l, a in factor(ZZ(level)):
        v_q1 = ZZ(q - 1).valuation(l)
        v_N = ZZ(N_pts).valuation(l)
        e1 = min(v_q1, v_N // 2)
        e2 = v_N - e1
        s1 = min(a, e1)
        s2 = min(a, e2)
        # exact-order-l^a count in Z/l^s1 x Z/l^s2
        result *= l ** (s1 + s2) - l ** (min(a - 1, s1) + min(a - 1, s2))
    return result


class EllFiniteFieldCatalogue:
    """In-memory catalogue of curves over a fixed finite field F_{p^n}.

    This class acts as the bridge between concrete curve objects and the more
    arithmetic, trace-indexed number-field catalogue.
    """
    def __init__(self, Fq: FqData, NF: Optional['NumberFieldCatalogue'] = None) -> None:
        self.p: int = Fq.p
        self.n: int = Fq.n
        self.q: int = Fq.q
        self.field: FqData = Fq
        self.size: int = 0
        self.NFC: 'NumberFieldCatalogue' = NF if NF is not None else NumberFieldCatalogue(self.p)
    
    def get_isogeny_class(self, t: int, auto_create: bool = True):
        ell_t = self.NFC.get_isogeny_class(t, n=self.n)
        if ell_t is None and auto_create:
            #print(f"{Colors.HEADER}Creating new isogeny class for trace t={t} at extension degree n={self.n}{Colors.ENDC}")
            ell_t = self.NFC.create_isogeny_class(t, n=self.n)
        return ell_t
        
    def add(self, curve: Curve) -> None:
        """Insert a curve into the isogeny class determined by its trace."""
        t = curve.t
        ell_t = self.get_isogeny_class(t)
        ell_t.add_curve(curve)
        self.size += 1
        
    def isogeny_classes(self) -> List:
        return self.NFC.get_isogeny_classes_by_n(self.n)

    def toJSON(self) -> Dict[str, List]:
        """Serialize the catalogue grouped by imaginary quadratic discriminant."""
        return {
            "number_fields": [{
                "D": int(nf_info.discriminant),
                "tree": [tree.toJSON(include_curves=True) for tree in nf_info.tree]
            } for dk, nf_info in self.NFC.data.items() ]
        }


class CurvesClassifier_Fq:
    """Main classifier for enumerating and organizing elliptic curves over F_q.
    
    Supports two enumeration methods:
    1. HCP/CM method: builds lightweight curve records from class polynomials.
    2. Direct method: enumerates j-invariants and constructs concrete twists.

    The classifier is the main entry point used by scripts that want either a
    full geometric catalogue or only the arithmetic data needed downstream.
    """
    def __init__(self, p: int, n: int, NF: Optional['NumberFieldCatalogue'] = None) -> None:
        self.field: FqData = FqData(p, n)
        self.catalogue: EllFiniteFieldCatalogue = EllFiniteFieldCatalogue(self.field, NF=NF)
        self.nth_roots_unity: List[Dict] = []
        size_G = self.field.q - 1
        mu_2 = RootsOfUnity.make_2(self.field.g, size_G)
        mu_4 = RootsOfUnity.make_4(self.field.g, size_G) if gcd(size_G, 4) == 4 else mu_2
        mu_6 = RootsOfUnity.make_6(self.field.g, size_G) if gcd(size_G, 6) == 6 else mu_2

        self.nth_roots_unity = [mu_2, mu_4, mu_6]
        self.HB = math.isqrt(4*self.field.q)

        self.tested_js = set()

    def _get_aut_group_for_j(self, j) -> Dict:
        """Return the relevant roots-of-unity data for the automorphism type of `j`."""
        if j.is_zero() and self.nth_roots_unity[2] is not None:
            return self.nth_roots_unity[2]
        elif (j - 1728).is_zero() and self.nth_roots_unity[1] is not None:
            return self.nth_roots_unity[1]
        else:
            return self.nth_roots_unity[0]

    def enumerate_curves(self, use_HCP: bool = False, use_CN: bool = False, add_SS: bool = True, add_curves:bool = True, special_only:bool = False) -> None:
        """Enumerate curves over the current finite field.

        `use_HCP` switches to the CM/Hilbert-class-polynomial pipeline.
        `use_CN` keeps only class-number counts where possible.
        `add_SS` manually inserts supersingular cases that are not produced by
        the HCP enumeration.
        """
        import time
        from tqdm import tqdm

        # if using CN, we never get j invariants hence no coefficients, this means we have to compute based on t and for t = 0 the uniqueness for twists breaks down, so we need to track which (A,B) pairs have already been claimed by t=0 curves to avoid duplicates
        reset_t0_cache()

        # print(f"{Colors.HEADER}------------------------------Starting curve enumeration for p={self.field.p}, n={self.field.n} F_{self.field.q} with use_HCP={use_HCP}, use_CN={use_CN}, add_SS={add_SS}, add_curves={add_curves}{Colors.ENDC}")

        if use_HCP or use_CN:

            if add_curves:
                nf_list = list(self.catalogue.NFC.data.values())
                _t0 = time.perf_counter()
                for nf in tqdm(nf_list, desc=f"HCP F_{self.field.q}", unit="nf", ncols=80, ascii=True):
                    D_K = nf.discriminant
                    # Orders are grouped by quadratic field and extension degree.
                    orders = nf.getOrders(self.field.n)
                    for order in orders:
                        f = order.conductor
                        # TODO: EDIT BACK IN FALSE TO ALWAYS USE
                        if True:#not use_CN or (D_K in [-3, -4] and int(f) == 1):
                            j_invs = get_j_invariants_from_order(D_K * f**2, f, self.field)
                            for j_inv in j_invs:
                                order.add_j_invariant(j_inv)
                                # The same order may contribute curves to several
                                # trace classes with the same CM field data.
                                for t in order.traces:
                                    self.add_nf_curve(j_inv, t=t, f_E=f)
                        else:
                            # In class-number mode we only track multiplicities.
                            self.catalogue.size += order.class_number*len(order.traces)

            # Supersingular j-invariants are handled separately: the HCP pass is
            # aimed at ordinary CM data and does not cover these cases by itself.
            if add_SS:
                if self.field.p % 3 == 2:
                    print(f"{Colors.FAIL}j-invariant 0 is SS, (6, q-1)={gcd(6, self.field.q-1)}{Colors.ENDC}")
                    self.add_ss_curve_by_j(self.field.F(0))

                if self.field.p % 4 == 3:
                    print(
                        f"{Colors.FAIL}j-invariant 1728 is SS, 1728 mod p = {1728 % self.field.p}, (4, q-1)={gcd(4, self.field.q-1)}{Colors.ENDC}"
                    )
                    self.add_ss_curve_by_j(self.field.F(1728))

                SS_poly = supersingular_j_polynomial(self.field.p)
                SS_poly_Fq = SS_poly.change_ring(self.field.F)

                for r in SS_poly_Fq.roots(multiplicities=False):
                    print(f"{Colors.HEADER}Adding supersingular curve with j={r}{Colors.ENDC}")
                    self.add_ss_curve_by_j(
                        self.field.F(r)
                    )  # we know these are SS, so t=0, and f_E=1 since they have maximal endomorphism ring
            # For each order, recover its j-invariants and attach them to all
            # trace classes compatible with that order.

        elif add_curves:
            _t0 = time.perf_counter()
            precompute_conductor = True

            if special_only:
                self.add_curves_by_j(
                    self.field.F(0), pre_compute_conductor=precompute_conductor
                )
                self.add_curves_by_j(
                    self.field.F(1728), pre_compute_conductor=precompute_conductor
                )
            else:
                n = 0
                for j in tqdm(self.field.F, total=self.field.q, desc=f"F_{self.field.q}", unit="j", ncols=80, ascii=True):
                    self.add_curves_by_j(j, pre_compute_conductor=precompute_conductor)
                    n+=1
                print(f"{Colors.HEADER}Finished enumerating curves over F_{self.field.q}, total size of catalogue: {self.catalogue.size}{Colors.ENDC}")
        global highest_ell
        # Standard count of isomorphism classes over F_q in characteristic > 3.
        NE = 2*(self.field.q -2) + gcd(4, self.field.q-1) + gcd(6, self.field.q-1)

        # print(f"{Colors.HEADER}Finished enumerating curves over F_{self.field.q}, total size of catalogue: {self.catalogue.size}, expected size from formula: {NE}{Colors.ENDC}")
        if self.catalogue.size != NE:
            print(f"{Colors.FAIL}Warning: total number of curves in catalogue ({self.catalogue.size}) does not match expected number from formula ({NE}), there may be duplicates or missing curves{Colors.ENDC}")
        else:
            print(f"{Colors.GREEN}Successfully enumerated curves, total size of catalogue: {self.catalogue.size}{Colors.ENDC}")

        # print(f"Enumeration done in {time.perf_counter() - _t0:.2f}s")

    def check_SS(self, level = 1)->None:
        N = 0
        NQ = 0
        all_js = []
        all_ts = []

        N_pts = 0
        for ell_t in self.catalogue.isogeny_classes():

            if not ell_t.ordinary:
                js = []
                is_Q = False

                t_color = Colors.BOLD
                if ell_t.t == 0:
                    t_color = Colors.GREEN
                elif abs(ell_t.t) == self.HB:
                    t_color = Colors.WARNING
                    is_Q = True
                elif abs(ell_t.t) == math.isqrt(2 * self.field.q) or abs(
                    ell_t.t) == math.isqrt(3 * self.field.q):
                    t_color = Colors.FAIL
                else:
                    t_color = Colors.BLUE

                '''print(
                    f"{t_color}ELL_T t={ell_t.t}, char={self.field.p}, f_pi = {ZZ(ell_t.f_pi).factor() if ell_t.f_pi else None}, D_K = {ZZ(ell_t.D_K).factor() if ell_t.D_K else 0}, O_K = {ell_t.O_K}, N(t)={ell_t.N_t}{Colors.ENDC}"
                )'''

                if level is not None:
                    q_color = Colors.BOLD
                    if self.field.q % level == 1:
                        q_color = Colors.GREEN
                    elif self.field.q % level == level - 1:
                        q_color = Colors.WARNING
                    else:
                        q_color = Colors.BOLD

                    fx_roots = ell_t.fx_pi.roots()
                    r_mod_ell_list = []
                    for root, multiplicity in fx_roots:
                        r_mod_ell = Zmod(level)(root)
                        r_mod_ell_list.append(r_mod_ell)
                    '''print(
                        f"level={level}: fx_roots={ell_t.fx_pi.roots()}, r_mod_ell_list={r_mod_ell_list}, {q_color}q mod level = {self.field.q % level}{Colors.ENDC}, Npts={ell_t.N_pts}, Npts mod level = {ell_t.N_pts % level}"
                    )'''

                # Append colored t value
                all_ts.append(f"{t_color}{ell_t.t}{Colors.ENDC}")

                for f, curves_list in ell_t.curves_by_order.items():
                    for c in curves_list:
                        N += 1
                        if is_Q:
                            NQ += 1
                        js.append(c.j)
                        if c.j not in all_js:
                            all_js.append(c.j)

                        if level is not None:    
                            torsion_subgroup = TorsionSubgroup(c, level)
                            torsion_subgroup.compute_rank(f_pi=-1, use_generators=False)
                            r = torsion_subgroup.rank

                            N_pts += level**r - 1

                            '''print(
                                f"{Colors.FAIL}SS CURVE rank = {r} for curve with j={c.j}, is_j0={(c.j).is_zero()}, is_j1728={(c.j-1728).is_zero()}{Colors.ENDC}"
                            )'''

        if self.field.n % 2 == 0:
            N_exp = num_supersingular_curves_q_square(self.field.p)
            NQ_exp = 2 * quaternion_class_number(self.field.p)
            ts_formatted = "[" + ", ".join(all_ts) + "]"
            if N != N_exp:
                print(
                    f"{Colors.FAIL}Warning: number of supersingular curves found ({N}) does not match expected count from formula ({N_exp}) for q={self.field.q}, js = {all_js}, ts = {ts_formatted}{Colors.ENDC}"
                )
            else:
                print(
                    f"Verified count of supersingular curves: {N} matches expected count from formula for q={self.field.q}, js = {all_js}, ts = {ts_formatted}"
                )

            if NQ != NQ_exp:
                print(
                    f"{Colors.FAIL}Warning: number of supersingular curves with t=±2√q found ({NQ}) does not match expected count from quaternion class number formula ({NQ_exp}) for q={self.field.q}, js = {all_js}, ts = {ts_formatted}{Colors.ENDC}"
                )
            else:
                print(
                    f"{Colors.WARNING}Supersingular curves with t=±2√q found: ({NQ}) matches expected count{Colors.ENDC}"
                )
        else:
            ts_formatted = "[" + ", ".join(all_ts) + "]"
            print(
                f"Found {N} supersingular curves across all traces for q={self.field.q}, js = {all_js}, ts = {ts_formatted}"
            )

        print(
            f"{Colors.GREEN if N_pts > 0 else Colors.BOLD}TOTAL N TORSION POINTS CONTRIBUTED BY SS CURVES: {N_pts}{Colors.ENDC}"
        )

    def add_curves_by_j(self, j, t: Optional[int] = None, f_E: Optional[int] = None, compute_twists: bool = True, pre_compute_conductor=False) -> None:
        """Create the geometric twist family attached to a given j-invariant."""
        # This path materializes full Sage elliptic curves and is therefore the
        # more expensive but more detailed enumeration route.
        aut_grp = self._get_aut_group_for_j(j)
        E = GeometricCurve(self.field, j, aut_grp=aut_grp, t=t, f_E=f_E)

        ell_t = self.catalogue.get_isogeny_class(E.t)
        # All twists in the family share the same conductor, so compute it once.
        if pre_compute_conductor:
            from utils.common import Config
            E.compute_conductor(ell_t.f_pi, use_true_height=Config.use_true_height)
        twists = E.compute_twists() if compute_twists else [E]

        if E.t == 0:
            print(f"{Colors.HEADER}Adding curve with j={j}, t=0, aut_size={E.aut_size}, f_E={E.f_E}, number of twists: {len(twists)}{Colors.ENDC}")

        '''if E.t == 0:
            print(f"{Colors.HEADER}Adding curve with j={j}, t=0, aut_size={E.aut_size}, f_E={E.f_E}, number of twists: {len(twists)}{Colors.ENDC}")'''

        for E_t in twists:
            self.catalogue.add(E_t)
            if E.t == 0:
                print(
                    f"{Colors.HEADER}Added twist with j={E_t.j}, t={E_t.t}, aut_size={E_t.aut_size}, f_E={E_t.f_E}, SS={E_t.is_supersingular}{Colors.ENDC}"
                )

    '''def add_ss_curve_by_j(self, j) -> None:
        """Insert the supersingular twist family attached to `j`."""
        aut_grp = self._get_aut_group_for_j(j)
        E = GeometricCurve(self.field, j, aut_grp=aut_grp, t=None, f_E=1)
        ell_t = self.catalogue.get_isogeny_class(E.t)
        twists = E.compute_twists()
        for E_t in twists:
            self.catalogue.add(E_t)'''

    def add_ss_curve_by_j(self, j) -> None:
        """Insert the supersingular twist family attached to `j`."""
        aut_grp = self._get_aut_group_for_j(j)
        # ell_t = self.catalogue.get_isogeny_class(0)
        E = GeometricCurve(self.field, j, aut_grp=aut_grp, t=None, f_E=None)
        ell_t = self.catalogue.get_isogeny_class(E.t)
        twists = E.compute_twists()
        for E_t in twists:
            self.catalogue.add(E_t)
            print(
                f"{Colors.HEADER}Adding curve with j={j}{Colors.ENDC}"
            )

    def add_nf_curve(self, j, t: int, f_E: int) -> None:
        """Insert a lightweight CM-derived curve record into the catalogue."""
        aut_grp = self._get_aut_group_for_j(j)
        self.catalogue.add(NFCurve(self.field, j, aut_grp=aut_grp, t=t, f_E=f_E))

    def compute_volcano(self, ell:int = -1, edges: bool = False) -> None:
        """Build ℓ-isogeny volcano data for each isogeny class."""
        from tqdm import tqdm
        import time
        primes = [ell] if ell != -1 else list(primerange(2, min(100,max_ell_from_HB(self.field.q)+1)))
        isogeny_classes = self.catalogue.isogeny_classes()
        _t0 = time.perf_counter()
        for ell_t in tqdm(isogeny_classes, desc="computing torsion", unit="ic", ncols=80, ascii=True):
            for ell in primes:
                if ell == self.field.p:
                    continue
                # print(f"{Colors.HEADER}Computing {ell}-isogeny volcano for trace t={ell_t.t}, n={self.field.n}...{Colors.ENDC}")
                ell_t.compute_volcano(ell, edges=edges)

    def count_EP(self, ell, use_CN : bool = False) -> int:
        """Count elliptic points contributing at the given level `ell`."""
        from tqdm import tqdm
        import time
        if ell == self.field.p:
            return 0
        isogeny_classes = self.catalogue.isogeny_classes()
        N_EP = 0
        _t0 = time.perf_counter()
        for ell_t in tqdm(isogeny_classes, desc=f"counting EP at ell={ell}", unit="ic", ncols=80, ascii=True):

            # Early exit: no rational ℓ-torsion can occur in this class.
            if ell_t.N_pts % ell != 0:
                continue

            if use_CN and ell_t.ordinary:
                # In class-number mode, use order counts instead of explicit curves.
                for f, o in ell_t.orders.items():

                    # TODO: ADDED FALSE HERE TO ALWAYS USE
                    if ell_t.D_K in [-3, -4] and int(f) == 1:
                        # The exceptional CM cases still need explicit curves.
                        curves = ell_t.curves_by_order.get(int(f))

                        for c in curves:
                            torsion_subgroup = TorsionSubgroup(c, ell)
                            torsion_subgroup.compute_rank(f_pi=ell_t.f_pi, use_generators=False)
                            n_orb = torsion_subgroup.count_orbits()
                            N_EP += n_orb
                    else:
                        # from apper, recall pl = 2^r-1 for ell = 2
                        m = 2 if ell > 2 else 1
                        r = 2 if ZZ(o.conductor).valuation(ell) < ZZ(ell_t.f_pi).valuation(ell) else 1
                        N_EP += o.class_number * (ell**r - 1) / m
            else:
                ## SS CURVES
                for f, curves_list in ell_t.curves_by_order.items():
                    for c in curves_list:
                        torsion_subgroup = TorsionSubgroup(c, ell)
                        torsion_subgroup.compute_rank(f_pi=ell_t.f_pi, use_generators=False)
                        n_orb = torsion_subgroup.count_orbits()
                        N_EP += n_orb

        return N_EP

    def compute_hecke(self, k, level, use_CN=False) -> Tuple[int, int, int, List[int], List[int]]:
        """Compute the trace contribution of the Hecke operator T_{level} in weight `k`."""
        from tqdm import tqdm
        import time
        hk_symbolic = Hk.construct(k)
        isogeny_classes = self.catalogue.isogeny_classes()
        T = 0
        NC = 0
        NSS = 0
        traces = []
        hk_evals = []
        vals = []
        full_r = False
        _t0 = time.perf_counter()

        for ell_t in isogeny_classes:

            if ell_t.N_pts % level != 0: # note, this holds true even if level not prime, converse nott rue tho, we cans till end up with 0 pts of order N
                continue

            traces.append(ell_t.t)

            hk = ell_t.eval_hk_mod_fx(level, hk_symbolic)
            hk_evals.append(hk)

            MAX_N = num_order_N(level, level, level)

            vl_fpi = ZZ(ell_t.f_pi).valuation(level)
            clr = Colors.BOLD if vl_fpi == 0 else Colors.WARNING

            '''if not ell_t.ordinary:
                print(
                    f"{clr}Computing EN structure for level={level}, vl_pi={ZZ(ell_t.f_pi).valuation(level)}, f_pi={ell_t.f_pi} | N_pts={ell_t.N_pts}, #E^2 mod ell = {ell_t.N_pts % level**2}, q mod level = {self.field.q % level}, MAX_N={MAX_N}{Colors.ENDC}"
                )'''

            # note, the contribution for each curve / order does not have to be integer value, for |aut| > 2 we can get n/3 or n/2 sums but these sum to integers over twists, uncomment below to see this in action
            if use_CN and ell_t.ordinary:
                has_full_r = False

                for f, o in ell_t.orders.items():
                    aut_size = 2
                    if ell_t.D_K in [-3, -4] and int(f) == 1:
                        aut_size = 6 if ell_t.D_K == -3 else 4

                    NP = num_P(level, ell_t.f_pi, f, ell_t.N_pts, self.field.q)

                    accum_val = hk * o.class_number * NP / aut_size
                    T -= accum_val
                    vals.append(accum_val)

                    '''print(f"t={ell_t.t}, hk={hk}, NC={o.class_number}, NP={NP}, contribution to count: {accum_val}")'''

                    NC += o.class_number

                    '''for f2, curves_list in ell_t.curves_by_order.items():
                        if ZZ(f2) == ZZ(f):
                            idx = 0
                            for c in curves_list:

                                # torsion_subgroup = TorsionSubgroup(c, level)
                                # NP_enum = torsion_subgroup.get_num_points_exact_order()

                                l_syl = invariants(level, ell_t.f_pi, f, ell_t.N_pts, self.field.q)

                                if not has_full_r and l_syl[0] == level and l_syl[1] == level:
                                    has_full_r = True
                                    # print(f"________________________DETECT FULL RANK {level}________________________")

                                NP_enum2 = num_order_N(level, l_syl[0], l_syl[1])

                                if idx == 0:
                                    print(
                                        f"l_sylow_structure = {l_syl}, ab invariants = {c.getSageCurve().abelian_group().invariants()}, NP_enum2={NP_enum2}, NP_enum={NP_enum}"
                                    )

                                idx += 1

                                if NP_enum != NP:
                                    print(
                                        f"\n {Colors.FAIL}COMPUTED {NP} does not match ENUMERATED {NP_enum}{Colors.ENDC}"
                                    )
                                    print(f"Details: t={ell_t.t}, D_K={ell_t.D_K}, f_pi={ell_t.f_pi}, N_pts={ell_t.N_pts}, class_number={o.class_number}, j0={c.is_j0}, j1728={c.is_j1728}")
                                    pts = c.getSageCurve().points()
                                    for P in pts:
                                        print(f"Point {P} has order {P.order()}")'''

                '''print(
                    f"{Colors.FAIL if vl_fpi > 0 and not has_full_r else (Colors.GREEN if vl_fpi > 0 and has_full_r else Colors.BOLD)}---------COMPUTED has_full_r={has_full_r}{Colors.ENDC}"
                )'''

            else:

                # TODO: we skip this for now to debug
                # if abs(ell_t.t) == self.HB:
                #    continue

                '''fx_roots = ell_t.fx_pi.roots()
                r_mod_ell_list = []
                for root, multiplicity in fx_roots:
                    r_mod_ell = Zmod(level)(root)
                    r_mod_ell_list.append(r_mod_ell)'''

                '''print(
                    f"SS ORDER, t={ell_t.t}, D_K={ell_t.D_K}, f_pi={ell_t.f_pi}, N_pts={ell_t.N_pts}, q mod level = {self.field.q % level}, K={ell_t.K}, orders={ell_t.orders}{Colors.ENDC}"
                )'''

                for f, curves_list in ell_t.curves_by_order.items():

                    #print(f"{Colors.HEADER}Processing curves of order f={f} for trace t={ell_t.t}, D_K={ell_t.D_K}, f_pi={ell_t.f_pi}, N_pts={ell_t.N_pts}, len(curves_list)={len(curves_list)}{Colors.ENDC}")
                    for c in curves_list:
                        torsion_subgroup = TorsionSubgroup(c, level)

                        NP = 0
                        r = torsion_subgroup.compute_rank(f_pi=ell_t.f_pi, use_generators=False)
                        if r is not None:
                            NP = level**r-1 #torsion_subgroup.get_num_points_exact_order()
                        else:
                            print(f"{Colors.FAIL}Failed to compute rank for curve with j={c.j}, t={c.t}, f_E={c.f_E}, f_pi={ell_t.f_pi}, N_pts={ell_t.N_pts}{Colors.ENDC}")
                        # NPSS = num_P(level, ell_t.f_pi, f, ell_t.N_pts, self.field.q)
                        # print(f"NP = {NP}, MAX_N ={MAX_N}, num_P_SS={NPSS}, f_pi={ell_t.f_pi}, f={f}")

                        accum_val = hk * NP / c.aut_size
                        T -= accum_val
                        NC += 1
                        NSS += 1
                        vals.append(accum_val)

                        '''print(
                            f"t={ell_t.t}, hk={hk}, NC={1}, NP={NPSS}, contribution to count: {accum_val}"
                        )'''

        # print(f"Computed hecke trace in {time.perf_counter() - _t0:.5f}s")
        return T, NC, NSS, traces, hk_evals, vals, full_r

    def toJSON(self) -> Dict[str, Any]:
        """Serialize the finite-field catalogue in the frontend/backend JSON format."""
        return {
            "char": int(self.field.p),
            "catalogue": self.catalogue.toJSON()
        }
