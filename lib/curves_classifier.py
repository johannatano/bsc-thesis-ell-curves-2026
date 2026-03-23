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
#from utils.mod_poly import _classical_modular_polynomial

highest_ell = 0

# =============================================================================
# Curve enumeration and classification
# =============================================================================

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
    
    def enumerate_curves(self, use_HCP: bool = False, use_CN: bool = False, add_SS: bool = True) -> None:
        """Enumerate curves over the current finite field.

        `use_HCP` switches to the CM/Hilbert-class-polynomial pipeline.
        `use_CN` keeps only class-number counts where possible.
        `add_SS` manually inserts supersingular cases that are not produced by
        the HCP enumeration.
        """
        import time
        from tqdm import tqdm
        if use_HCP or use_CN:
            # Supersingular j-invariants are handled separately: the HCP pass is
            # aimed at ordinary CM data and does not cover these cases by itself.
            if add_SS:
                if self.field.p % 3 == 2:
                    print(f"{Colors.FAIL}j-invariant 0 is SS{Colors.ENDC}")
                    self.add_ss_curve_by_j(self.field.F(0))
                    
                if self.field.p % 4 == 3:
                    print(f"{Colors.FAIL}j-invariant 1728 is SS{Colors.ENDC}")
                    self.add_ss_curve_by_j(self.field.F(1728))
                
                SS_poly = supersingular_j_polynomial(self.field.p)
                for r in SS_poly.roots(multiplicities=False):
                    print(f"{Colors.HEADER}Found supersingular j-invariant root in F_{self.field.q}: j={r}, is_1728={(r-1728).is_zero()}{Colors.ENDC}")
                    self.add_ss_curve_by_j(self.field.F(r))  # we know these are SS, so t=0, and f_E=1 since they have maximal endomorphism ring
            # For each order, recover its j-invariants and attach them to all
            # trace classes compatible with that order.
            nf_list = list(self.catalogue.NFC.data.values())
            _t0 = time.perf_counter()
            for nf in tqdm(nf_list, desc=f"HCP F_{self.field.q}", unit="nf", ncols=80, ascii=True):
                D_K = nf.discriminant
                # Orders are grouped by quadratic field and extension degree.
                orders = nf.getOrders(self.field.n)
                for order in orders:
                    f = order.conductor
                    if not use_CN or (D_K in [-3, -4] and int(f) == 1):
                        j_invs = get_j_invariants_from_order(D_K * f**2, f, self.field.q)
                        for j_inv in j_invs:
                            order.add_j_invariant(j_inv)
                            # The same order may contribute curves to several
                            # trace classes with the same CM field data.
                            for t in order.traces:
                                self.add_nf_curve(j_inv, t=t, f_E=f)
                    else:
                        # In class-number mode we only track multiplicities.
                        self.catalogue.size += order.class_number*len(order.traces)
                        
        else:
            _t0 = time.perf_counter()
            precompute_conductor = True
            for j in tqdm(self.field.F, total=self.field.q, desc=f"F_{self.field.q}", unit="j", ncols=80, ascii=True):
                self.add_curves_by_j(j, pre_compute_conductor=precompute_conductor)
        
        global highest_ell
        # Standard count of isomorphism classes over F_q in characteristic > 3.
        NE = 2*(self.field.q -2) + gcd(4, self.field.q-1) + gcd(6, self.field.q-1)
        if self.catalogue.size != NE:
            print(f"{Colors.FAIL}Warning: total number of curves in catalogue ({self.catalogue.size}) does not match expected number from formula ({NE}), there may be duplicates or missing curves{Colors.ENDC}")
        else:
            print(f"{Colors.GREEN}Successfully enumerated curves, total size of catalogue: {self.catalogue.size}, highest_ell={highest_ell}{Colors.ENDC}")
        
        print(f"Enumeration done in {time.perf_counter() - _t0:.2f}s")
                
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
        for E_t in twists:
            self.catalogue.add(E_t)

    def add_ss_curve_by_j(self, j) -> None:
        """Insert the supersingular twist family attached to `j`."""
        aut_grp = self._get_aut_group_for_j(j)
        E = GeometricCurve(self.field, j, aut_grp=aut_grp, t=None, f_E=1)
        ell_t = self.catalogue.get_isogeny_class(E.t)
        twists = E.compute_twists()
        for E_t in twists:
            self.catalogue.add(E_t)
               
    def add_nf_curve(self, j, t: int, f_E: int) -> None:
        """Insert a lightweight CM-derived curve record into the catalogue."""
        aut_grp = self._get_aut_group_for_j(j)
        self.catalogue.add(NFCurve(self.field, j, aut_grp=aut_grp, t=t, f_E=f_E))
    
    def compute_volcano(self, ell:int = -1, edges: bool = False) -> None:
        """Build ℓ-isogeny volcano data for each isogeny class."""
        from tqdm import tqdm
        import time
        primes = [ell] if ell != -1 else list(primerange(2, max_ell_from_HB(self.field.q)+1))
        isogeny_classes = self.catalogue.isogeny_classes()
        _t0 = time.perf_counter()
        for ell_t in tqdm(isogeny_classes, desc="computing torsion", unit="ic", ncols=80, ascii=True):
            for ell in primes:
                if ell == self.field.p:
                    continue
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
                    
                    if ell_t.D_K in [-3, -4] and int(f) == 1:
                        # The exceptional CM cases still need explicit curves.
                        curves = ell_t.curves_by_order.get(int(f))
                        for c in curves:
                            torsion_subgroup = TorsionSubgroup(c, ell)
                            torsion_subgroup.compute_rank(f_pi=ell_t.f_pi, use_generators=False)
                            N_EP += torsion_subgroup.count_orbits()
                    else:    
                        r = 2 if ZZ(o.conductor).valuation(ell) < ZZ(ell_t.f_pi).valuation(ell) else 1   
                        N_EP += o.class_number*(ell**r - 1) // 2
            else:
                for f, curves_list in ell_t.curves_by_order.items():
                    for c in curves_list:
                        torsion_subgroup = TorsionSubgroup(c, ell)
                        torsion_subgroup.compute_rank(f_pi=ell_t.f_pi, use_generators=False)
                        N_EP += torsion_subgroup.count_orbits()
        print(f"Computed torsion in {time.perf_counter() - _t0:.2f}s")
        return N_EP
        
    def compute_hecke(self, k, level) -> int:
        """Compute the trace contribution of the Hecke operator $T_{level}$ in weight `k`."""
        from tqdm import tqdm
        import time
        hk_symbolic = Hk.construct(k)
        isogeny_classes = self.catalogue.isogeny_classes()
        T = 0
        _t0 = time.perf_counter()
        for ell_t in tqdm(isogeny_classes, desc="computing torsion", unit="ic", ncols=80, ascii=True):
            hk = ell_t.eval_hk_mod_fx(level, hk_symbolic)
            v = ell_t.volcanoes.get(level)
            if not v.hasStructure():
                continue
            for f, curves_list in ell_t.curves_by_order.items():
                r = 2 if ZZ(f).valuation(level) < ZZ(ell_t.f_pi).valuation(level) else 1
                for c in curves_list:
                    T -= hk * (level**r-1) // c.aut_size
        print(f"Computed hecke trace in {time.perf_counter() - _t0:.2f}s")
        return T
    
    def toJSON(self) -> Dict[str, Any]:
        """Serialize the finite-field catalogue in the frontend/backend JSON format."""
        return {
            "char": int(self.field.p),
            "catalogue": self.catalogue.toJSON()
        }