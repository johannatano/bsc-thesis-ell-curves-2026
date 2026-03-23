"""Curve objects and finite-field curve catalogues used by the project.

This module defines both lightweight CM/HCP curve records and fully geometric
elliptic-curve models over finite fields. It also provides the catalogue layer
that places these curves into trace-indexed isogeny classes supplied by the
number-field infrastructure.
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
from lib.nr_fields import *
#from utils.mod_poly import _classical_modular_polynomial


highest_ell = 0

def get_j_invariants_from_order(D: int, f: int, q: int) -> List:
    """Return the j-invariants attached to an order via its Hilbert class polynomial.
    
    Args:
        D: Discriminant of the order $D_K f^2$.
        f: Conductor of the order. Included for API clarity.
        q: Size of the finite field.
    
    Returns:
        Roots of the Hilbert class polynomial in $F_q$.
    """
    j_invs = []
    try:
        H = hilbert_class_polynomial(D)
        F_q = GF(q)
        # Find roots of H(x) in F_q
        for j in H.change_ring(F_q).roots(multiplicities=False):
            j_invs.append(j)
    except Exception as e:
        print(f"Warning: Could not compute HCP for D={D}: {e}")
    return j_invs


class Curve:
    """Common interface for elliptic curves tracked by the project.

    The class stores the arithmetic data shared by both lightweight
    CM/Hilbert-class-polynomial records and fully instantiated Sage curves.
    """
    @staticmethod
    def ABFromJ(j) -> Tuple:
        """Construct a short Weierstrass model with the given j-invariant.

        This helper assumes the usual characteristic restrictions for the
        formula in the generic case.
        """
        if((j.is_zero())):
            return (0, 1)
        elif((j - 1728).is_zero()):
            return (1, 0)
        else:
            return ((3*j) / (1728 - j), (2*j) / (1728 - j)) #formula for A,B with given j (char not 2,3 implicit)
        
    def __init__(self, Fq: FqData, j, aut_grp: Optional[Dict] = None, t: Optional[int] = None, f_E: Optional[int] = None) -> None:
        self.field: FqData = Fq
        self.p: int = Fq.p
        self.n: int = Fq.n
        self.q: int = Fq.q
        self.F = Fq.F
        
        self.A = None
        self.B = None
        
        self.j = self.F(j)
        self.j_invariant = self._inv(self.j)
        self.j_invariant_flat = self._inv_flatten(self.j_invariant[1])
        
        self.is_j0 = (self.j.is_zero())
        self.is_j1728 = (self.j - 1728).is_zero()
        
        self.t = t
        if t is not None:
            self.is_supersingular = self.t % self.p == 0
        
        self.f_E = f_E
        
        self.N_pts = self.q + 1 - self.t if t is not None else None
        self.is_generic = not self.is_j0 and not self.is_j1728
        
        self.aut_grp = aut_grp
        self.aut_size = aut_grp["size"] if aut_grp is not None else 2
        
        self.E = None  # placeholder for Sage elliptic curve object, to be initialized when coefficients are known
             
    def getCoefficients(self):
        """Recover Weierstrass coefficients compatible with the stored trace.

        For generic curves this searches the twist family determined by the
        automorphism group until it finds a model with Frobenius trace `self.t`.
        """
        if self.A is not None and self.B is not None:
            return (self.A, self.B)
        
        #print(f"{Colors.WARNING}Warning: Coefficients A,B not set for curve with j={self.j}, t={self.t}. Attempting to compute from j-invariant...{Colors.ENDC}")
        # Search the relevant twist family until the Frobenius trace matches.
        if self.aut_grp is None:
            print(f"{Colors.FAIL}Error: Cannot compute twists without automorphism group info{Colors.ENDC}")
            return (self.A, self.B)
        
        if self.is_j0:
            e = (0, 1)
        elif self.is_j1728:
            e = (1, 0)
        else:
            e = (2, 3)
        
        coset_reps = [self.field.g ** k for k in range(0, self.aut_size)]
        for u in coset_reps:
            A, B = self.ABFromJ(self.j)
            E = EllipticCurve(self.F, [u**e[0] *A, u**e[1] *B])
            t = E.trace_of_frobenius()
            if t == self.t:
                self.A, self.B = u**e[0] *A, u**e[1] *B
                #print(f"{Colors.GREEN}Successfully computed coefficients A,B from j-invariant for curve with j={self.j}, coeff t={t}, my t={self.t}{Colors.ENDC}")
                break
    
        return (self.A, self.B)
    
    def getSageCurve(self):
        """Return the cached Sage `EllipticCurve`, creating it on demand."""
        if self.E is not None:
            return self.E
        A, B = self.getCoefficients()
        self.E = EllipticCurve(self.F, [A, B])
        return self.E
     
    def _inv(self, el) -> Tuple[int, Tuple]:
        poly = el.minpoly()
        #print(f"Element {el} has minimal polynomial {poly} of degree {poly.degree()} with coefficients {poly.list()}")
        return poly.degree(), tuple(int(c) for c in poly.list())
    
    def _inv_flatten(self, inv: Tuple) -> int:
        index = 0
        for c in reversed(inv):
            index = index * self.field.p + c
        return index
    
    def _rank_by_group_structure(self, ell) -> int:
        if self.E is None:
            print(f"{Colors.WARNING}Warning: Curve not initialized with Weierstrass form, cannot compute rank by division polynomial{Colors.ENDC}")
            return self._rank_by_modular_poly(ell)
        invariants = self.E.abelian_group().invariants()
        r = 0
        for inv in invariants:
            if inv % ell == 0:
                r += 1
        return r
    
    def _rank_by_div_poly(self, ell):
        if self.E is None:
            print(f"{Colors.WARNING}Warning: Curve not initialized with Weierstrass form, cannot compute rank by division polynomial{Colors.ENDC}")
            return self._rank_by_modular_poly(ell)
        psi_n = self.E.division_polynomial(ell)
        if ell > 50:
            print(f"{Colors.WARNING}Warning: Division polynomial for ell={ell} may be expensive to compute, consider using modular polynomial method instead{Colors.ENDC}")
        n_roots = sum(m for _, m in psi_n.roots(multiplicities=True))
        return 2 if n_roots > 2 else (1 if n_roots > 0 else 0)
            
    def _rank_by_modular_poly(self, ell):
        # Modular polynomials are slower per call, but avoid constructing full
        # group structure data and are easier to cache across many curves.
        x = polygen(self.field.F)
        X, Y = polygens(self.field.F, 'X,Y')
        phi = classical_modular_polynomial(ell)(X, Y)
        phi_j = phi([x, self.j]) 
        #phi_j = classical_modular_polynomial(ell, self.j)
        n_roots = sum(m for _, m in phi_j.roots(multiplicities=True))
        return 2 if n_roots > 2 else (1 if n_roots > 0 else 0)
    
    def _above_floor(self, ell: int) -> bool:
        """Heuristic/test for whether the curve lies above the volcano floor."""
        global highest_ell
        if ell > highest_ell:
            highest_ell = ell
        j = self.j
        # Special j-values correspond to maximal-order cases in this model.
        if j.is_zero() or (j - 1728).is_zero():
            return True
        from utils.common import Config
        if Config.rank_method == "auto":
            return self._rank_by_div_poly(ell) == 2 if ell < 13 else self._rank_by_modular_poly(ell) == 2
        elif Config.rank_method == "mod_poly":
            return self._rank_by_modular_poly(ell) == 2
        elif Config.rank_method == "invariants":
            return self._rank_by_group_structure(ell) == 2
        return self._rank_by_div_poly(ell) == 2
    
    def height_above_floor(self, ell, e, use_true_height=False):
        """Compute height above floor in ℓ-isogeny volcano.
        
        Args:
            ell: Prime defining the isogeny degree
            e: ℓ-adic valuation of f_π
            use_true_height: If True, compute exact height via BFS in isogeny graph
        
        Returns:
            Height above volcano floor (0 means on floor, e means at surface)
        """
        if not use_true_height:
            return e if self._above_floor(ell) else 0
        
        e = ZZ(e)
        j = self.j  
        if not e:
            return ZZ.zero()
        if self.is_j0 or self.is_j1728:
            return e
        x = polygen(self.F)
        X, Y = polygens(self.F, 'X,Y')
        phi = classical_modular_polynomial(ell)(X, Y)
        phi_j = phi([x, j]) 
        j_roots_mult = phi_j.roots(multiplicities=True)
        roots = []
        for r, m in j_roots_mult:
            for _ in range(m):
                roots.append(r)
                
        nj1 = len(roots)
        on_floor = self.E.two_torsion_rank() < 2 if ell == 2 else nj1 <= ell
        if on_floor:
            return ZZ.zero()
        if e == 1 or nj1 != ell + 1:
            return e
        if nj1 < 3:
            return ZZ.zero()
        current_nodes = roots[:3] 
        previous_nodes = [j, j, j]
        h = ZZ.one()
        while True:
            for i in range(3):
                val = current_nodes[i]
                prev = previous_nodes[i]
                next_poly = phi([x, val]) // (x - prev)
                next_roots = next_poly.roots(multiplicities=False)
                if not next_roots:
                    return h
                previous_nodes[i] = val
                current_nodes[i] = next_roots[0] # Take the first valid neighbor
            h += 1
            if h > e: 
                return h
            
    def compute_conductor(self, f_pi: int, use_true_height: bool = False) -> int:
        """Compute the conductor of the endomorphism ring of the curve.

        The computation factors the Frobenius conductor `f_pi` and removes the
        contribution coming from the curve's height in each ℓ-volcano.
        """
        if self.f_E is not None:
            return self.f_E
        self.f_E = 1
        if self.is_j0 or self.is_j1728:
            return self.f_E  # Special j-values stay at maximal order here.
        if not self.is_supersingular:
            for _l, e in f_pi.factor():
                if not use_true_height and (self.q-1) % _l != 0:
                    self.f_E *= _l**e
                    # No relevant rational torsion for this prime, so the height is 0.
                    continue
                h = self.height_above_floor(_l, e, use_true_height=use_true_height)
                self.f_E *= _l**(e-h)
        return self.f_E
    
    def toID(self) -> str:
        """Override in subclasses"""
        raise NotImplementedError("Subclasses must implement toID()")
    
    def toJSON(self, include_points: bool = False) -> Dict[str, Any]:
        """Override in subclasses"""
        raise NotImplementedError("Subclasses must implement toJSON()")

        
class NFCurve(Curve):
    """Lightweight curve representation for number field enumeration.
    
    Stores only the j-invariant and minimal metadata, without constructing a
    full Weierstrass model. This is the representation used during CM/HCP-based
    enumeration before concrete curve models are needed.
    """
    
    def __init__(self, Fq: FqData, j, aut_grp: Optional[Dict] = None, t: Optional[int] = None, f_E: Optional[int] = None) -> None:
        super().__init__(Fq, j, aut_grp, t, f_E)
        self.ID: str = self.toID()
        
        


    def toID(self) -> str:
        return ''.join(str(self.j).split())
            
    def toJSON(self, include_points: bool = False) -> Dict[str, Any]:
        return {
            "ID": self.ID,
            "j": element_to_tuple(self.j),
            "j_minpoly": self.j_invariant[1],
        }

class GeometricCurve(Curve):
    """Full geometric curve with Weierstrass form.
    
    Creates a Sage `EllipticCurve` object and supports the computations needed
    later in the pipeline, such as twists, torsion data, and volcano placement.
    """
    

        
    def __init__(self, Fq: FqData, j, aut_grp: Optional[Dict] = None, A = None, B = None, t: Optional[int] = None, f_E: Optional[int] = None) -> None:
        super().__init__(Fq, j, aut_grp, t, f_E)
        self.g = Fq.g
        
        self.torsion = {}
        # If not initialized with A,B, compute from j
        if(A is None or B is None):
            A, B = GeometricCurve.ABFromJ(self.j)
        self.A = self.F(A)
        self.B = self.F(B)
        self._create()
        self.ID = self.toID() # short hash ID for quick reference, not guaranteed unique but should be good enough for our scale
 
    def _create(self):
        # Materialize the actual Sage curve once the coefficients are fixed.
        self.E = EllipticCurve(GF(self.q), [self.A, self.B])
        if self.t is None:
            self.t = self.E.trace_of_frobenius() # TODO, reuse t = -t twist? for j = 0,1728 ??
        self.N_pts = self.q + 1 - self.t
        ss = self.E.is_supersingular()
        self.is_supersingular = self.t % self.p == 0 # self.E.is_supersingular() --> bit nmroe involved but bcs avoids calc cardinality, so just use t method here isntead
        self.is_generic = not self.is_j0 and not self.is_j1728
        if(ss != self.is_supersingular):
            print(f"{Colors.WARNING}Discrepancy in supersingularity check for j={self.j}: by t={self.is_supersingular}, by is_supersingular()={ss}{Colors.ENDC}")

    def toID(self) -> str:
        import hashlib
        # Short stable identifier for serialized graph/JSON references.
        a_tuple = element_to_tuple(self.A)
        b_tuple = element_to_tuple(self.B)
        hash_input = f"{a_tuple}{b_tuple}".encode('utf-8')
        return hashlib.sha256(hash_input).hexdigest()[:8]
        
    def compute_twists(self) -> List['GeometricCurve']:
        """Compute all twists of this curve with the same j-invariant.
        
        Returns list including self and all non-trivial twists.
        For special j-invariants (0, 1728), computes sextic/quartic twists.
        Populates self.twists on every curve in the family.
        """
        if self.aut_grp is None:
            print(f"{Colors.FAIL}Error: Cannot compute twists without automorphism group info{Colors.ENDC}")
            return []
        
        if self.is_j0:
            e = (0, 1)
        elif self.is_j1728:
            e = (1, 0)
        else:
            e = (2, 3)
        
        coset_reps = [self.g ** k for k in range(1, self.aut_size)]
        t = self.t * -1 if not self.is_j0 and not self.is_j1728 else None
        result = [self] + [GeometricCurve(self.field, self.j, self.aut_grp, u**e[0] * self.A, u**e[1] * self.B, t=t, f_E=self.f_E) for u in coset_reps]
        return result
        
    def weierstrass_polynomial(self):
        R = PolynomialRing(self.F, 'x')
        x = R.gen()
        return x**3 + self.A * x + self.B
    
    def invariants(self) -> Tuple[int, int]:
        inv = self.E.abelian_group().invariants()
        return (int(inv[0]), int(inv[1]) if len(inv) > 1 else 0)
    
    #IMPORTED FROM SAGE
    def endomorphism_order(self):
        pi = self.E.frobenius()
        if pi in ZZ:
            raise NotImplementedError('the rank-4 case is not supported yet')
        O = self.E.frobenius_order()
        f0 = O.conductor()
        f = 1
        for l,e in f0.factor():
            h = self.height_above_floor(l, e)
            f *= l**(e-h)
        K = O.number_field()
        return K.order_of_conductor(f)
            
    def toJSON(self, include_points: bool = False) -> Dict[str, Any]:
        return {
            "ID": self.ID,
            "j": element_to_tuple(self.j),
            "A": element_to_tuple(self.A),
            "B": element_to_tuple(self.B),
            "inv": self.invariants()
        }

    
class EllFiniteFieldCatalogue:
    """In-memory catalogue of curves over a fixed finite field $F_{p^n}$.

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
