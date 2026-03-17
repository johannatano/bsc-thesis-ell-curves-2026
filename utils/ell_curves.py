"""
Elliptic curve enumeration and classification over finite fields.

This module provides tools for:
- Enumerating elliptic curves over F_q using Hilbert class polynomials (CM theory)
- Organizing curves by isogeny classes and endomorphism orders
- Computing torsion subgroups and isogeny volcano structures
- Separating t-independent data (orders, j-invariants) from t-dependent data (isogeny classes)
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
from utils.ell_nr_field import *
#from utils.mod_poly import _classical_modular_polynomial


highest_ell = 0

def get_j_invariants_from_order(D: int, f: int, q: int) -> List:
    """Compute j-invariants for an endomorphism order using Hilbert class polynomial.
    
    Args:
        D: Discriminant of the order (D_K * f^2)
        f: Conductor of the order
        q: Size of the finite field
    
    Returns:
        List of j-invariants (roots of HCP in F_q)
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


class EPOrbit:
    def __init__(self, aut_size: int) -> None:
        self.aut_size: int = aut_size
       
class Point:
    def __init__(self, P) -> None:
        self.P = P
        self._orbit: Optional[List] = None
        self._stab_size: int = 1

    def phi(self, u):
        if self.P.is_zero(): return self.P
        x, y = self.P.xy()
        return (u^2 * x, u^3 * y)

    def orbit(self, G) -> List:
        isXZero = False
        if self._orbit is None:
            self._orbit = [(g**2 * self.P.x(), g**3 * self.P.y()) for g in G]
            unique_orbit = []
            seen = set()
            for pt in self._orbit:
                if pt not in seen:
                    unique_orbit.append(pt)
                    seen.add(pt)
            self._orbit = unique_orbit
        return self._orbit
    
    def point(self):
        return self.P
    
    def xy(self) -> Tuple:
        return self.P.xy()
    
    def toTuple(self) -> Tuple[Tuple, Tuple]:
        return (element_to_tuple(self.P.x()),element_to_tuple(self.P.y()))
    
class TorsionSubgroup:
    """Computes and stores ℓ-torsion subgroup data for an elliptic curve."""
    
    def __init__(self, curve, l: int) -> None:
        self.curve = curve
        self.l: int = l
        self.gens: Optional[List] = None
        self.orbits: List[Point] = []
        self.rank: int = 0
        self.n_orbits: int = 0
        self.points: List[Point] = []
        self.special_aut6_orbit: bool = False

    def compute_rank(self, use_generators: bool = False) -> None:
        if self.curve.N_pts % self.l != 0:
            self.rank = 0
            return
        if use_generators:
            self.gens = self._get_generators()
            self.rank = len(self.gens)
        else:
            self.rank = self._compute_rank()

    def generate_orbits(self):
        if self.rank == 2:
            for i in range(self.l):
                for j in range(self.l):
                    self._add(i * self.gens[0] + j * self.gens[1])
        else:
            for i in range(self.l):
                self._add(i * self.gens[0])
                
    def count_orbits(self) -> int:
        return self._count_orbits_level_2() if self.l == 2 else (self._count_orbits_level_3() if self.l == 3 else self._count_orbits_general())
    
    def _count_orbits_level_2(self) -> int:
        if(self.curve.aut_size == 2):
            return 2**self.rank - 1
        elif(self.curve.aut_size == 4):
            return 2**(self.rank -1)
        else:
            return 1
        
    def _count_orbits_level_3(self) -> int:
        fixed = 3**self.rank - 1
        if(self.curve.aut_size == 6):
            fixed += 2*(1+self.curve.B.is_square())
        return (fixed // self.curve.aut_size)
    
    def _count_orbits_general(self) -> int:
        has_fixed_by_6 = self._get_fixed_points_count_aut6()
        if self.curve.aut_size == 6 and has_fixed_by_6 > 0:
            print(f"{Colors.FAIL}Has fixed points from aut size 6 BUT SHOULD NOT HAVE{Colors.ENDC}, for self.aut_size={self.curve.aut_size}, l={self.l}")
        return (self.l**self.rank-1) // self.curve.aut_size
    
    def _get_num_fixed(self) -> int:
        if(self.l == 2):
            if(self.curve.aut_size < 4):
                return self.l**self.rank*self.curve.aut_size
            elif(self.curve.aut_size == 4):
                div_poly = self.curve.E.division_polynomial(self.l)
                if(div_poly(0) == 0):
                    return (self.l**self.rank-2)*2 + self.curve.aut_size*2 # 0, (0,0) always fixed + all others fixed by inverse and id
                else:
                    return (self.l**self.rank-1)*2 + self.curve.aut_size # 0, (0,0) always fixed + all nonzero fixed by inverse and id
            else:
                return ((self.l**self.rank-1)*2 + self.curve.aut_size) # 0 always fixed + all nonzero fixed by inverse and id
        fixed = self.curve.aut_size # id fixed by all
        fixed += ( self.l**self.rank - 1) # all nonzero points fixed by id
    
        if(self.curve.aut_size == 6):
            n_extra = self._get_fixed_points_count_aut6()
            fixed += n_extra # add extra fixed from special case
            self.special_aut6_orbit = True if n_extra > 0 else False
        return fixed

    def _add(self, ell_P) -> None:
        if ell_P.is_zero():
            return
        P = Point(ell_P)
        self.points.append(P)
        if self._check_unique_orbit(P):
            self.orbits.append(P)
            self.n_orbits += 1
            
    def _get_fixed_points_count_aut6(self) -> int:
        #if(self.curve.aut_size != 6):
        #    return 0
        # we need to get all l-torsion points with x = 0, hence, check wheter x = 0 is a root of the division polynomial
        div_poly = self.curve.E.division_polynomial(self.l)
        return 4 if div_poly(0) == 0 else 0
    
    def _check_unique_orbit(self, torsion_point: Point) -> bool:
        for S in self.orbits:
            for P in S.orbit(self.curve.aut_grp):
                if torsion_point.point().xy() == P:
                    return False
        return True
    
    def _rank_2_test(self) -> bool:
        if(self.curve.N_pts % self.l**2 != 0): # necessary: group order multiple of l^2
           return False
        #if((self.curve.t - 2) % self.l != 0): # necessary: alpha + beta = t for a double root x = 1 gives t = 2 mod l
        #    return False
        if((self.curve.q - 1) % self.l != 0): # necessary: alpha*beta = q for a double root x = 1 gives q = 1 mod l
            return False
        return True
    
    def _rank_by_group_structure(self) -> int:
        invariants = self.curve.E.abelian_group().invariants()
        r = 0
        for inv in invariants:
            if inv % self.l == 0:
                r += 1
        return r
    
    def _rank_by_modular_poly(self) -> int:
        mod_poly = classical_modular_polynomial(self.l, self.curve.j)
        factors = mod_poly.factor()
        nr_linear_factors = 0
        for (f, m) in factors:
            if f.degree() == 1:
                nr_linear_factors += m
        return 2 if nr_linear_factors > 1 else 1
    
    def _rank_by_enum_points(self) -> int:
        torsion_pts = [P for P in self.curve.E.points() if not P.is_zero() and (self.l * P).is_zero()]
        return 1 if len(torsion_pts) == (self.l-1) else 2  # rough estimate, may overcount if not full rank
    
    def _two_torsion_rank(self) -> int:
        if (self.curve.N_pts % 2 != 0):
            return 0
        w_poly = self.curve.weierstrass_polynomial()
        splits = all(f.degree() == 1 for f, m in w_poly.factor())
        return 2 if splits else 1
        
    def _compute_rank(self) -> int:
        if self.curve.N_pts % self.l != 0:
            return 0
        if not self._rank_2_test():
            return 1
        if self.curve.is_j0 or self.curve.is_j1728 or self.curve.is_supersingular:
            rk2 = self._rank_by_group_structure()
            return rk2
        return self._rank_by_modular_poly()

    def _get_generators(self) -> List:
        """
        Returns a 2-tuple of torsion generators as curve points.
        Guarantees (P, Q) even if subgroup rank < 2.
        Missing entries are filled with None.
        """
        # TAKE FROM SAGE: torsion_basis SOURCE CODE, MODIFIED TO ALLOW RANK 1
        #return tuple(P.element() for P in self.E.abelian_group().torsion_subgroup(n).gens())
        return [P.element() for P in self.curve.E.abelian_group().torsion_subgroup(self.l).gens()]
        
    def toJSON(self) -> Dict[str, int]:
        return {"rank": self.rank}


class Curve:
    """Base class for elliptic curves containing common data fields"""
    
    def __init__(self, Fq: FqData, j, aut_grp: Optional[Dict] = None, t: Optional[int] = None, f_E: Optional[int] = None) -> None:
        self.field: FqData = Fq
        self.p: int = Fq.p
        self.n: int = Fq.n
        self.q: int = Fq.q
        self.F = Fq.F
        
        self.j = self.F(j)
        self.j_invariant = self._inv(self.j)
        self.j_invariant_flat = self._inv_flatten(self.j_invariant[1])
        
        self.is_j0 = (self.j.is_zero())
        self.is_j1728 = (self.j - 1728).is_zero()
        
        self.t = t
        self.f_E = f_E
        
        self.N_pts = self.q + 1 - self.t if t is not None else None
        self.is_generic = not self.is_j0 and not self.is_j1728
        
        self.aut_grp = aut_grp
        self.aut_size = aut_grp["size"] if aut_grp is not None else 2
    
    def _inv(self, el) -> Tuple[int, Tuple]:
        poly = el.minpoly()
        print(f"Element {el} has minimal polynomial {poly} of degree {poly.degree()} with coefficients {poly.list()}")
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
    
    
    '''def _rank_by_modular_poly(self, ell):
        # more cache firendly
        x = polygen(self.field.F)
        X, Y = polygens(self.field.F, 'X,Y')
        phi = classical_modular_polynomial(ell)(X, Y)
        phi_j = phi([x, self.j]) 
        #phi_j = classical_modular_polynomial(ell, self.j)
        n_roots = len(phi_j.roots(multiplicities=True))
        return 2 if n_roots > 2 else (1 if n_roots > 0 else 0)
    
    def _above_floor(self, ell: int) -> bool:
        return self._rank_by_modular_poly(ell) == 2
        
    def compute_prime_conductor(self, f_pi: int) -> int:
        self.f_E = 1
        if not self.is_supersingular: # all ss curves have maximal order ie fe = 1
            for _l, e in f_pi.factor(): # since we require ell | f_pi, only need to consider the prime factors of f_pi
                if (self.q-1) % _l != 0:
                    self.f_E *= _l**e # volcano has no height at this l, set ladic height at max
                    continue
                h = e if self._above_floor(_l) else 0
                self.f_E *= _l**(e-h)
        return self.f_E    
        '''
            
            
    def _rank_by_modular_poly(self, ell):
        # more cache firendly
        x = polygen(self.field.F)
        X, Y = polygens(self.field.F, 'X,Y')
        phi = classical_modular_polynomial(ell)(X, Y)
        phi_j = phi([x, self.j]) 
        #phi_j = classical_modular_polynomial(ell, self.j)
        n_roots = sum(m for _, m in phi_j.roots(multiplicities=True))
        return 2 if n_roots > 2 else (1 if n_roots > 0 else 0)
    
    def _above_floor(self, ell: int) -> bool:
        global highest_ell
        if ell > highest_ell:
            highest_ell = ell
        j = self.j
        # if the max order has height > 1 then all curves of special j lies at surface, ie max order
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
        # else, we simply factor the modular polynomial and check the multiplicity of the root corresponding to j, which gives us the number of neighbors in the isogeny graph, which is > 2 iff we are above the floor
        #phi_j = classical_modular_polynomial(ell)(X, Y)
        #j_roots_mult = phi_j.roots(multiplicities=True)
        #return True if len(j_roots_mult) > 2 else False
    
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
        if self.f_E is not None:
            return self.f_E
        self.f_E = 1
        if self.is_j0 or self.is_j1728:
            return self.f_E # all curves with special j have maximal order, so f_E = 1
        if not self.is_supersingular:
            for _l, e in f_pi.factor():
                if not use_true_height and (self.q-1) % _l != 0:
                    self.f_E *= _l**e
                    #NO TORSION FOR THIS PRIME; SO HEIGHT IS AUTOMATICALLY 0
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
    
    Stores only j-invariant and metadata, without computing Weierstrass form.
    Used when enumerating via Hilbert class polynomials.
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
    
    Creates a Sage EllipticCurve object and computes all geometric properties.
    Used for direct enumeration and when computing torsion/isogeny structure.
    """
    
    @staticmethod
    def ABFromJ(j) -> Tuple:
        if((j.is_zero())):
            return (0, 1)
        elif((j - 1728).is_zero()):
            return (1, 0)
        else:
            return ((3*j) / (1728 - j), (2*j) / (1728 - j)) #formula for A,B with given j (char not 2,3 implicit)
        
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
        # init a real E curve object in sage
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
        # Create a short hash from A and B coefficients
        a_tuple = element_to_tuple(self.A)
        b_tuple = element_to_tuple(self.B)
        hash_input = f"{a_tuple}{b_tuple}".encode('utf-8')
        return hashlib.sha256(hash_input).hexdigest()[:8]
        
    ''' PUBLIC METHODS '''
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


    def addTorsionSubgroup(self, l: int, subgroup: TorsionSubgroup) -> None:
        self.torsion[str(l)] = subgroup.toJSON()
        
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
            ell_t = self.NFC.create_isogeny_class(t, n=self.n)
        return ell_t
        
    def add(self, curve: Curve) -> None:
        t = curve.t
        ell_t = self.get_isogeny_class(t)
        ell_t.add_curve(curve)
        self.size += 1
        
    def isogeny_classes(self) -> List:
        return self.NFC.get_isogeny_classes_by_n(self.n)

    def toJSON(self) -> Dict[str, List]:
        return {
            "number_fields": [{
                "D": int(nf_info.discriminant),
                "tree": [tree.toJSON(include_curves=True) for tree in nf_info.tree]
            } for dk, nf_info in self.NFC.data.items() ]
        }

# =============================================================================
# Curve enumeration and classification
# =============================================================================

class CurvesClassifier_Fq:
    """Main classifier for enumerating and organizing elliptic curves over F_q.
    
    Supports two enumeration methods:
    1. HPC method: Uses Hilbert class polynomials and CM theory (fast, lightweight)
    2. Direct method: Enumerates all j-invariants and computes twists (slower, full geometric data)
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
        if j.is_zero() and self.nth_roots_unity[2] is not None:
            return self.nth_roots_unity[2]
        elif (j - 1728).is_zero() and self.nth_roots_unity[1] is not None:
            return self.nth_roots_unity[1]
        else:
            return self.nth_roots_unity[0]
    
    def enumerate_curves(self, use_HPC: bool = False, add_SS: bool = True) -> None:
        import time
        from tqdm import tqdm
        
        '''z4 = self.field.F.gen()
        #j_test = z3**2 + 4*z3 + 1  # field element in F_{5^3} #83
        j_test = 2*z4**3 + z4**2 + 4*z4 + 3 #29
        j_test = 3*z4**3 + 3*z4**2 + z4 + 3 #97
        ell = 71
        import time
        E = EllipticCurve(GF(self.field.q), j=j_test)
        _t = time.perf_counter()
        psi_n = E.division_polynomial(ell)
        pts = psi_n.roots(multiplicities=False)
        print(f"div_poly  ell={ell} q={self.field.q}: {(time.perf_counter()-_t)*1000:.2f} ms  ({len(pts)} roots)")
        _t = time.perf_counter()
        R = parent(j_test)['Y']
        Y = R.gen()
        #phi_l = _classical_modular_polynomial(ell, j_test)
        phi_l = classical_modular_polynomial(ell)(j_test, Y)
        pts_2 = phi_l.roots(multiplicities=True)'''
        
        '''print(floor(log(self.field.q)))
        print(f"mod_poly  ell={ell} q={self.field.q}: {(time.perf_counter()-_t)*1000:.2f} ms  ({sum(m for _,m in pts_2)} roots)")
        #return
        print(f"O({ell**2})+O({floor((ell+1)*log(self.field.q))}) mod_poly vs O({ell})+O({floor(ell**2*log(self.field.q))}) div_poly for ell={ell}, q={self.field.q}, number of roots")
        return
        print(f"Precomputed 59-division polynomial roots for later use in rank tests, number of roots: {len(pts_2)}")
        print(f"{Colors.HEADER}Enumerating Curves use_HPC={use_HPC} {Colors.ENDC}")'''
        if use_HPC:
            # first we have to manually add the curves which has supersingular j invariants since these are not included in HPC
            # check whether we are in ss case for special j-invariants, if so add them manually since these are not included in HPC enumeration
            
            if add_SS:
                if self.field.p % 3 == 2:
                    print(f"{Colors.FAIL}j-invariant 0 is SS{Colors.ENDC}")
                    for r in range(0, self.nth_roots_unity[2]['size']):
                        print(f"  Adding twist for j=0 with root of unity {r}")
                        self.add_nf_curve(self.field.F(0), t=0, f_E=1)
                    
                if self.field.p % 4 == 3:
                    print(f"{Colors.FAIL}j-invariant 1728 is SS{Colors.ENDC}")
                    for r in range(0, self.nth_roots_unity[1]['size']):
                        print(f"  Adding twist for j=1728 with root of unity {r}")
                        self.add_nf_curve(self.field.F(1728), t=0, f_E=1)
                
                SS_poly = supersingular_j_polynomial(self.field.p)
                for r in SS_poly.roots(multiplicities=False):
                    print(f"{Colors.HEADER}Found supersingular j-invariant root in F_{self.field.q}: j={r}, is_1728={(r-1728).is_zero()}{Colors.ENDC}")
                    self.add_nf_curve(self.field.F(r), t=0, f_E=1)  # we know these are SS, so t=0, and f_E=1 since they have maximal endomorphism ring
                    self.add_nf_curve(self.field.F(r), t=0, f_E=1)
            
            # Iterate orders: get j-invariants, then add curves for matching isogeny classes
            nf_list = list(self.catalogue.NFC.data.values())
            _t0 = time.perf_counter()
            for nf in tqdm(nf_list, desc=f"HPC F_{self.field.q}", unit="nf", ncols=80, ascii=True):
                D_K = nf.discriminant
                # Get all orders for this number field at the current extension degree
                orders = nf.getOrders(self.field.n)
                for order in orders:
                    f = order.conductor
                    # Get j-invariants for this order
                    j_invs = get_j_invariants_from_order(D_K * f**2, f, self.field.q)
                    
                    print(f"{Colors.HEADER}\nOrder with D_K={D_K}, f={f} has j-invariants: {j_invs}{Colors.ENDC}")
                    for j_inv in j_invs:
                        order.add_j_invariant(j_inv)
                        # For each isogeny class with matching D_K, add curves if conductor divides f_pi
                        for t in order.traces:
                            print(f"{Colors.HEADER}t ={t}, Adding curves for order with D_K={D_K}, f={f}, j={j_inv}{Colors.ENDC}")
                            self.add_nf_curve(j_inv, t=t, f_E=f)
        else:
            
            _t0 = time.perf_counter()
            precompute_conductor = True
            
            #valid = [0, 1728] + [50, 667, 39]
            for j in tqdm(self.field.F, total=self.field.q, desc=f"F_{self.field.q}", unit="j", ncols=80, ascii=True):
                #if j == 0 or (j - 1728).is_zero() or j in valid:
                self.add_curves_by_j(j, pre_compute_conductor=precompute_conductor)
                
        #for ic in self.catalogue.isogeny_classes():
            #if ic.empty:
                #print(f"{Colors.WARNING}Warning: isogeny class with t={ic.t} is empty after enumeration{Colors.ENDC}")
        global highest_ell
        NE = 2*(self.field.q -2) + gcd(4, self.field.q-1) + gcd(6, self.field.q-1)  # Hasse bound for number of points, so we have at most NE curves per isogeny class, but usually much less
        if self.catalogue.size != NE:
            print(f"{Colors.FAIL}Warning: total number of curves in catalogue ({self.catalogue.size}) does not match expected number from Hasse bound ({NE}), there may be duplicates or missing curves{Colors.ENDC}")
        else:
            print(f"{Colors.GREEN}Successfully enumerated curves, total size of catalogue: {self.catalogue.size}, highest_ell={highest_ell}{Colors.ENDC}")
        
        print(f"Enumeration done in {time.perf_counter() - _t0:.2f}s")
                
    def add_curves_by_j(self, j, t: Optional[int] = None, f_E: Optional[int] = None, compute_twists: bool = True, pre_compute_conductor=False) -> None:
        # WE CREATE A TRUE GEOMETRIC CURVE WITH WEIERSTRASS FORM, INCLUDING A FULL ELLIPTICCURVE INSTANCE FROM SAGE
        aut_grp = self._get_aut_group_for_j(j)
        E = GeometricCurve(self.field, j, aut_grp=aut_grp, t=t, f_E=f_E)
        ell_t = self.catalogue.get_isogeny_class(E.t)
        # we compute the conductor for the first curve (they all share this, hence will get reused)
        if pre_compute_conductor:
            from utils.common import Config
            E.compute_conductor(ell_t.f_pi, use_true_height=Config.use_true_height)
        twists = E.compute_twists() if compute_twists else [E]  # note, compute_twists also returns a list with E included
        #print(f"f_E={E.f_E}, Adding curve with j={j}, is_j1728={E.is_j1728}, t={E.t}, f_E={E.f_E}, is_ss={E.is_supersingular}, twists={len(twists)}")
        for E_t in twists:
            self.catalogue.add(E_t)
            
    def add_nf_curve(self, j, t: int, f_E: int) -> None:
        # A MUCH SIMPLIFIED CURVE ONLY CONSISITING OF J AND T
        aut_grp = self._get_aut_group_for_j(j)
        self.catalogue.add(NFCurve(self.field, j, aut_grp=aut_grp, t=t, f_E=f_E))
    
    def compute_torsion(self, max_ell: int = 50, compute_volcano: bool = False) -> None:
        from tqdm import tqdm
        import time
        primes = list(primerange(2, max_ell))
        isogeny_classes = self.catalogue.isogeny_classes()
        _t0 = time.perf_counter()
        for ell_t in tqdm(isogeny_classes, desc="computing torsion", unit="ic", ncols=80, ascii=True):
            for ell in primes:
                if ell == self.field.p:
                    continue
                ell_t.eval_torsion_at_ell(ell, compute_volcano=compute_volcano)
        print(f"Computed torsion in {time.perf_counter() - _t0:.2f}s")
    
    
    
            
    def toJSON(self) -> Dict[str, Any]:
        return {
            "char": int(self.field.p),
            "catalogue": self.catalogue.toJSON()
        }