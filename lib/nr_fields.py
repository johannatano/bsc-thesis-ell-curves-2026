
"""Arithmetic data structures for isogeny classes over finite fields.

This module contains the number-field side of the project: finite-field helper
data, endomorphism orders, trace-indexed isogeny classes, and ℓ-isogeny
volcano information. These structures are used by the curve classifiers but can
also be serialized independently.
"""

from sage.all import *
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Set, Any
import math
from sympy import primerange
from utils.common import Colors

def max_ell_from_HB(q: int) -> int:
    HB = math.isqrt(4*q)
    return q + 1 + HB
    
def toID(id) -> str:
    return ''.join(str(id).split())
    
checked_js = set()  # global set to track which j-invariants have been processed
@dataclass
class FqData:
    """Basic container for the finite field F_{p^n}."""
    p: int
    n: int
    q: int = field(init=False)
    F: GF = field(init=False)
    
    def __post_init__(self) -> None:
        self.q = self.p ** self.n
        self.F = GF(self.q)
        self.g = self.F.multiplicative_generator()
        
    def multiplicative_generator(self):
        return self.F.multiplicative_generator()
    
    def extend(self, n: int) -> 'FqData':
        return FqData(self.p, self.n * n)
    
#################### GENERIC HELPER ####################
def element_to_tuple(x) -> Tuple:
    """Return coordinates of an element of GF(p^n) in the standard basis."""
    poly = x.polynomial()
    coeffs = tuple(int(c) for c in poly.list())
    return coeffs + (0,)*(x.parent().degree() - len(coeffs))

def poly_to_tuple(poly) -> Tuple:
    """Return the coefficient tuple of a polynomial."""
    coeffs = tuple(int(c) for c in poly.list())
    return coeffs + (0,)*(poly.degree() - len(coeffs))

def flatten(x, base: int) -> int:
    """Encode a field element as a single integer via its base-`p` coordinates."""
    coeffs = element_to_tuple(x)
    index = 0
    for c in reversed(coeffs):
        index = index * base + c
    return index

#################### CURVE CLASS W HELPER METHODS ####################
def get_DK(D: int) -> int:
    """Extract the fundamental discriminant part from an order discriminant."""
    squarefree_part = D
    #print(D)
    for p, e in factor(abs(D)):
        if(p == 2):
            e = max(e-2, 0)
        #print(f"Prime factor {p} with exponent {e}")
        if e >= 2:
            squarefree_part //= p**(e - (e % 2))  # divide by even powers
    D0 = squarefree_part
    #print(f"Square-free part of D is {D0}")
    return D0
from math import gcd


class MOD_POLY:
    """Session-local cache for modular polynomials.

    Caches:
      - bivariate Phi_ell(X,Y) reduced to GF(q), keyed by (ell, q)
      - univariate Phi_ell(x, j) in GF(q)[x], keyed by (ell, q, j_id)
    """

    _bivariate_cache: Dict[Tuple[int, int], Any] = {}
    _eval_cache: Dict[Tuple[int, int, str], Any] = {}

    @classmethod
    def construct(cls, ell: int, q: int):
        """Return `Phi_ell(X,Y)` reduced to `GF(q)` and cached by `(ell,q)`."""
        key = (int(ell), int(q))
        if key not in cls._bivariate_cache:
            R = PolynomialRing(GF(int(q)), ['X', 'Y'])
            phi_ZZ = classical_modular_polynomial(int(ell))
            cls._bivariate_cache[key] = R(phi_ZZ)
        return cls._bivariate_cache[key]

    @classmethod
    def eval(cls, ell: int, j):
        """Specialize the modular polynomial at the given j-invariant."""
        F = j.parent()
        q = int(F.cardinality())
        j_id = toID(j)
        key = (int(ell), q, j_id)
        if key not in cls._eval_cache:
            x = polygen(F)
            phi = cls.construct(ell, q)
            cls._eval_cache[key] = phi([x, j])
        return cls._eval_cache[key]

    @classmethod
    def clear(cls) -> None:
        cls._bivariate_cache.clear()
        cls._eval_cache.clear()

class Hk:
    """Helpers for the symmetric polynomials appearing in Hecke computations."""
    @staticmethod
    def construct(k: int):
        """Build the polynomial h_k(X,Y)=sum_{i=0}^k X^{k-i}Y^i."""
        R = PolynomialRing(ZZ, ['X', 'Y'])
        X, Y = R.gens()
        terms = [X**(k-i) * Y**i for i in range(k+1)]
        return R(sum(t for t in terms))
    @staticmethod
    def dickson_recursive(k: int, t: int, q: int) -> int:
        """Evaluate the Dickson-style recurrence for `h_k` without building the full polynomial."""
        # Base cases for h_k
        if k == 0: return 1
        if k == 1: return t
        # Recurrence: h_k = t * h_{k-1} - n * h_{k-2}
        # This avoids building the huge polynomial object
        hk_prev, hk_curr = 1, t
        for _ in range(k - 1):
            hk_prev, hk_curr = hk_curr, t * hk_curr - q * hk_prev
        return hk_curr
      
class RootsOfUnity:
    """Factories for the root-of-unity groups governing curve automorphisms."""
    @staticmethod
    def _build(g, size_G: int, i: int) -> Dict:
        """Build the subgroup of `i`-th roots of unity inside `F_q^*`."""
        grp = [(g ** (size_G // i)) ** k for k in range(i)]
        return dict(
            size=i,
            grp=grp,
        )

    @staticmethod
    def toJSON(dict: Dict) -> Optional[Dict]:
        return {
            "size": dict['size'],
            "grp": [element_to_tuple(u) for u in dict['grp']],
        } if dict is not None else None

    # --- Static entry points per j-case ---

    @staticmethod
    def make_2(g, size_G: int) -> Dict:
        """Generic case: exponent i = 2 (A' = u^4A, B' = u^6B)"""
        return RootsOfUnity._build(g, size_G, i=2)

    @staticmethod
    def make_4(g, size_G: int) -> Dict:
        """j = 1728 case: exponent i = 4  (A' = u^4A, B = 0)"""
        return RootsOfUnity._build(g, size_G, i=4)

    @staticmethod
    def make_6(g, size_G: int) -> Dict:
        """j = 0 case: exponent i = 6  (A = 0, B' = u^6B)"""
        return RootsOfUnity._build(g, size_G, i=6)

@dataclass   
class Isogeny():
    _type: str
    domain: str
    codomain: str
    def edge(self) -> Tuple[str, str]:
        return (self.domain, self.codomain)
    
@dataclass
class IsogenyVolcanoLevel:
    """Single level of an ℓ-isogeny volcano graph."""
    h: int
    vertrices: list = field(default_factory=list)  # list of j-invariants at this level
    edges: list = field(default_factory=list)      # list of edges (j1, j2) at this level
    def toJSON(self) -> Dict[str, Any]:
        return {
           "h": int(self.h),
           "vertrices": [v for v in self.vertrices],
           "edges": [ e for e in self.edges ]
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any]) -> 'IsogenyVolcanoLevel':
        return cls(
            h=int(data.get("h", 0)),
            vertrices=list(data.get("vertrices", data.get("vertices", [])) or []),
            edges=[tuple(edge) for edge in (data.get("edges", []) or [])],
        )
    
class IsogenyVolcano:
    """Compressed representation of the ℓ-isogeny volcano for one trace class."""
    def __init__(self, ell: int, height: int, fx_pi = None, N:int = 0) -> None:
        self.ell: int = ell
        self.levels: List[IsogenyVolcanoLevel] = []
        for h in range(height+1):
            self.addLevel(h)  # height above floor, so floor is h=0
        self.fx_pi = fx_pi
        #self.fx_pi_2 = fx_pi_2
        self.N = N
        #self.N_2 = N_2
        self.roots = {}
        
    def hasStructure(self) -> bool:
        """Return whether the class has nontrivial ℓ-adic volcano structure."""
        return self.N % self.ell == 0# or self.N_2 % self.ell == 0
        
    def addLevel(self, h: int = 0) -> None:
        self.levels.append(IsogenyVolcanoLevel(h=h))

    def addVertrices(self, level: int, curves: List) -> None:
        """Attach curves to a given volcano level."""
        for c in curves:
            self.levels[level].vertrices.append(c.ID)
            
    def _get_level_by_ID(self, ID: str) -> Optional[int]:
        for i, level in enumerate(self.levels):
            if ID in level.vertrices:
                return i
        return None 
                    
    def addIsogeny(self, E1: str, E2: str) -> None:
        """Add an edge between two serialized curve identifiers."""
        level = self._get_level_by_ID(E1)
        if level is not None:
            iso = Isogeny(_type="horizontal", domain=E1, codomain=E2)
            self.levels[level].edges.append(iso.edge())
        
    def toJSON(self) -> Dict[str, Any]:
        """Serialize a single-sign volcano in the older JSON format."""
        roots = self.fx_pi.roots(multiplicities=True) if self.fx_pi is not None else None
        #roots_2 = self.fx_pi_2.roots(multiplicities=True) if self.fx_pi_2 is not None else None
        # if self.hasStructure() else None
        return {
            "ell": int(self.ell),
            "fx_roots": [int(r) for r, mult in roots for _ in range(mult)] if roots is not None else [],
            #"-fx_roots": [int(r) for r, mult in roots for _ in range(mult)] if roots is not None else [],
            #"+fx_roots": [int(r) for r, mult in roots_2 for _ in range(mult)] if roots_2 is not None else [],
            #"-": self.N % self.ell == 0,
            #"+": self.N_2 % self.ell == 0,
            #"roots": self.roots,
            "levels": [level.toJSON() for level in self.levels]
        }

    @staticmethod
    def _roots_to_list(fx_pi) -> List[int]:
        roots = fx_pi.roots(multiplicities=True) if fx_pi is not None else None
        return [int(r) for r, mult in roots for _ in range(mult)] if roots is not None else []

    @staticmethod
    def toJSON_pair(positive_trace_volcano: Optional['IsogenyVolcano'], negative_trace_volcano: Optional['IsogenyVolcano']) -> Dict[str, Any]:
        """Serialize the pair of volcanoes for `t` and `-t` into one compressed entry."""
        template = positive_trace_volcano if positive_trace_volcano is not None else negative_trace_volcano
        return {
            "ell": int(template.ell),
            "-fx_roots": IsogenyVolcano._roots_to_list(positive_trace_volcano.fx_pi) if positive_trace_volcano is not None else [],
            "+fx_roots": IsogenyVolcano._roots_to_list(negative_trace_volcano.fx_pi) if negative_trace_volcano is not None else [],
            "-": positive_trace_volcano.hasStructure() if positive_trace_volcano is not None else False,
            "+": negative_trace_volcano.hasStructure() if negative_trace_volcano is not None else False,
            "levels": [level.toJSON() for level in template.levels],
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], sign: str = '-') -> 'IsogenyVolcano':
        """Restore a single signed volcano view from serialized JSON."""
        ell = int(data.get("ell", 0))
        levels_payload = data.get("levels", []) or []
        height = max((int(level.get("h", 0)) for level in levels_payload), default=0)

        roots_key = f"{sign}fx_roots"
        roots_payload = data.get(roots_key)
        if roots_payload is None:
            roots_payload = data.get("fx_roots", []) or []

        R = PolynomialRing(Zmod(ell), 'x')
        x = R.gen()
        fx_pi = R.one()
        for root in roots_payload:
            fx_pi *= (x - Zmod(ell)(int(root)))

        has_structure = bool(data.get(sign, data.get("has_structure", True)))
        volcano = cls(ell=ell, height=height, fx_pi=fx_pi, N=(ell if has_structure else 1))
        volcano.levels = [IsogenyVolcanoLevel.fromJSON(level_data) for level_data in levels_payload]
        return volcano


@dataclass
class SerializedCurve:
    """Lightweight curve placeholder used when reading JSON back into memory."""
    payload: Dict[str, Any]
    conductor: int = -1
    ID: str = field(init=False)
    j: Any = field(init=False)
    is_j0: bool = field(default=False, init=False)
    is_j1728: bool = field(default=False, init=False)
    aut_size: int = field(default=2, init=False)

    def __post_init__(self) -> None:
        self.ID = str(self.payload.get("ID", ""))
        self.j = self.payload.get("j")
        self.f_E = int(self.conductor)

    def toJSON(self) -> Dict[str, Any]:
        return dict(self.payload)

@dataclass
class EndomorphismOrder:
    """CM order data shared across trace classes with the same quadratic field."""
    conductor: int
    class_number: int
    D_K: int = 0  # fundamental discriminant of the number field
    j_invariants: list = field(default_factory=list)  # list of j-invariants that belong to this order
    #ic: set = field(default_factory=set)  # set of isogeny class ids t that belong to this order
    
    def add_j_invariant(self, j) -> None:
        """Add a j-invariant to this order if not already present"""
        # Compare by converting to tuple for hashability
        j_tuple = element_to_tuple(j) if hasattr(j, 'minpoly') else j
        for existing_j in self.j_invariants:
            existing_j_tuple = element_to_tuple(existing_j) if hasattr(existing_j, 'minpoly') else existing_j
            if j_tuple == existing_j_tuple:
                return
        self.j_invariants.append(j)

    
    def toJSON(self) -> Dict[str, Any]:
        return {
            "f": int(self.conductor),
            "cn": int(self.class_number),
            "D": int(self.conductor**2 * self.D_K)
            #"j": [element_to_tuple(j) for j in self.j_invariants],
            #"isogeny_classes": list(self.ic)
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any]) -> 'EndomorphismOrder':
        f = int(data.get("f", 1))
        cn = int(data.get("cn", 0))
        D_order = int(data.get("D", 0))
        D_K = int(D_order // (f * f)) if f != 0 else 0
        return cls(conductor=f, class_number=cn, D_K=D_K)


@dataclass
class IsogenyClass:
    """Isogeny class of elliptic curves with a given trace of Frobenius.

    This is the main trace-dependent object. It stores the Frobenius data, the
    compatible endomorphism orders, the curves in the class, and any precomputed
    volcano information.
    """
    def __init__(self, t: int, q: int) -> None:
        self.t: int = t
        self.q: int = q
        self.N_pts = q + 1 - t
        #self.N_pts_2 = q + 1 + t
        self.n: int = int(ZZ(q).factor()[0][1]) if q > 1 else 1  # extract n from q = p^n
        self.p: int = int(ZZ(q).factor()[0][0])
        self.D_pi = ZZ(t*t - 4*q)
        self.ordinary = (self.t % self.p) != 0
        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        self.fx_pi = x**2 - t*x + q
        #self.fx_pi_2 = x**2 + t*x + q
        self.volcanoes = {}
        self.O_K = None
        self.f_pi = 1
        self.K = None
        self.D_K = 0
        self.orders: Dict[str, EndomorphismOrder] = {}  # dict of conductor -> EndomorphismOrder
        # Compute number field info but don't create order objects here
        if self.ordinary:
            self.K = NumberField(self.fx_pi, 'x')          # same field as Q(sqrt(Dpi))
            self.D_K = ZZ(self.K.discriminant())        # fundamental discriminant
            f_pi2 = ZZ(self.D_pi // self.D_K)
            self.f_pi = f_pi2.isqrt()
            for d in sorted(self.f_pi.divisors(), reverse=False):
                order = self.ensure_order_exists(self.D_K, self.K, d)
            self.O_K = self.K.maximal_order()
        self.generic = True
        self.empty = True
        self.curves_by_order = {}  # dict: conductor (f_E) -> list of Curve objects
        self.curves = {}  # dict: curve_ID -> Curve (for quick lookup)
        
        self.ID = f"t{self.t}_n{self.n}"
        
    def ensure_order_exists(self, D_K: int, K, conductor: int) -> EndomorphismOrder:
        """Ensure that the order of the given conductor has been constructed."""
        conductor_key = str(conductor)
        if conductor_key not in self.orders:
            # Create the order from the number field
            order_sage = K.order_of_conductor(conductor)
            h_O = order_sage.class_number()
            self.orders[conductor_key] = EndomorphismOrder(
                conductor=int(conductor),
                class_number=int(h_O),
                D_K=int(D_K)
            )
        return self.orders[conductor_key]
    
    def l_adic_height(self, l: int) -> int:
        """Return the ℓ-adic valuation of the Frobenius conductor `f_pi`."""
        return self.f_pi.valuation(l)
    
    def full_rank_primes(self) -> int:
        """Return primes for which the class can have full rational ℓ-torsion."""
        primes = []
        for _l, e in self.f_pi.factor():
            if self.N_pts % _l**2 == 0:
                primes.append(_l)
        return primes
    
    def D_K_Legedre(self, l: int) -> int:
        return kronecker(self.D_K, l)
        
    def H(self) -> int:
        """Return the Hurwitz class number attached to the Frobenius discriminant."""
        return pari(4*self.q - self.t*self.t).qfbhclassno()
    
    def eval_hk_mod_fx(self, ell: int, hk_symbolic):
        """Reduce the Hecke polynomial expression modulo the Frobenius polynomial."""
        fx = self.fx_pi
        R = fx.parent()
        x = R.gen()
        frob = x
        frob_dual = self.t - x
        result_multi = hk_symbolic.subs(X=frob, Y=frob_dual)
        result_uni = R(result_multi)
        final_value = result_uni.quo_rem(fx)[1]
        return final_value    
    
    def compute_volcano(self, ell: int, edges:bool = False) -> IsogenyVolcano:
        """Build the ℓ-isogeny volcano for this isogeny class."""
        fx_l = self.fx_pi.change_ring(Zmod(ell))
        total_height = self.l_adic_height(ell) if self.ordinary else 0
        volcano = IsogenyVolcano(ell, total_height, fx_l, self.q + 1 - self.t)
        self.volcanoes[ell] = volcano
        if not volcano.hasStructure():
            return volcano
        if self.ordinary:
            for f, curves_list in self.curves_by_order.items():
                if f > 0:
                    h = total_height - ZZ(f).valuation(ell)
                    self.volcanoes[ell].addVertrices(h, curves_list)
                    if edges:
                        for c in curves_list:
                            self.compute_volcano_edges(ell, c)
                else:
                    for c in curves_list:
                        h = c.height_above_floor(ell, total_height, use_true_height=False)
                        self.volcanoes[ell].addVertrices(h, [c])
        else:
            self.volcanoes[ell].addVertrices(0, self.getCurves())
        return volcano
    
    def compute_volcano_edges(self, ell, curve) -> None:
        """Add outgoing isogeny edges for one curve using the modular polynomial."""
        mod_poly = MOD_POLY.eval(ell, curve.j)
        roots = mod_poly.roots(multiplicities=True)
        for r, m in roots:
            for k in range(m):
                target_curve = self.getCurveByJ(r)
                if target_curve is not None:
                    self.volcanoes[ell].addIsogeny(curve.ID, target_curve.ID)
                
    def add_curve(self, curve) -> None:
        """Insert a curve into the class, grouped by endomorphism conductor."""
        level = curve.f_E if curve.f_E is not None else -1
        
        if level not in self.curves_by_order:
            self.curves_by_order[level] = []
        curr_len = len(self.curves_by_order[level])
        curve.ID = f"f{level}_{curr_len}"
        self.curves_by_order[level].append(curve)
        self.curves[curve.ID] = curve
        self.empty = False
        self.generic = self.generic and (not curve.is_j0) and (not curve.is_j1728)
        
    def getCurves(self, conductor: Optional[int] = None) -> List:
        if conductor is not None:
            return self.curves_by_order.get(conductor, [])
        else:
            return list(self.curves.values())
        
    def getCurveByJ(self, j) -> Optional[Any]:
        for curve in self.curves.values():
            if curve.j == j:
                return curve
        return None
         
    def toJSON(self, include_curves: bool = False) -> Dict[str, Any]:
        """Serialize a single signed isogeny class."""
        result = {
            "t": int(self.t),
            "f_pi": int(self.f_pi),
            "O": [order.toJSON() for order in self.orders.values()]
        }
        if include_curves:
            result["curves"] = {
                int(f): [c.toJSON() for c in curves]
                for f, curves in self.curves_by_order.items()
            }
            result["volcanoes"] = [vol.toJSON() for ell, vol in self.volcanoes.items() if vol.hasStructure()]
        return result

    @staticmethod
    def toCompressedJSON(positive_trace_ic: Optional['IsogenyClass'], negative_trace_ic: Optional['IsogenyClass'], include_curves: bool = False) -> Dict[str, Any]:
        """Serialize the pair of classes for `t` and `-t` into one JSON entry."""
        template = positive_trace_ic if positive_trace_ic is not None else negative_trace_ic
        curve_source = template
        if curve_source is not None and curve_source.empty and negative_trace_ic is not None and not negative_trace_ic.empty:
            curve_source = negative_trace_ic

        result = {
            "t": int(abs(template.t)),
            "f_pi": int(template.f_pi),
            "O": [order.toJSON() for order in template.orders.values()]
        }
        if include_curves:
            result["curves"] = {
                int(f): [c.toJSON() for c in curves]
                for f, curves in curve_source.curves_by_order.items()
            } if curve_source is not None else {}

            positive_ells = set(positive_trace_ic.volcanoes.keys()) if positive_trace_ic is not None else set()
            negative_ells = set(negative_trace_ic.volcanoes.keys()) if negative_trace_ic is not None else set()
            result["volcanoes"] = [
                IsogenyVolcano.toJSON_pair(
                    positive_trace_ic.volcanoes.get(ell) if positive_trace_ic is not None else None,
                    negative_trace_ic.volcanoes.get(ell) if negative_trace_ic is not None else None,
                )
                for ell in sorted(positive_ells | negative_ells)
                if ((positive_trace_ic is not None and ell in positive_ells and positive_trace_ic.volcanoes[ell].hasStructure()) or
                    (negative_trace_ic is not None and ell in negative_ells and negative_trace_ic.volcanoes[ell].hasStructure()))
            ]
        return result

    def _restore_curves_from_json(self, curves_payload: Dict[str, Any]) -> None:
        """Restore lightweight curve payloads from serialized JSON."""
        self.curves_by_order = {}
        self.curves = {}
        for conductor_key, curves_list in (curves_payload or {}).items():
            conductor = int(conductor_key)
            restored_curves = [SerializedCurve(dict(curve_data), conductor=conductor) for curve_data in (curves_list or [])]
            self.curves_by_order[conductor] = restored_curves
            for curve in restored_curves:
                self.curves[curve.ID] = curve
        self.empty = len(self.curves) == 0

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int, n: int) -> 'IsogenyClass':
        """Lightweight deserializer for a single signed isogeny class."""
        t = int(data.get("t", 0))
        q = int(p) ** int(n)

        ic = cls.__new__(cls)
        ic.t = t
        ic.q = q
        ic.N_pts = q + 1 - t
        ic.n = int(n)
        ic.p = int(p)
        ic.D_pi = ZZ(t * t - 4 * q)
        ic.ordinary = (ic.t % ic.p) != 0

        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        ic.fx_pi = x**2 - t * x + q

        ic.volcanoes = {}
        ic.O_K = None
        ic.K = None
        ic.D_K = 0
        ic.f_pi = 1
        ic.orders = {}

        orders_payload = data.get("O", []) or []
        f_candidates = []
        for order_data in orders_payload:
            order = EndomorphismOrder.fromJSON(order_data)
            ic.orders[str(order.conductor)] = order
            f_candidates.append(int(order.conductor))
            if ic.D_K == 0 and int(order.conductor) != 0:
                ic.D_K = int(order.D_K)

        if f_candidates:
            ic.f_pi = int(max(f_candidates))

        if ic.D_K == 0 and ic.f_pi != 0:
            ic.D_K = int(ic.D_pi // (ic.f_pi * ic.f_pi))

        ic.generic = True
        ic.empty = True
        ic.curves_by_order = {}
        ic.curves = {}
        ic._restore_curves_from_json(data.get("curves", {}) or {})
        active_sign = '-' if ic.t >= 0 else '+'
        for volcano_data in data.get("volcanoes", []) or []:
            volcano = IsogenyVolcano.fromJSON(volcano_data, sign=active_sign)
            ic.volcanoes[volcano.ell] = volcano
        ic.ID = f"t{ic.t}_n{ic.n}"
        return ic

    @classmethod
    def expandFromJSON(cls, data: Dict[str, Any], p: int, n: int) -> List['IsogenyClass']:
        """Expand one compressed `|t|` JSON entry into signed in-memory classes."""
        t = int(data.get("t", 0))
        volcanoes_payload = data.get("volcanoes", []) or []
        has_signed_volcanoes = any(("-fx_roots" in volcano_data) or ("+fx_roots" in volcano_data) for volcano_data in volcanoes_payload)

        if t == 0 or not has_signed_volcanoes:
            return [cls.fromJSON(data, p=p, n=n)]

        positive_payload = dict(data)
        positive_payload["t"] = abs(t)

        negative_payload = dict(data)
        negative_payload["t"] = -abs(t)

        return [
            cls.fromJSON(positive_payload, p=p, n=n),
            cls.fromJSON(negative_payload, p=p, n=n),
        ]
    


class NumberFieldTree:
    """Trace-indexed isogeny-class data for one extension degree `n`."""
    
    def __init__(self, n: int) -> None:
        self.n: int = n
        self.isogeny_classes = []  # t-dependent: dict of (t, n) -> IsogenyClass
        #self.orders: Dict[str, EndomorphismOrder] = {}  # dict of conductor -> EndomorphismOrder

    def _group_by_abs_trace(self) -> Dict[int, Dict[str, Optional[IsogenyClass]]]:
        """Group signed isogeny classes by absolute trace for compact serialization."""
        grouped: Dict[int, Dict[str, Optional[IsogenyClass]]] = {}
        for ic in self.isogeny_classes:
            abs_t = abs(int(ic.t))
            if abs_t not in grouped:
                grouped[abs_t] = {"positive": None, "negative": None}
            sign_key = "negative" if int(ic.t) < 0 else "positive"
            grouped[abs_t][sign_key] = ic
        return grouped
        
    def toJSON(self, include_curves: bool = True) -> Dict[str, Any]:
        """Serialize the tree, compressing `t` and `-t` into one entry per `|t|`."""
        grouped = self._group_by_abs_trace()
        return {
            "n": int(self.n),
            "I_t": [
                IsogenyClass.toCompressedJSON(group["positive"], group["negative"], include_curves=include_curves)
                for _, group in sorted(grouped.items(), key=lambda item: item[0])
            ]
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int) -> 'NumberFieldTree':
        n = int(data.get("n", 1))
        tree = cls(n=n)
        for ic_data in data.get("I_t", []) or []:
            tree.isogeny_classes.extend(IsogenyClass.expandFromJSON(ic_data, p=p, n=n))
        return tree


class NumberFieldData:
    """Container for a single imaginary quadratic field Q(√D_K).
    
    Separates:
    - t-independent data: NumberFieldTree (orders and j-invariants by extension degree)
    - t-dependent data: isogeny classes indexed by (trace, n)
    """
    
    def __init__(self, dk: Optional[int] = None) -> None:
        from utils.common import Colors
        self.discriminant: Optional[int] = dk
        d = ZZ(self.discriminant)
        if d != 0:
            self.K = QuadraticField(d, 'a')
            self.O_K = self.K.maximal_order()
        else:
            #print(f"{Colors.WARNING}Warning: trivial number field with discriminant 0 created{Colors.ENDC}")
            self.K = None
            self.O_K = None
        self.tree = []  # NumberFieldTree objects (t-independent: orders and j-invariants by n)
        
    
    def getOrders(self, n: int) -> List[EndomorphismOrder]:
        """Get unique endomorphism orders for a given n.

        Also annotates each returned order with `traces` for HPC usage.
        """
        tree = self.getTreeByN(n)
        merged: Dict[int, EndomorphismOrder] = {}

        for ic in tree.isogeny_classes:
            for order in ic.orders.values():
                key = int(order.conductor)
                if key not in merged:
                    merged[key] = EndomorphismOrder(
                        conductor=int(order.conductor),
                        class_number=int(order.class_number),
                        D_K=int(order.D_K),
                    )
                    merged[key].traces = set()
                merged[key].traces.add(int(ic.t))

        out = [merged[k] for k in sorted(merged.keys())]
        for order in out:
            if hasattr(order, "traces") and isinstance(order.traces, set):
                order.traces = sorted(order.traces)
        return out
    
    def getIsogenyClass(self, t: int, n:int) -> Optional[IsogenyClass]:
        """Return the isogeny class for the given trace and extension degree."""
        tree = self.getTreeByN(n)
        for ic in tree.isogeny_classes:
            if ic.t == t:
                return ic
    
    def addIsogenyClass(self, ic: IsogenyClass) -> None:
        """Attach an isogeny class to the appropriate degree-`n` tree."""
        tree = self.getTreeByN(ic.n)
        tree.isogeny_classes.append(ic)
    
    def getTreeByN(self, n: int) -> NumberFieldTree:
        for tree in self.tree:
            if tree.n == n:
                return tree
        return self.createNTree(n)
    
    def createNTree(self, n: int) -> NumberFieldTree:
        tree = NumberFieldTree(n)
        self.tree.append(tree)
        return tree

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int) -> 'NumberFieldData':
        nf = cls(dk=int(data.get("D", 0)))
        nf.tree = []
        for tree_data in data.get("tree", []) or []:
            nf.tree.append(NumberFieldTree.fromJSON(tree_data, p=p))
        return nf


