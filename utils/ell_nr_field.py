
from sage.all import *
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Set, Any
import math
from sympy import primerange
from utils.common import Colors


def toID(id) -> str:
    return ''.join(str(id).split())
    
checked_js = set()  # global set to track which j-invariants have been processed
@dataclass
class FqData:
    """Finite field data container for F_q = F_{p^n}"""
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
    """Return coefficients of x ∈ GF(p^n) in the base basis (length = n)."""
    poly = x.polynomial()
    coeffs = tuple(int(c) for c in poly.list())
    return coeffs + (0,)*(x.parent().degree() - len(coeffs))

def poly_to_tuple(poly) -> Tuple:
    """Return coefficients of x ∈ GF(p^n) in the base basis (length = n)."""
    coeffs = tuple(int(c) for c in poly.list())
    return coeffs + (0,)*(poly.degree() - len(coeffs))

def flatten(x, base: int) -> int:
    """Return a flattened index for x in GF(p^n) by interpreting its coefficients as digits in base p."""
    coeffs = element_to_tuple(x)
    index = 0
    for c in reversed(coeffs):
        index = index * base + c
    return index

#################### CURVE CLASS W HELPER METHODS ####################
def get_DK(D: int) -> int:
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
        key = (int(ell), int(q))
        if key not in cls._bivariate_cache:
            R = PolynomialRing(GF(int(q)), ['X', 'Y'])
            phi_ZZ = classical_modular_polynomial(int(ell))
            cls._bivariate_cache[key] = R(phi_ZZ)
        return cls._bivariate_cache[key]

    @classmethod
    def eval(cls, ell: int, j):
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
    @staticmethod
    def construct(k: int):
        R = PolynomialRing(ZZ, ['X', 'Y'])
        X, Y = R.gens()
        terms = [X**(k-i) * Y**i for i in range(k+1)]
        return R(sum(t for t in terms))
    @staticmethod
    def dickson_recursive(k: int, t: int, q: int) -> int:
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
    @staticmethod
    def _build(g, size_G: int, i: int) -> Dict:
        """Generic helper: builds orbit, stabilizer, and coset info for exponent i."""
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


# =============================================================================
# Number Field and Endomorphism Order structures
# =============================================================================



# =============================================================================
# Isogeny volcano structures
# =============================================================================

@dataclass   
class Isogeny():
    _type: str
    domain: str
    codomain: str
    def edge(self) -> Tuple[str, str]:
        return (self.domain, self.codomain)
    
@dataclass
class IsogenyVolcanoLevel:
    h: int
    vertrices: list = field(default_factory=list)  # list of j-invariants at this level
    edges: list = field(default_factory=list)      # list of edges (j1, j2) at this level
    def toJSON(self) -> Dict[str, Any]:
        return {
           "h": int(self.h),
           "vertrices": [v for v in self.vertrices],
           "edges": [ e for e in self.edges ]
        }
    
class IsogenyVolcano:
    def __init__(self, ell: int, height: int, fx_pi = None, N:int = 0, fx_pi_2 = None, N_2:int = 0) -> None:
        self.ell: int = ell
        self.levels: List[IsogenyVolcanoLevel] = []
        for h in range(height+1):
            self.addLevel(h)  # height above floor, so floor is h=0
        self.fx_pi = fx_pi
        self.fx_pi_2 = fx_pi_2
        self.N = N
        self.N_2 = N_2
        self.roots = {}
        
    def hasStructure(self) -> bool:
        return self.N % self.ell == 0 or self.N_2 % self.ell == 0
        
    def addLevel(self, h: int = 0) -> None:
        self.levels.append(IsogenyVolcanoLevel(h=h))

    def addVertrices(self, level: int, curves: List) -> None:
        for c in curves:
            self.levels[level].vertrices.append(c.ID)
    def _get_level_by_ID(self, ID: str) -> Optional[int]:
        for i, level in enumerate(self.levels):
            if ID in level.vertrices:
                return i
        return None 
    

                    
                    
    def addIsogeny(self, E1: str, E2: str) -> None:
        level = self._get_level_by_ID(E1)
        if level is not None:
            iso = Isogeny(_type="horizontal", domain=E1, codomain=E2)
            self.levels[level].edges.append(iso.edge())
        
    def toJSON(self) -> Dict[str, Any]:
        roots = self.fx_pi.roots(multiplicities=True) if self.fx_pi is not None else None
        roots_2 = self.fx_pi_2.roots(multiplicities=True) if self.fx_pi_2 is not None else None
        
        # if self.hasStructure() else None
        return {
            "ell": int(self.ell),
            "-fx_roots": [int(r) for r, mult in roots for _ in range(mult)] if roots is not None else [],
            "+fx_roots": [int(r) for r, mult in roots_2 for _ in range(mult)] if roots_2 is not None else [],
            "-": self.N % self.ell == 0,
            "+": self.N_2 % self.ell == 0,
            #"roots": self.roots,
            "levels": [level.toJSON() for level in self.levels]
        }

@dataclass
class EndomorphismOrder:
    """t-independent: stores orders and j-invariants from CM theory"""
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
    
    '''def add_ic(self, ic:str) -> None:
        """Add a isogeny class ID to this order"""
        self.ic.add(ic)'''
    
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
    
    This is t-dependent data: stores all curves with trace t over F_q,
    organized by conductor of their endomorphism ring.
    """
    
    def __init__(self, t: int, q: int) -> None:
        self.t: int = t
        self.q: int = q
        self.N_pts = q + 1 - t
        self.N_pts_2 = q + 1 + t
        self.n: int = int(ZZ(q).factor()[0][1]) if q > 1 else 1  # extract n from q = p^n
        self.p: int = int(ZZ(q).factor()[0][0])
        self.D_pi = ZZ(t*t - 4*q)
        self.ordinary = (self.t % self.p) != 0
        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        self.fx_pi = x**2 - t*x + q
        self.fx_pi_2 = x**2 + t*x + q
        self.volcanoes = {}
        self.O_K = None
        self.f_pi = 1
        self.K = None
        self.D_K = 0
        
        self.orders: Dict[str, EndomorphismOrder] = {}  # dict of conductor -> EndomorphismOrder
        
        # For ordinary curves, ensure all orders exist for this number field
        # and add trace t to each order where conductor divides f_pi
        '''if ell_t.ordinary and ell_t.f_pi is not None:
            tree = field.getTreeByN(n)
            for d in sorted(ell_t.f_pi.divisors(), reverse=False):
                order = tree.ensure_order_exists(D_K, field.K, d)
                order.add_ic(ell_t.ID)'''
                
        # Compute number field info but don't create order objects here
        if self.ordinary:
            self.K = NumberField(self.fx_pi, 'x')          # same field as Q(sqrt(Dpi))
            self.D_K = ZZ(self.K.discriminant())        # fundamental discriminant
            f_pi2 = ZZ(self.D_pi // self.D_K)
            #if(self.D_K == -3):
            #    print(f"Creating order for D_K={self.D_K}, f_pi={self.f_pi}, t={t}")
            self.f_pi = f_pi2.isqrt()
            
            for d in sorted(self.f_pi.divisors(), reverse=False):
                order = self.ensure_order_exists(self.D_K, self.K, d)
                #order.add_ic(ell_t.ID)
                
            self.O_K = self.K.maximal_order()
        
        self.generic = True
        self.empty = True
        self.curves_by_order = {}  # dict: conductor (f_E) -> list of Curve objects
        self.curves = {}  # dict: curve_ID -> Curve (for quick lookup)
        
        self.ID = f"t{self.t}_n{self.n}"
        
    def ensure_order_exists(self, D_K: int, K, conductor: int) -> EndomorphismOrder:
        """Ensure an order exists for the given conductor, create if needed"""
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
        return self.f_pi.valuation(l)
    
    def full_rank_primes(self) -> int:
        primes = []
        for _l, e in self.f_pi.factor():
            if self.N_pts % _l**2 == 0:
                primes.append(_l)
        return primes
    
    def D_K_Legedre(self, l: int) -> int:
        return kronecker(self.D_K, l)
        
    def H(self) -> int:
        return pari(4*self.q - self.t*self.t).qfbhclassno()
    
    def eval_hk_mod_fx(self, k, fx, t):
        # h_k(Frob, Frob_dual) mod fx
        # 1. Grab the ring (R) and the variable (x) from the existing fx_pi
        R = fx.parent()
        x = R.gen()
        # 2. Now define your values using THIS x
        frob = x
        frob_dual = t - x
        # Substitute
        result_multi = hk_symbolic.subs(X=frob, Y=frob_dual)
        # Force it into the correct Ring (R) before dividing
        result_uni = R(result_multi)
        # Now the division works
        final_value = result_uni.quo_rem(fx)[1]
        #print("H_k(Frob, Frob_dual) mod P:", final_value)
        return final_value
    
    
    def eval_torsion_at_ell(self, ell: int, compute_volcano:bool = False) -> Tuple:
        fx_l = self.fx_pi.change_ring(Zmod(ell))
        fx_l_2 = self.fx_pi_2.change_ring(Zmod(ell))
        total_height = self.l_adic_height(ell) if self.ordinary else 0
        
        volcano = IsogenyVolcano(ell, total_height, fx_l, self.q + 1 - self.t, fx_l_2, self.q + 1 + self.t)
        self.volcanoes[ell] = volcano
        
        #N_pts = self.q + 1 - self.t
        #N_pts_2 = self.q + 1 + self.t
        
        if not volcano.hasStructure():
            return fx_l, 0
        #print(f"{Colors.HEADER}Evaluating torsion at ell={ell} for t={self.t}, q={self.q}, N_pts={N_pts}, ordinary={self.ordinary}{Colors.ENDC}")
        #for r, m in fx_l.roots(multiplicities=True):
            #dual_is_q = (r-self.q) % ell == 0
            #print(f"{Colors.HEADER if (r == 1 or dual_is_q) else Colors.FAIL}dual_is_q={dual_is_q}, r={r}, ell={ell}{Colors.ENDC}")
        #D_K_Legendre={self.D_K_Legedre(ell) tells me if its a square mod l?
        #print(f"{Colors.GREEN if total_height > 0 else Colors.WARNING}ell={ell} t={self.t}, q={self.q}, N_pts={N_pts}, SS={not self.ordinary}, f_pi={fx_l.roots(multiplicities=True)}, fx_l={fx_l.factor()} l-adic height={total_height} D_K={self.D_K}, D_K_Legendre={self.D_K_Legedre(ell)}{Colors.ENDC}")
                
        # and total_height > 0
        if self.ordinary:
            '''#for c in self.getCurves():
                #if ell > 50:
                    #print(f"{Colors.GREEN if c.is_j0 or c.is_j1728 else Colors.BLUE}ell={ell} Curve {c.j} has {len(c.E.points())} points, is_ss={c.E.is_supersingular()}, inv={c.E.abelian_group().invariants()}{Colors.ENDC}")
                if(c.j == 0 or c.is_j1728):
                    pts = [P for P in c.E.points() if P.order() == ell]
                    npts = len(pts)
                    rk2 = 2 if self.f_pi % ell == 0 else 1
                    rk_true = 2 if npts > ell else 1
                    #x = PolynomialRing(c.F, 'x').gen()
                    #poly = x**3 + c.E.a6()
                    #print(poly, poly.roots())
                    print(f"{Colors.FAIL if rk2 != rk_true else Colors.GREEN}Curve {c.j} has rank {rk2 if npts > ell else (1 if npts > 0 else 1)}, npts={npts} points of order {ell}, is_ss={c.E.is_supersingular()}, inv={c.E.abelian_group().invariants()}{Colors.ENDC}")'''
            #early exit, not needed for correctness but computationally efficient for batch jobs, ie no volcano height, surface = floor at ell
            #if (self.q - 1) % ell != 0 or total_height == 0: 
            #    self.volcanoes[ell].addVertrices(0, self.getCurves())
            #    return fx_l, 0
            for f, curves_list in self.curves_by_order.items():
                if f > 0:
                    h = total_height - ZZ(f).valuation(ell)
                    #print(f"t={self.t}, f_pi={self.f_pi}, f_E={f}, ell={ell}, valuation={ZZ(f).valuation(ell)}, height={h}, j={[c.j for c in curves_list]}")
                    #rk = 2 if (self.f_pi / f) % ell == 0 else 1 
                    self.volcanoes[ell].addVertrices(h, curves_list)
                    
                    if compute_volcano:
                        for c in curves_list:
                            self.compute_volcano(ell, c)
                    #print(f"{Colors.GREEN if rk == 2 else Colors.BLUE}Conductor f={f} gives height h={h} and predicted rank {rk} at ell={ell}{Colors.ENDC}")
                    #for c in curves_list:
                    #    verification_pts = [P for P in c.E.points() if P.order() == ell]
                    #    n_pts = len(verification_pts)
                    #    rank = 2 if n_pts > ell else (1 if n_pts > 0 else 0)  # each point has a negative, so divide by 2
                    #    clr = Colors.GREEN if c.is_j0 or c.is_j1728 else Colors.BLUE
                    #    if rank != rk:
                    #        print(f"{clr}h={h}, Curve {c.j} has {n_pts} points of order {ell}, giving rank {rank}{Colors.ENDC}")
                else:
                    for c in curves_list:
                        h = c.height_above_floor(ell, total_height, use_true_height=False)
                        self.volcanoes[ell].addVertrices(h, [c])
        else:
            print(f"{Colors.WARNING}Supersingular case: no volcano structure, all curves at floor level{Colors.ENDC}")
            self.volcanoes[ell].addVertrices(0, self.getCurves())
        return fx_l, total_height
    
    def compute_volcano(self, ell, curve) -> None:
        # TODO: Implement isogeny edge computation
        mod_poly = MOD_POLY.eval(ell, curve.j)
        roots = mod_poly.roots(multiplicities=True)
        #self.roots[curve.ID] = [(str(r), m) for r, m in roots]
        for r, m in roots:
            #print(f"At ell={self.ell}, curve {curve.j} has root {r} with multiplicity {m} in the modular polynomial")    
            for k in range(m):
                target_curve = self.getCurveByJ(r)
                if target_curve is not None:
                    self.volcanoes[ell].addIsogeny(curve.ID, target_curve.ID)
                
    def add_curve(self, curve) -> None:
        # Compute conductor if not already set
        level = curve.f_E if curve.f_E is not None else -1
        
        # Initialize curves list for this conductor if needed
        if level not in self.curves_by_order:
            self.curves_by_order[level] = []
        
        # Generate curve ID based on conductor and position in list
        curr_len = len(self.curves_by_order[level])
        curve.ID = f"f{level}_{curr_len}"
        
        # Add curve to both data structures
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
        result = {
            
            #"ID": self.ID,
            "t": int(self.t),
            #"D_pi": int(self.D_pi),
            "f_pi": int(self.f_pi),
            #"ordinary": self.ordinary,
            #"generic": self.generic
            "O": [order.toJSON() for order in self.orders.values()]
        }
        
        '''if self.f_pi is not None:
            result["D_pi"] = int(self.D_pi)
            result["f_pi"] = int(self.f_pi)'''
        # 
        if include_curves:
            # Export curves organized by conductor
            result["curves"] = {
                int(f): [c.toJSON() for c in curves]
                for f, curves in self.curves_by_order.items()
            }
            result["volcanoes"] = [vol.toJSON() for ell, vol in self.volcanoes.items() if vol.hasStructure()]
        
        return result

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int, n: int) -> 'IsogenyClass':
        """Lightweight deserializer that avoids expensive CM/Sage order construction."""
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
        ic.ID = f"t{ic.t}_n{ic.n}"
        return ic
    


class NumberFieldTree:
    """Stores t-independent data: endomorphism orders and j-invariants for degree n.
    
    Organizes orders by conductor, each containing j-invariants from Hilbert class polynomial
    and traces that correspond to curves with that endomorphism ring.
    """
    
    def __init__(self, n: int) -> None:
        self.n: int = n
        self.isogeny_classes = []  # t-dependent: dict of (t, n) -> IsogenyClass
        #self.orders: Dict[str, EndomorphismOrder] = {}  # dict of conductor -> EndomorphismOrder
    
    '''def ensure_order_exists(self, D_K: int, K, conductor: int) -> EndomorphismOrder:
        """Ensure an order exists for the given conductor, create if needed"""
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
        return self.orders[conductor_key]'''

        
    def toJSON(self, include_curves: bool = True) -> Dict[str, Any]:
        return {
            "n": int(self.n),
            "I_t": [ic.toJSON(include_curves=include_curves) for ic in self.isogeny_classes]
            #"orders": [order.toJSON() for order in self.orders.values()]
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int) -> 'NumberFieldTree':
        n = int(data.get("n", 1))
        tree = cls(n=n)
        for ic_data in data.get("I_t", []) or []:
            tree.isogeny_classes.append(IsogenyClass.fromJSON(ic_data, p=p, n=n))
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
        """Get isogeny class by trace and extension degree"""
        tree = self.getTreeByN(n)
        for ic in tree.isogeny_classes:
            if ic.t == t:
                return ic
    
    def addIsogenyClass(self, ic: IsogenyClass) -> None:
        """Add isogeny class to registry"""
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


class NumberFieldCatalogue:
    """Global catalogue of all imaginary quadratic fields and their associated data.
    
    Central registry that manages:
    - Number field creation and lookup by discriminant
    - Isogeny class creation and trace assignment to orders
    - Cross-referencing between orders and isogeny classes
    """
    
    def __init__(self, p: int) -> None:
        self.p: int = p
        self.data: Dict[int, NumberFieldData] = {}

    def create_isogeny_class(self, t: int, n: int) -> IsogenyClass:
        """Create isogeny class and populate corresponding endomorphism orders.
        
        Note: This method has a circular dependency with IsogenyClass from elliptic.py.
        IsogenyClass is imported within the method to avoid circular imports.
        """
        
        
        # Check if isogeny class for t already exists
        ell_t = self.get_isogeny_class(t, n)
        if ell_t:
            return ell_t
        
        # Import here to avoid circular dependency
        
        # Create new isogeny class
        ell_t = IsogenyClass(t=t, q=self.p**n)
        D_K = ell_t.D_K
        # Get or create number field
        field = self.getFieldByDiscriminant(D_K)
        #tree = self.getTreeByN(n)
        # Add isogeny class to the number field's registry
        field.addIsogenyClass(ell_t)
        return ell_t
    
    def getFieldByDiscriminant(self, D: int) -> NumberFieldData:
        field = self.data.get(D)
        return field if field is not None else self.addField(D)
    
    def addField(self, D: int) -> NumberFieldData:
        self.data[D] = NumberFieldData(dk=D)
        return self.data[D]
    
    def get_isogeny_class(self, t: int, n: int) -> Optional[IsogenyClass]:
        """Get isogeny class by searching all number fields"""
        for nf_info in self.data.values():
            ell_t = nf_info.getIsogenyClass(t, n)
            if ell_t is not None:
                return ell_t
        return None
    
    def get_isogeny_classes_by_n(self, n: int) -> List[IsogenyClass]:
        """Return all isogeny classes for a given n"""
        classes = []
        for nf_info in self.data.values():
            tree = nf_info.getTreeByN(n)
            classes.extend(tree.isogeny_classes)
        return classes
    
    def getCurvesByJ(self, j, n: Optional[int] = None) -> List:
        """Get all curves with given j-invariant across all isogeny classes"""
        if n is None:
            n = self.N
        curves = []
        for nf_info in self.data.values():
            tree = nf_info.getTreeByN(n)
            for ic in tree.isogeny_classes:
                curve = ic.getCurveByJ(j)
                if curve:
                    curves.append(curve)
        return curves

    def sort(self) -> None:
        """Sort number fields by discriminant (smallest absolute value first)"""
        self.data = dict(sorted(self.data.items(), key=lambda item: abs(item[0])))

    def toJSON(self) -> Dict[str, List]:
        return {
            "nf": [{
                "D": int(nf.discriminant),
                "tree": [tree.toJSON() for tree in nf.tree],
                #"isogeny_classes": [ic.toJSON(include_curves=False) for ic in nf.isogeny_classes.values()]
            } for dk, nf in self.data.items() ]
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int) -> 'NumberFieldCatalogue':
        catalogue = cls(p=int(p))
        catalogue.data = {}
        for nf_data in data.get("nf", []) or []:
            nf = NumberFieldData.fromJSON(nf_data, p=int(p))
            catalogue.data[int(nf.discriminant)] = nf
        catalogue.sort()
        return catalogue


class NumberFieldsClassifier_Fq:
    """Pre-computes isogeny class structure for multiple finite fields.
    
    Generates all possible traces within Hasse bound and creates corresponding
    isogeny classes with their endomorphism orders.
    """
    def __init__(self, char: int) -> None:
        self.char: int = char
        self.nr_fields: NumberFieldCatalogue = NumberFieldCatalogue(self.char)

    def generate(self, p_powers: List[int], q_max: int = 1000000, t_list: Optional[List[int]] = None) -> NumberFieldCatalogue:
        from utils.common import Colors
        
        for n in p_powers:
            self.nr_fields.N = n
            used_ts = set()
            q = self.char ** n
            print(f"{Colors.HEADER}Generating isogeny classes for F_{q}{Colors.ENDC}")
            
            if t_list is not None:
                print(f"{Colors.HEADER}Using custom trace list with {len(t_list)} traces{Colors.ENDC}")
                for t in t_list:
                    if abs(t) > 2*math.isqrt(q):
                        print(f"{Colors.WARNING}Warning: skipping trace t={t} as it exceeds the Hasse bound for F_{q}{Colors.ENDC}")
                        continue
                    if t % self.char == 0:
                        print(f"{Colors.WARNING}Warning: skipping trace t={t} as it is divisible by the characteristic{Colors.ENDC}")
                        continue
                    self.nr_fields.create_isogeny_class(t, n)
                    #self.nr_fields.create_isogeny_class(-t, n)
                continue
            
            if(q > q_max):
                print(f"{Colors.WARNING}Warning: skipping initialization for F_{q} as it exceeds the current limit{Colors.ENDC}")
                continue
            
            HB = math.isqrt(4*q)
            max_prime = q + 1 + HB
            primes = list(primerange(2, min(max_prime + 1, q_max)))
            
            from tqdm import tqdm
            for ell in tqdm(primes, desc=f"F_{q} primes", unit="ell", leave=False, ncols=80, ascii=True):
                i_min = (q + 1 - HB + ell - 1) // ell
                i_max = (q + 1 + HB) // ell
                for i in tqdm(range(i_min, i_max + 1), desc=f"ell={ell}", unit="i", leave=False, ncols=80, ascii=True):
                    t = q + 1 - i*ell
                    if( t % self.char == 0):
                        continue
                    if abs(t) in used_ts:
                        continue
                    used_ts.add(abs(t))
                    self.nr_fields.create_isogeny_class(abs(t), n)
                    #self.nr_fields.create_isogeny_class(-t, n)
        
        print(f"{Colors.GREEN}Finished generating isogeny classes for all specified p-powers{Colors.ENDC}")
        self.nr_fields.sort()
        return self.nr_fields
    
    def toJSON(self) -> Dict[str, Any]:
        return {
            "char": int(self.char),
            "nr_fields": self.nr_fields.toJSON()
        }
        
    @classmethod
    def fromJson(cls, data: Dict[str, Any]) -> 'NumberFieldsClassifier_Fq':
        """Build a classifier instance from serialized JSON data.

        Expected shape:
        {
            "char": <int>,
            "nr_fields": { "nf": [...] }
        }
        """
        if not isinstance(data, dict):
            raise TypeError("fromJson expects a dict payload")

        char = int(data.get("char", 0))
        if char <= 0:
            raise ValueError("Invalid or missing 'char' in payload")

        obj = cls(char=char)
        nr_fields_data = data.get("nr_fields", {}) or {}
        obj.nr_fields = NumberFieldCatalogue.fromJSON(nr_fields_data, p=char)
        return obj