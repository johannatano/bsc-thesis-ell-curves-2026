"""Helpers for ℓ-torsion subgroups of elliptic curves over finite fields.

This module focuses on the local computations needed after a curve has already
been constructed: torsion rank, orbit counting under automorphisms, and simple
serialization of the resulting subgroup data.
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

class EPOrbit:
    """Minimal container for an elliptic-point orbit size."""
    def __init__(self, aut_size: int) -> None:
        self.aut_size: int = aut_size
       
class Point:
    """Wrapper around a Sage point with cached automorphism-orbit data."""
    def __init__(self, P) -> None:
        self.P = P
        self._orbit: Optional[List] = None
        self._stab_size: int = 1

    def phi(self, u):
        """Apply the standard automorphism action $(x,y) \mapsto (u^2x,u^3y)$."""
        if self.P.is_zero(): return self.P
        x, y = self.P.xy()
        return (u^2 * x, u^3 * y)

    def orbit(self, G) -> List:
        """Return the orbit of the point under the supplied root-of-unity group."""
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
    """Compute basic ℓ-torsion data for a single elliptic curve.

    The class supports two main tasks used elsewhere in the project:
    determining the rank of `E[ℓ](F_q)` and counting automorphism orbits of
    nonzero torsion points.
    """
    
    def __init__(self, curve, l: int) -> None:
        self.curve = curve
        self.l: int = l
        self.gens: Optional[List] = None
        self.orbits: List[Point] = []
        self.rank: int = 0
        self.n_orbits: int = 0
        self.points: List[Point] = []
        self.special_aut6_orbit: bool = False

    def compute_rank(self, f_pi: int, use_generators: bool = False) -> None:
        """Compute the ℓ-torsion rank, either from generators or valuation data."""
        if self.curve.N_pts % self.l != 0:
            self.rank = 0
            return
        if use_generators:
            self.gens = self._get_generators()
            self.rank = len(self.gens)
        else:
            self.rank = self._compute_rank(f_pi)

    def generate_orbits(self):
        """Enumerate nonzero torsion points from the chosen generators."""
        if self.rank == 2:
            for i in range(self.l):
                for j in range(self.l):
                    self._add(i * self.gens[0] + j * self.gens[1])
        else:
            for i in range(self.l):
                self._add(i * self.gens[0])
                
    def count_orbits(self) -> int:
        """Count automorphism orbits of nonzero ℓ-torsion points."""
        return self._count_orbits_level_2() if self.l == 2 else (self._count_orbits_level_3() if self.l == 3 else self._count_orbits_general())
    
    def _count_orbits_level_2(self) -> int:
        """Closed-form orbit count for ℓ = 2."""
        if(self.curve.aut_size == 2):
            return 2**self.rank - 1
        elif(self.curve.aut_size == 4):
            return 2**(self.rank -1)
        else:
            return 1
        
    def _count_orbits_level_3(self) -> int:
        """Closed-form orbit count for ℓ = 3, including the extra `j=0` symmetry."""
        fixed = 3**self.rank - 1
        if(self.curve.aut_size == 6):
            B = self.curve.getCoefficients()[1]
            fixed += 2*(1+B.is_square())
        return (fixed // self.curve.aut_size)
    
    def _count_orbits_general(self) -> int:
        """Generic orbit count when no exceptional low-level case is needed."""
        return (self.l**self.rank-1) // self.curve.aut_size

    def _add(self, ell_P) -> None:
        """Add a torsion point if it represents a new automorphism orbit."""
        if ell_P.is_zero():
            return
        P = Point(ell_P)
        self.points.append(P)
        if self._check_unique_orbit(P):
            self.orbits.append(P)
            self.n_orbits += 1
    
    def _check_unique_orbit(self, torsion_point: Point) -> bool:
        for S in self.orbits:
            for P in S.orbit(self.curve.aut_grp):
                if torsion_point.point().xy() == P:
                    return False
        return True
    
    def _rank_by_group_structure(self) -> int:
        """Read the ℓ-rank from the abelian group invariants of `E(F_q)`."""
        invariants = self.curve.getSageCurve().abelian_group().invariants()
        r = 0
        for inv in invariants:
            if inv % self.l == 0:
                r += 1
        return r
    
    def _rank_by_modular_poly(self) -> int:
        """Estimate the rank from linear factors of the modular polynomial specialization."""
        mod_poly = classical_modular_polynomial(self.l, self.curve.j)
        factors = mod_poly.factor()
        nr_linear_factors = 0
        for (f, m) in factors:
            if f.degree() == 1:
                nr_linear_factors += m
        return 2 if nr_linear_factors > 1 else 1
    
    def _rank_by_enum_points(self) -> int:
        """Fallback rank estimate by explicit enumeration of rational points."""
        torsion_pts = [P for P in self.curve.getSageCurve().points() if not P.is_zero() and (self.l * P).is_zero()]
        return 1 if len(torsion_pts) == (self.l-1) else 2  # rough estimate, may overcount if not full rank
    
    def _two_torsion_rank(self) -> int:
        """Specialized 2-torsion rank test via splitting of the Weierstrass cubic."""
        if (self.curve.N_pts % 2 != 0):
            return 0
        w_poly = self.curve.weierstrass_polynomial()
        splits = all(f.degree() == 1 for f, m in w_poly.factor())
        return 2 if splits else 1
        
    def _compute_rank(self, f_pi) -> int:
        """Compute the expected torsion rank from the isogeny-class conductor data."""
        if self.curve.is_supersingular:
            return self._rank_by_group_structure()
        return 2 if ZZ(self.curve.f_E).valuation(self.l) < ZZ(f_pi).valuation(self.l) else 1

    def _get_generators(self) -> List:
        """Return generators of the rational ℓ-torsion subgroup from Sage."""
        # Adapted from Sage's torsion basis logic, but allowing rank 1 output.
        return [P.element() for P in self.curve.getSageCurve().abelian_group().torsion_subgroup(self.l).gens()]
        
    def toJSON(self) -> Dict[str, int]:
        return {"rank": self.rank}
