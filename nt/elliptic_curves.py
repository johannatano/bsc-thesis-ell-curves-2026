from functools import cached_property
from fractions import Fraction
import math

from nt.rings import ImaginaryQuadraticField, quaternion_class_number
from nt.common import legendre, factorize

class IsogenyClass:
    def __init__(self, q: int, t: int, N_t:int = None) -> None:
        self.q = q
        self.t = t
        self.__N_t = N_t

    @cached_property
    def D_pi(self) -> int:
        return self.t ** 2 - 4 * self.q

    @cached_property
    def N_pts(self) -> int:
        return self.q + 1 - self.t

    @cached_property
    def field(self) -> ImaginaryQuadraticField:
        return ImaginaryQuadraticField(self.D_pi)

    @cached_property
    def is_quaternion(self) -> bool:
        return self.D_pi == 0

    @cached_property
    def N_t(self) -> int:
        """Hurwitz-weighted count of curves in this isogeny class."""
        # TODO: add the SS cases from Schoof
        return (
            self.__N_t
            if self.__N_t is not None
            else ImaginaryQuadraticField.H(self.D_pi)
        )


class CurvesRecordFq:
    """All isogeny classes of elliptic curves over F_q."""

    def __init__(self, p: int, n:int) -> None:
        self.q = p**n
        self.p = p
        self.n = n
        self.j0_SS = (legendre(-3, p) == -1) # note: true for p = 2, not true for p = 3
        self.j1728_SS = (legendre(-4, p) == -1) # note: true for p = 2, not true for p = 3

        self.classes: dict[int, IsogenyClass] = {
            t: IsogenyClass(self.q, t, N_t)
            for t, N_t in self.__ordinary_ts + self.__supersingular_ts
        }

    @cached_property
    def HB(self) -> int:
        return math.isqrt(4 * self.q)

    @cached_property
    def __ordinary_ts(self) -> list[(int, int)]:

        """
        Ordinary traces t satisfy |t| < 2*sqrt(q) and p does not divide t.
        N(t) = H(t^2 - 4q) counts the number of curves in these isogeny classes, no inert factor scaling for ordinary curves.
        see: Schoof 1987, Thm 4.6 for reference

        Return a list of tuples (t, None) for ordinary traces, where N(t) will fallback to H(t^2 - 4q). NOTE: this is NOT weighted by 1/|Aut(E)|

        """

        return [(s * t, None) for t in range(1, self.HB + 1) if t % self.p != 0 for s in (1, -1)]

    @cached_property
    def __supersingular_ts(self) -> list[(int, int)]:
        """
        Supersingular traces t satisfy specific conditions depending on the characteristic and degree of the field.
        see: Schoof 1987, Thm 4.6 for reference on each case
        
        Return a list of tuples (t, N(t)) for supersingular traces, where N(t) is the number of curves in the isogeny class with trace t, including any inert factor scaling. NOTE: this is NOT weighted by 1/|Aut(E)|
        """

        ts: list[(int, int)] = []
        if self.n % 2 == 1:
            """
            For char(Fq) > 3, the only permitted SS trace is 0.
            This IsogenyClass will lie in field of D_K = -4p, and all curves will lie in max order, hence: N(t) = H(-4p)
            """
            ts.append((0, ImaginaryQuadraticField.H(-4 * self.p)))
            """
            For char(Fq) = 2, 3, we get t = pm sqrt(2q) or t = pm sqrt(3q) respectively, producing exactly one curve per such trace, N(t) = 1.
            """
            if self.p == 2 or self.p == 3:
                t = self.p ** ((self.n + 1) // 2)
                ts += [(s * t, 1) for s in (1, -1)]
        else:

            # N(t) = 1 - (-3/p), will yield D_K = -3 and only max order allowed, only curve here is j0 with inert factor 2
            if self.j0_SS:
                sqrt_q = self.p ** (self.n // 2)
                ts += [(s * sqrt_q, 2) for s in (1, -1)]

            # N(t) = 1 - (-4/p), will yield D_K = -4 and only max order allowed, only curve here is j1728 with inert factor 2
            if self.j1728_SS:
                ts.append((0, 2))

            # Quaternions, we compute the full class number as given in Schoof 1987, Thm 4.6
            num_curves = quaternion_class_number(self.p)
            ts += [(s * self.HB, num_curves) for s in (1, -1)]
        return ts

    @cached_property
    def num_curves_total(self) -> Fraction:
        return Fraction(0) + sum(c.N_t for c in self.classes.values())

    def check(self) -> bool:
        return self.num_curves_total == self.q
