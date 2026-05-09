from __future__ import annotations
from functools import cached_property, lru_cache
from fractions import Fraction
import math
from sage.libs.pari import pari


# TODO: make small db for this
@lru_cache(maxsize=4096)
def _class_number(D: int) -> int:
    return int(pari(D).qfbclassno())


class Order:
    def __init__(self, K: NumberField, f: int):
        self.K = K
        self.f = f

    @property
    def is_maximal(self) -> bool:
        return self.f == 1

    @cached_property
    def D(self) -> int:
        return self.f ** 2 * self.K.D

    @cached_property
    def unit_index(self) -> int:
        """[O_K^× : O_f^×] — equals 1 for f>1 except field has extra roots of unity."""
        if self.is_maximal:
            return 1
        return self.K._maximal_unit_index

    @cached_property
    def h(self) -> int:
        raise NotImplementedError

    def order_of(self, f: int) -> Order:
        """Return the sub-order of conductor f (must divide self.f)."""
        return self.K.order(f)

class QuadraticOrder(Order):
    @cached_property
    def h(self) -> int:
        if self.is_maximal:
            return _class_number(self.K.D)
        result = Fraction(self.f * self.K.maximal_order.h, self.unit_index)
        for p, _ in _factorize(self.f):
            result *= Fraction(p - legendre(self.K.D, p), p)
        assert (
            result.denominator == 1
        ), f"Non-integer class number: D_K={self.K.D}, f={self.f}, result={result}"
        return int(result)

    @cached_property
    def hw(self) -> int:
        """Returns weight by 1/|Aut| (Schoof 1987, Thm 4.6)."""
        return Fraction(
            self.h, self.K._maximal_unit_index * 2 if self.is_maximal else 2
        )


class QuaternionOrder(Order):
    @cached_property
    def h(self) -> int:
        """Curve count in Q_{∞,p} (Schoof 1987, Thm 4.6).
        """
        if self.f != 1:
            raise NotImplementedError("Non-maximal quaternion orders not implemented")
        return Fraction(
            self.K.p + 6 - 4 * legendre(-3, self.K.p) - 3 * legendre(-4, self.K.p), 12
        )

    @cached_property
    def hw(self) -> int:
        """ Curve count in Q_{∞,p} weighted with 1/|Aut| (Schoof 1987, Thm 4.6)."""
        if self.f != 1:
            raise NotImplementedError("Non-maximal quaternion orders not implemented")
        return Fraction(self.K.p - 1, 12)


class NumberField:
    @staticmethod
    def __D0(D: int) -> int:
        if D == 0:
            return 0
        sign = -1 if D < 0 else 1
        n = abs(D)
        sf = 1
        d = 2
        while d * d <= n:
            exp = 0
            while n % d == 0:
                n //= d
                exp += 1
            if exp % 2 == 1:
                sf *= d
            d += 1
        if n > 1:
            sf *= n
        sf *= sign
        return sf if sf % 4 == 1 else 4 * sf

    def __init__(self, D: int) -> None:
        self._disc = NumberField.__D0(D)
        self._order_cache: dict[int, Order] = {}

    @cached_property
    def D(self) -> int:
        return self._disc

    @cached_property
    def O_K(self) -> Order:
        return self.order(1)

    def h(self) -> int:
        return self.O_K.h

    def order(self, f: int) -> Order:
        if f not in self._order_cache:
            self._order_cache[f] = Order(self, f)
        return self._order_cache[f]

# --- Imaginary quadratic fields and quaternion algebras ---

class ImaginaryQuadraticField(NumberField):
    @staticmethod
    @lru_cache(maxsize=4096)
    def H(D: int) -> int:
        """Kronecker class number H(Δ) — unweighted count of SL2(Z)-orbits (Schoof 1987)."""
        if D >= 0:
            return 0
        total = 0
        abs_D = abs(D)
        for d in range(1, math.isqrt(abs_D) + 1):
            if abs_D % (d * d) == 0:
                delta = D // (d * d)
                if delta % 4 == 0 or delta % 4 == 1:
                    total += _class_number(delta)
        return total

    def __init__(self, D: int) -> None:
        super().__init__(D)
        self._maximal_unit_index = (
            3 if self._disc == -3 else (2 if self._disc == -4 else 1)
        )


class QuaternionAlgebra(NumberField):

    @staticmethod
    @lru_cache(maxsize=4096)
    def H(self, total_mass: bool = False) -> Fraction:
        """Curve count in Q_{∞,p} (Schoof 1987, Thm 4.6).
        total_mass=False: weighted count with 1/|Aut| for j=0 and j=1728.
        total_mass=True:  full mass, simplifies to (p-1)/12.
        """
        if total_mass:
            return Fraction(p - 1, 12)
        chi3 = _legendre(-3, p)
        chi4 = _legendre(-4, p)
        return Fraction(p + 6 - 4 * chi3 - 3 * chi4, 12)

    def __init__(self, p:int) -> None:
        super().__init__(0)
        self.p = p
