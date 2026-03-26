"""Construction of number-field and isogeny-class data over finite fields.

This module prepares the arithmetic catalogue indexed by Frobenius traces. It
does not enumerate concrete curve models; instead, it builds the number-field,
order, and isogeny-class scaffolding later used by the curve classifiers.
"""

from sage.all import *
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Set, Any
import math
from sympy import primerange
from utils.common import Colors
from lib.nr_fields import *
import time
from tqdm import tqdm

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
        """Create the signed isogeny class for trace `t` if it does not yet exist."""

        ell_t = self.get_isogeny_class(t, n)
        if ell_t:
            return ell_t

        ell_t = IsogenyClass(t=t, p=self.p, n=n)
        D_K = ell_t.D_K
        field = self.getFieldByDiscriminant(D_K)
        field.addIsogenyClass(ell_t)
        return ell_t

    def getFieldByDiscriminant(self, D: int) -> NumberFieldData:
        field = self.data.get(D)
        return field if field is not None else self.addField(D)

    def addField(self, D: int) -> NumberFieldData:
        self.data[D] = NumberFieldData(dk=D)
        return self.data[D]

    def get_isogeny_class(self, t: int, n: int) -> Optional[IsogenyClass]:
        """Search the whole catalogue for the isogeny class with trace `t`."""
        for nf_info in self.data.values():
            ell_t = nf_info.getIsogenyClass(t, n)
            if ell_t is not None:
                return ell_t
        return None

    def get_isogeny_classes_by_n(self, n: int) -> List[IsogenyClass]:
        """Return all signed isogeny classes for the given extension degree."""
        classes = []
        for nf_info in self.data.values():
            tree = nf_info.getTreeByN(n)
            classes.extend(tree.isogeny_classes)
        return classes

    def getCurvesByJ(self, j, n: Optional[int] = None) -> List:
        """Return all curves with the given j-invariant across the selected degree."""
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
        """Serialize the full catalogue grouped by quadratic discriminant."""
        return {
            "nf": [
                {
                    "D": int(nf.discriminant),
                    "tree": [tree.toJSON() for tree in nf.tree],
                    # "isogeny_classes": [ic.toJSON(include_curves=False) for ic in nf.isogeny_classes.values()]
                }
                for dk, nf in self.data.items()
            ]
        }

    @classmethod
    def fromJSON(cls, data: Dict[str, Any], p: int) -> "NumberFieldCatalogue":
        catalogue = cls(p=int(p))
        catalogue.data = {}
        for nf_data in data.get("nf", []) or []:
            nf = NumberFieldData.fromJSON(nf_data, p=int(p))
            catalogue.data[int(nf.discriminant)] = nf
        catalogue.sort()
        return catalogue



class NumberFieldsClassifier_Fq:
    """Build the arithmetic catalogue of isogeny classes over fields of characteristic `p`.
    
    For each extension degree `n`, the classifier determines relevant traces and
    creates the associated number fields, endomorphism orders, and signed
    isogeny classes.
    """
    def __init__(self, char: int) -> None:
        self.char: int = char
        self.nr_fields: NumberFieldCatalogue = NumberFieldCatalogue(self.char)

    def generate(self, p_powers: List[int], q_max: int = 1000000, t_list: Optional[List[int]] = None) -> NumberFieldCatalogue:
        """Populate the catalogue for the requested extension degrees.

        If `t_list` is provided, only those traces are used. Otherwise the code
        scans traces indirectly through primes in the Hasse interval and keeps
        one representative for each absolute trace.
        """
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
                    # Store both sign choices in memory; JSON serialization may
                    # later compress them back to a single |t| entry.
                    self.nr_fields.create_isogeny_class(abs(t), n)
                    self.nr_fields.create_isogeny_class(-abs(t), n)
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
                    if abs(t) in used_ts:
                        continue
                    used_ts.add(abs(t))
                    # Again, keep both signed trace classes available in memory.
                    self.nr_fields.create_isogeny_class(abs(t), n)
                    self.nr_fields.create_isogeny_class(-abs(t), n)

        print(f"{Colors.GREEN}Finished generating isogeny classes for all specified p-powers{Colors.ENDC}")
        self.nr_fields.sort()
        return self.nr_fields

    def toJSON(self) -> Dict[str, Any]:
        """Serialize the classifier and its number-field catalogue."""
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
