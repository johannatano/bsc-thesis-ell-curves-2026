
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