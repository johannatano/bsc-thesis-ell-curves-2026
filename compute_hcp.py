#!/usr/bin/env python
"""
Standalone script to compute curves for a field using HCP with Classifier_Fq.
Usage: sage -python compute_hcp.py <p> <n> <D>
"""
import sys
import json
import os
from pathlib import Path
from sage.all import *
from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial

# Import your custom classes
sys.path.insert(0, str(Path(__file__).parent))

# Suppress stdout during import and Classifier_Fq creation to avoid polluting JSON output
from contextlib import redirect_stdout
with redirect_stdout(sys.stderr):
    from utils.elliptic import Classifier_Fq

def get_j_invariants_from_order(D_0, f, q):
    """Compute j-invariants from HCP for order with discriminant D_0 and conductor f."""
    D = D_0 * f**2
    j_invs = []
    try:
        H = hilbert_class_polynomial(D)
        F_q = GF(q)
        for j in H.change_ring(F_q).roots(multiplicities=False):
            j_invs.append(j)  # Keep as field element, not string
    except Exception as e:
        print(f"Warning: Could not compute HCP for D={D}: {e}", file=sys.stderr)
    return j_invs

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: sage -python compute_hcp.py <p> <n> <D>")
        sys.exit(1)
    
    p = int(sys.argv[1])
    n = int(sys.argv[2])
    D_0 = int(sys.argv[3])
    q = p**n
    
    # Load nr_fields data
    DATA_ROOT = Path(__file__).parent / "data"
    nr_fields_file = DATA_ROOT / str(p) / "nr_fields.json"
    
    if not nr_fields_file.exists():
        print(json.dumps({'success': False, 'error': f'No data for p={p}'}))
        sys.exit(1)
    
    with open(nr_fields_file) as f:
        data = json.load(f)
        nr_fields_data = data['nr_fields']['nf']
    
    # Find the number field with this discriminant
    nf = None
    for field in nr_fields_data:
        if field['D'] == D_0:
            nf = field
            break
    
    if not nf:
        print(json.dumps({'success': False, 'error': f'D={D_0} not found'}))
        sys.exit(1)
    
    # Find tree for this n
    tree = None
    for t in nf.get('tree', []):
        if t['n'] == n:
            tree = t
            break
    
    if not tree:
        print(json.dumps({'success': False, 'error': f'n={n} not found'}))
        sys.exit(1)
    
    try:
        # Create Classifier_Fq using YOUR custom class
        # Redirect stdout to stderr to avoid polluting JSON output with debug prints
        print(f"Creating Classifier_Fq for p={p}, n={n}", file=sys.stderr)
        with redirect_stdout(sys.stderr):
            CFq = Classifier_Fq(p, n)
        
        # Enumerate curves using HCP
        added_js = set()
        for ic in tree.get('ic', []):
            trace = ic['t']
            if 'orders' in ic:
                for order in ic['orders']:
                    f = order['f']
                    # Get j-invariants for this order
                    j_invs = get_j_invariants_from_order(D_0, f, q)
                    
                    print(f"Found {len(j_invs)} j-invariants for order with D_0={D_0}, f={f}", file=sys.stderr)
                    for j_inv in j_invs:
                        if j_inv not in added_js:
                            # Use YOUR add_curve method
                            CFq.add_curve(j_inv, trace, f)
                            added_js.add(j_inv)
        CFq.compute_volcanoes(min_height=1, max_ell=50)  # Adjust parameters as needed
        # Output using YOUR toJSON() format
        catalogue_json = CFq.toJSON()
        
        result = {
            'success': True,
            'catalogue': catalogue_json
        }
        
        print("Serializing JSON...", file=sys.stderr)
        output = json.dumps(result, default=str)
        print(f"JSON size: {len(output)} bytes", file=sys.stderr)
        print(output)
        
    except Exception as e:
        import traceback
        print(json.dumps({
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc()
        }), file=sys.stderr)
        sys.exit(1)
