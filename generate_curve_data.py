import argparse
import sys
import os
import time
import json
import sqlite3
from pathlib import Path

from utils.ell_curves import CurvesClassifier_Fq
from utils.ell_nr_field import NumberFieldsClassifier_Fq
from utils.common import Logger, Colors, Data, Config
from math import gcd
from sympy import primerange
from sage.all import *
import numpy as np
import requests
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from sage.all import DirichletGroup

def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument("-p", "--p", type=int, required=True, default=-1, help="Field char p")
    p.add_argument("-n", "--n", type=int, required=True, default=-1, help="Field extension degree n")
    p.add_argument("-t", "--t_list", type=int, nargs='+', required=True, help="List of traces t")
    p.add_argument("--to-db", action="store_true", default=False, help="Also persist output to SQLite cache (one row per discriminant D)")
    p.add_argument("--db-path", type=str, default="", help="Optional SQLite path (default: ./data/cache.sqlite3)")
    return p.parse_args()


def _cache_db_path(user_db_path: str = "") -> Path:
    if user_db_path:
        return Path(user_db_path)
    return Path(__file__).parent / "data" / "cache.sqlite3"


def _init_cache_db(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS curve_cache (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                p INTEGER NOT NULL,
                n INTEGER NOT NULL,
                d INTEGER NOT NULL,
                status TEXT NOT NULL DEFAULT 'ready',
                payload_json TEXT NOT NULL,
                created_at INTEGER NOT NULL,
                updated_at INTEGER NOT NULL,
                UNIQUE(p, n, d)
            )
            """
        )


def _save_catalogue_to_cache_db(p: int, n: int, payload: dict, db_path: Path) -> int:
    now = int(time.time())
    nfs = payload.get("catalogue", {}).get("number_fields", []) or []
    written = 0
    with sqlite3.connect(db_path) as conn:
        for nf in nfs:
            d = int(nf.get("D"))
            row_payload = {
                "char": int(p),
                "catalogue": {
                    "number_fields": [nf]
                }
            }
            payload_json = json.dumps(row_payload, separators=(",", ":"))
            conn.execute(
                """
                INSERT INTO curve_cache (p, n, d, status, payload_json, created_at, updated_at)
                VALUES (?, ?, ?, 'ready', ?, ?, ?)
                ON CONFLICT(p, n, d)
                DO UPDATE SET
                    status = excluded.status,
                    payload_json = excluded.payload_json,
                    updated_at = excluded.updated_at
                """,
                (int(p), int(n), int(d), payload_json, now, now),
            )
            written += 1
    return written

def generate(p: int, n: int, t_list: list[int], to_db: bool = False, db_path: str = ""):
    q_max = 10**20  # Set a maximum q to avoid long computations
    nf = None
    # we prepopulate a catalogue of possible traces, hence discriminants ie nr fields for this char, up to a max n
    print(f"\n=== Processing prime p={p}, t={t_list}, loading json.... ===")
    q = p**n
    NFC = NumberFieldsClassifier_Fq(p)
    nf = NFC.generate([n], q_max=q+1, t_list=t_list)
        # if we do not load precomputed json, first generate the number fields catalogue for this p and n
        # create the classifier for curves over F_q, add the precomputed nr fields data, and enumerate curves using HCP
    CC = CurvesClassifier_Fq(p, n, NF=nf)
    CC.enumerate_curves(use_HPC=True, add_SS=False) #here we create curves by j invariant, and generate the twists
    CC.compute_torsion(max_ell=50, compute_volcano=True) #
        #CC.enumerate_curves()
        #print(f"\nTotal unique j-invariants found: {len(added_js)}")
    #CC.compute_volcano_edges()
        #print(f"Saving results for F_{q}...")
    out = CC.toJSON()
    #Data.saveJSON(f"./data/{p}", f"curves_TEST_{q}.json", out, readable=False)

    if to_db:
        db = _cache_db_path(db_path)
        _init_cache_db(db)
        written = _save_catalogue_to_cache_db(p, n, out, db)
        print(f"Saved {written} number-field payload(s) to cache DB: {db}")
        #print(f"Done with F_{q}. {CC.catalogue.size} curves found.\n") 
       
if __name__ == "__main__":
    args = parse_args()
    generate(args.p, args.n, t_list=args.t_list, to_db=args.to_db, db_path=args.db_path)
    

        
        
    
    
    