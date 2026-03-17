# app.py
from flask import Flask, jsonify, render_template, abort,request
from pathlib import Path
from sage.all import GF, EllipticCurve, ZZ
import json
import subprocess
import sqlite3
import time

app = Flask(__name__)
DATA_ROOT = Path(__file__).parent / "data"   # Points to WSL.BscThesis/data
CACHE_DB = Path(__file__).parent / "data" / "cache.sqlite3"


def _get_db_conn():
    conn = sqlite3.connect(CACHE_DB)
    conn.row_factory = sqlite3.Row
    return conn


def init_cache_db():
    CACHE_DB.parent.mkdir(parents=True, exist_ok=True)
    with _get_db_conn() as conn:
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


def upsert_test_cache_payload(p: int, n: int, d: int, payload: dict) -> dict:
    now = int(time.time())
    payload_json = json.dumps(payload, separators=(",", ":"))
    with _get_db_conn() as conn:
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
        row = conn.execute(
            "SELECT p, n, d, status, payload_json, created_at, updated_at FROM curve_cache WHERE p=? AND n=? AND d=?",
            (int(p), int(n), int(d)),
        ).fetchone()
    return {
        "p": int(row["p"]),
        "n": int(row["n"]),
        "D": int(row["d"]),
        "status": row["status"],
        "payload": json.loads(row["payload_json"]),
        "created_at": int(row["created_at"]),
        "updated_at": int(row["updated_at"]),
    }


def _load_traces_for_field(p: int, n: int, d: int):
    """Load trace list t from nr_fields.json for a specific (p,n,D)."""
    path = DATA_ROOT / str(p) / "nr_fields.json"
    if not path.exists():
        return []

    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    nf_list = data.get("nr_fields", {}).get("nf", [])
    for nf in nf_list:
        if int(nf.get("D")) != int(d):
            continue
        for tree in nf.get("tree", []):
            if int(tree.get("n", -1)) != int(n):
                continue
            return [int(ic.get("t")) for ic in (tree.get("I_t", []) or []) if "t" in ic]
    return []


def _get_cached_curve_payload(p: int, n: int, d: int):
    init_cache_db()
    with _get_db_conn() as conn:
        row = conn.execute(
            "SELECT payload_json FROM curve_cache WHERE p=? AND n=? AND d=?",
            (int(p), int(n), int(d)),
        ).fetchone()
    if not row:
        return None
    try:
        return json.loads(row["payload_json"])
    except Exception:
        return None


def _generate_and_cache_field_curves(p: int, n: int, d: int):
    traces = _load_traces_for_field(p, n, d)
    if not traces:
        return False, {
            "success": False,
            "error": f"No traces found for p={p}, n={n}, D={d}"
        }, 404

    init_cache_db()

    cmd = [
        "sage", "-python", "generate_curve_data.py",
        "-p", str(p),
        "-n", str(n),
        "-t", *[str(t) for t in traces],
        "--to-db",
        "--db-path", str(CACHE_DB),
    ]

    result = subprocess.run(
        cmd,
        cwd=str(Path(__file__).parent),
        capture_output=True,
        text=True,
        timeout=300,
    )

    if result.returncode != 0:
        return False, {
            "success": False,
            "error": "curve generation failed",
            "returncode": result.returncode,
            "stderr": result.stderr[-4000:],
            "stdout": result.stdout[-1000:],
        }, 500

    return True, {
        "success": True,
        "message": "Generation completed and cached",
        "key": {"p": p, "n": n, "D": d},
        "traces": traces,
    }, 200


def vec_to_F(v, F, a):
    # v = [v0, v1, ..., v_{n-1}] (coeffs mod p), returns sum v_i * a^i
    s = F(0)
    pwr = F(1)
    for c in v:
        s += F(int(c)) * pwr
        pwr *= a
    return s

def F_to_vec(x, n):
    # coefficients in the {1,a,a^2,...} power basis, padded to length n
    coeffs = list(x.polynomial())  # ascending powers
    coeffs += [0] * (n - len(coeffs))
    return [int(c) for c in coeffs[:n]]

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/cache/test_generate", methods=["GET"])
def cache_test_generate():
    """Hello-world endpoint for on-demand cache generation in SQLite.

    Query params:
      - p (default: 5)
      - n (default: 1)
      - D (default: -4)
    """
    try:
        p = int(request.args.get("p", 5))
        n = int(request.args.get("n", 1))
        d = int(request.args.get("D", -4))
    except Exception as e:
        abort(400, f"Bad query params: {e}")

    init_cache_db()

    # Basic fake payload to prove out on-demand write path.
    fake_payload = {
        "hello": "world",
        "kind": "test-cache-payload",
        "generated_for": {"p": p, "n": n, "D": d},
        "generated_at": int(time.time()),
        "sample": {
            "number_fields": [
                {
                    "D": d,
                    "tree": [
                        {
                            "n": n,
                            "I_t": [
                                {
                                    "t": 0,
                                    "O": [{"f": 1, "cn": 1, "D": d}]
                                }
                            ]
                        }
                    ]
                }
            ]
        }
    }

    stored = upsert_test_cache_payload(p, n, d, fake_payload)
    return jsonify({"success": True, "cache": stored})


@app.route("/cache/test_get", methods=["GET"])
def cache_test_get():
    """Read back cached test payload from SQLite by (p,n,D)."""
    try:
        p = int(request.args.get("p", 5))
        n = int(request.args.get("n", 1))
        d = int(request.args.get("D", -4))
    except Exception as e:
        abort(400, f"Bad query params: {e}")

    init_cache_db()
    with _get_db_conn() as conn:
        row = conn.execute(
            "SELECT p, n, d, status, payload_json, created_at, updated_at FROM curve_cache WHERE p=? AND n=? AND d=?",
            (int(p), int(n), int(d)),
        ).fetchone()

    if row is None:
        return jsonify({"success": False, "error": "cache miss", "key": {"p": p, "n": n, "D": d}}), 404

    return jsonify({
        "success": True,
        "cache": {
            "p": int(row["p"]),
            "n": int(row["n"]),
            "D": int(row["d"]),
            "status": row["status"],
            "payload": json.loads(row["payload_json"]),
            "created_at": int(row["created_at"]),
            "updated_at": int(row["updated_at"]),
        }
    })


@app.route("/cache/generate_field_curves", methods=["GET"])
def cache_generate_field_curves():
    """Generate curves for selected field (p,n,D) and persist into SQLite cache.

    This triggers generate_curve_data.py with the trace list loaded from nr_fields.json.
    """
    try:
        p = int(request.args.get("p"))
        n = int(request.args.get("n"))
        d = int(request.args.get("D"))
    except Exception as e:
        abort(400, f"Bad query params: {e}")

    ok, payload, status = _generate_and_cache_field_curves(p, n, d)
    if not ok:
        return jsonify(payload), status

    # Return cache row count for this (p,n,D) as a quick verification.
    with _get_db_conn() as conn:
        row = conn.execute(
            "SELECT COUNT(1) AS cnt FROM curve_cache WHERE p=? AND n=? AND d=?",
            (int(p), int(n), int(d)),
        ).fetchone()

    return jsonify({
        "success": True,
        "message": "Generation completed and cached",
        "key": {"p": p, "n": n, "D": d},
        "traces": payload.get("traces", []),
        "cache_rows_for_key": int(row["cnt"]) if row else 0,
    })

@app.route("/curves/<int:p>/<int:n>")
def curves(p, n):
    # Load curves data from data/<p>/curves_{q}.json
    path = DATA_ROOT / str(p) / f"curves_{p**n}.json"
    if path.exists():
        with path.open("r", encoding="utf-8") as f:
            return jsonify(json.load(f))
    abort(404, description="curves.json not found")

@app.route("/fields/<int:p>")
def fields(p):
    # Load number fields data
    path = DATA_ROOT / str(p) / "nr_fields.json"
    if path.exists():
        with path.open("r", encoding="utf-8") as f:
            return jsonify(json.load(f))
    abort(404, description="nr_fields.json not found")


@app.route('/get_curves/<int:p>/<int:n>', methods=['GET'])
def enumerate_order(p, n):
    """Load curves for a specific number field D over F_{p^n}.

    Cache-first: return from SQLite if present; otherwise generate via
    generate_curve_data.py (for traces belonging to this D) and then return cache.
    """
    try:
        D = int(request.args.get('D'))

        cached = _get_cached_curve_payload(p, n, D)
        if cached is None:
            ok, payload, status = _generate_and_cache_field_curves(p, n, D)
            if not ok:
                return jsonify(payload), status
            cached = _get_cached_curve_payload(p, n, D)

        if cached is None:
            return jsonify({
                'success': False,
                'error': 'Generated data was not found in cache afterwards'
            }), 500

        return jsonify({
            'success': True,
            'catalogue': cached.get('catalogue', {})
        })
        
    except Exception as e:
        import traceback
        print(f"[DEBUG] Exception: {traceback.format_exc()}")
        return jsonify({'success': False, 'error': str(e), 'trace': traceback.format_exc()}), 500

@app.post("/torsion_points")
def torsion_points():
    """
    Request JSON:
    {
      "p": 5, "n": 3, "l": 7,
      "A": [a0,...], "B": [b0,...],
      "generators": [ {"x":[...], "y":[...]}, {"x":[...], "y":[...]} ]  # second optional
    }
    Returns:
    {
      "points": [ {"x":[...], "y":[...]}, ... ],   # coefficient vectors length n
      "has_infinity": true/false
    }
    """
    data = request.get_json(force=True, silent=False)
    try:
        p  = int(data["p"])
        n  = int(data["n"])
        l  = int(data["l"])
        A  = list(map(int, data["A"]))
        B  = list(map(int, data["B"]))
        gens = data.get("generators", [])
    except Exception as e:
        abort(400, f"Bad payload: {e}")

    if l <= 0 or not gens:
        abort(400, "Need positive l and at least one generator")

    # Finite field and curve
    F = GF(p**n, name='a')
    a = F.gen()
    Ael = vec_to_F(A, F, a)
    Bel = vec_to_F(B, F, a)
    E = EllipticCurve(F, [0, 0, 0, Ael, Bel])   # y^2 = x^3 + A x + B

    # Parse generators
    def to_point(g):
        x = vec_to_F(g["x"], F, a)
        y = vec_to_F(g["y"], F, a)
        return E(x, y)
    P = to_point(gens[0])
    Q = to_point(gens[1]) if len(gens) > 1 else None

    pts = set()
    add_inf = False

    if Q is None:
        # cyclic: <P> of size dividing l
        for k in range(l):
            R = k * P
            if R == E(0):
                add_inf = True
            else:
                pts.add(R)
    else:
        # product: <P,Q> ~ Z_l × Z_l (typical), compute grid
        for i in range(l):
            for j in range(l):
                R = i * P + j * Q
                if R == E(0):
                    add_inf = True
                else:
                    pts.add(R)

    out = []
    for R in pts:
        xvec = F_to_vec(R[0], n)
        yvec = F_to_vec(R[1], n)
        out.append({"x": xvec, "y": yvec})

    return jsonify({"points": out, "has_infinity": add_inf})

if __name__ == "__main__":
    app.run(debug=True)
