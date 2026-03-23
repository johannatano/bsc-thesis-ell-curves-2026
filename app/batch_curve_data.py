import argparse
import json
import sqlite3
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "data"
DEFAULT_DB = DATA_ROOT / "cache.sqlite3"
DEBUG_DUMP = ROOT / "out.json"


@dataclass(frozen=True)
class Job:
    p: int
    n: int
    discriminants: list[int]
    traces: list[int]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch prepopulate curve_cache by invoking generate_curve_data.py per (p,n,D)."
    )
    parser.add_argument(
        "-p",
        "--p",
        type=int,
        nargs="*",
        default=None,
        help="Prime list to process (default: all primes found in data/*/nr_fields.json)",
    )
    parser.add_argument(
        "-n",
        "--n",
        type=int,
        nargs="*",
        default=None,
        help="Optional n filter",
    )
    parser.add_argument(
        "-D",
        "--discriminants",
        type=int,
        nargs="*",
        default=None,
        help="Optional discriminant filter",
    )
    parser.add_argument(
        "--db-path",
        type=str,
        default=str(DEFAULT_DB),
        help="SQLite cache path (default: ./data/cache.sqlite3)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=900,
        help="Per-job timeout in seconds (default: 900)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Regenerate even if cache row already exists",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Only print planned jobs; do not run generation",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        default=False,
        help="Stop batch at first failed job",
    )
    return parser.parse_args()


def init_cache_db(db_path: Path) -> None:
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


def cache_has_row(db_path: Path, p: int, n: int, d: int) -> bool:
    with sqlite3.connect(db_path) as conn:
        row = conn.execute(
            "SELECT 1 FROM curve_cache WHERE p=? AND n=? AND d=? LIMIT 1",
            (int(p), int(n), int(d)),
        ).fetchone()
    return row is not None


def cache_has_rows_for_job(db_path: Path, job: Job) -> bool:
    if not job.discriminants:
        return False

    placeholders = ",".join("?" for _ in job.discriminants)
    params = [int(job.p), int(job.n), *[int(d) for d in job.discriminants]]
    with sqlite3.connect(db_path) as conn:
        row = conn.execute(
            f"SELECT COUNT(DISTINCT d) FROM curve_cache WHERE p=? AND n=? AND d IN ({placeholders})",
            params,
        ).fetchone()
    return int(row[0] or 0) == len(job.discriminants)


def discover_primes() -> list[int]:
    primes = []
    if not DATA_ROOT.exists():
        return primes
    for child in DATA_ROOT.iterdir():
        if not child.is_dir():
            continue
        if not child.name.isdigit():
            continue
        if (child / "nr_fields.json").exists():
            primes.append(int(child.name))
    return sorted(primes)


def load_jobs_for_prime(p: int, n_filter: set[int] | None, d_filter: set[int] | None) -> list[Job]:
    path = DATA_ROOT / str(p) / "nr_fields.json"
    if not path.exists():
        print(f"[skip] no nr_fields.json for p={p}: {path}")
        return []

    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    grouped: dict[int, dict[str, set[int]]] = {}
    for nf in data.get("nr_fields", {}).get("nf", []):
        d = int(nf.get("D"))
        if d_filter is not None and d not in d_filter:
            continue

        for tree in nf.get("tree", []):
            n = int(tree.get("n"))
            if n_filter is not None and n not in n_filter:
                continue

            traces = [abs(int(ic.get("t"))) for ic in (tree.get("I_t", []) or []) if "t" in ic]
            if not traces:
                continue

            if n not in grouped:
                grouped[n] = {"discriminants": set(), "traces": set()}
            grouped[n]["discriminants"].add(d)
            grouped[n]["traces"].update(traces)

    jobs: list[Job] = [
        Job(
            p=p,
            n=n,
            discriminants=sorted(grouped[n]["discriminants"], key=lambda d: (abs(d), d)),
            traces=sorted(grouped[n]["traces"]),
        )
        for n in sorted(grouped.keys())
    ]
    jobs.sort(key=lambda j: (j.p, j.n))
    return jobs


def write_debug_dump_from_db(job: Job, db_path: Path) -> bool:
    if not job.discriminants:
        return False

    placeholders = ",".join("?" for _ in job.discriminants)
    with sqlite3.connect(db_path) as conn:
        rows = conn.execute(
            f"SELECT d, payload_json FROM curve_cache WHERE p=? AND n=? AND d IN ({placeholders})",
            [int(job.p), int(job.n), *[int(d) for d in job.discriminants]],
        ).fetchall()

    if not rows:
        return False

    number_fields = []
    for _, payload_json in sorted(rows, key=lambda row: (abs(int(row[0])), int(row[0]))):
        if not payload_json:
            continue
        payload = json.loads(payload_json)
        number_fields.extend(payload.get("catalogue", {}).get("number_fields", []) or [])

    payload = {
        "char": int(job.p),
        "catalogue": {
            "number_fields": number_fields,
        },
    }
    print(f"Writing DB payload debug dump to {DEBUG_DUMP}")
    with DEBUG_DUMP.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    return True


def run_job(job: Job, db_path: Path, timeout: int) -> tuple[bool, str]:
    
    print(f"Running job: p={job.p}, n={job.n}, D_count={len(job.discriminants)}, traces={job.traces}")
    cmd = [
        "sage",
        "-python",
        "generate_curve_data.py",
        "-p",
        str(job.p),
        "-n",
        str(job.n),
        "-t",
        *[str(t) for t in job.traces],
        "--to-db",
        "--db-path",
        str(db_path),
    ]

    start = time.time()
    result = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.time() - start
    if result.returncode == 0:
        write_debug_dump_from_db(job, db_path)

    key = f"p={job.p} n={job.n} D_count={len(job.discriminants)}"
    if result.returncode == 0:
        return True, f"[ok] {key} ({elapsed:.1f}s)"

    tail_out = (result.stdout or "")[-1200:]
    tail_err = (result.stderr or "")[-1200:]
    msg = (
        f"[fail] {key} (code={result.returncode}, {elapsed:.1f}s)\n"
        f"stdout_tail:\n{tail_out}\n"
        f"stderr_tail:\n{tail_err}"
    )
    return False, msg


def main():
    args = parse_args()
    db_path = Path(args.db_path).resolve()
    init_cache_db(db_path)

    primes = sorted(set(args.p)) if args.p else discover_primes()
    if not primes:
        print("No primes selected/found. Nothing to do.")
        return 0

    n_filter = set(args.n) if args.n else None
    d_filter = set(args.discriminants) if args.discriminants else None

    jobs: list[Job] = []
    for p in primes:
        jobs.extend(load_jobs_for_prime(p, n_filter=n_filter, d_filter=d_filter))

    if not jobs:
        print("No jobs found for provided filters.")
        return 0

    planned = len(jobs)
    print(f"Planned jobs: {planned}")

    if args.dry_run:
        for j in jobs:
            print(f"[plan] p={j.p} n={j.n} D_count={len(j.discriminants)} traces={j.traces}")
        return 0

    ok_count = 0
    skip_count = 0
    fail_count = 0

    started = time.time()
    for idx, job in enumerate(jobs, start=1):
        key = f"p={job.p} n={job.n} D_count={len(job.discriminants)}"
        if not args.force and cache_has_rows_for_job(db_path, job):
            skip_count += 1
            write_debug_dump_from_db(job, db_path)
            print(f"[{idx}/{planned}] [skip-cache] {key}")
            continue

        print(f"[{idx}/{planned}] [run] {key} traces={job.traces}")
        ok, message = run_job(job, db_path=db_path, timeout=args.timeout)
        print(message)

        if ok:
            ok_count += 1
        else:
            fail_count += 1
            if args.stop_on_error:
                break

    elapsed = time.time() - started
    print("=" * 72)
    print(
        f"Done in {elapsed:.1f}s | ok={ok_count} skip={skip_count} fail={fail_count} "
        f"total={planned} db={db_path}"
    )
    print("=" * 72)

    return 1 if fail_count > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
    





