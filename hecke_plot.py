import argparse
import sys
import os
import time

from lib.curves_classifier import CurvesClassifier_Fq
from lib.nr_fields_classifier import NumberFieldsClassifier_Fq

from utils.common import Logger, Colors, Data, Config
from sympy import primerange
from sage.all import *
import numpy as np
import requests
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument(
        "-p", "--p", type=int, required=False, default=-1, help="Field char p"
    )
    p.add_argument(
        "-n",
        "--n",
        type=int,
        required=False,
        default=1,
        help="Field extension degree n",
    )
    p.add_argument("-l", "--l", type=int, required=False, default=-1, help="Level ℓ")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")
    p.add_argument(
        "--use-hcp",
        action="store_true",
        default=False,
        help="Use HCP (Hilbert class polynomial) enumeration instead of direct method",
    )
    p.add_argument(
        "--use-cn",
        action="store_true",
        default=False,
        help="Use Class Numbers ie no j invariants instead of direct method",
    )
    p.add_argument(
        "--rank-method",
        choices=["auto", "div_poly", "mod_poly", "invariants"],
        default="mod_poly",
        help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)",
    )
    p.add_argument(
        "--true-height",
        action="store_true",
        default=False,
        help="Use exact BFS height in isogeny volcano instead of floor test",
    )
    p.add_argument(
        "--interactive",
        action="store_true",
        default=False,
        help="Show interactive slider UI (caches enumeration, recomputes Hecke on k/ell change)",
    )
    p.add_argument("--pmax", type=int, default=100, help="Upper prime bound for interactive mode")
    return p.parse_args()


def enumerate_only(primes, n, use_HCP=False, use_CN=False, q_max=10**20):
    from lib.curves import reset_t0_cache
    cached = []
    for p in primes:
        reset_t0_cache()
        nf = None
        if use_HCP or use_CN:
            NFC = NumberFieldsClassifier_Fq(p)
            nf = NFC.generate([n], q_max=q_max)
        q = p**n
        if q > q_max:
            continue
        CC = CurvesClassifier_Fq(p, n, NF=nf)
        CC.enumerate_curves(use_HCP=use_HCP, use_CN=use_CN, add_curves=True, add_SS=True)
        cached.append((p, CC))
        print(f"Enumerated F_{p}")
    return cached


def compute_traces(cached, ell, k, use_CN=False):
    ps, Ts, good, infos = [], [], [], []
    for p, CC in cached:
        if p % ell == 0:
            continue
        T, NC, NSS, traces, hk_evals, vals, full_r = CC.compute_hecke(k=k, level=ell, use_CN=use_CN)
        ps.append(float(p))
        Ts.append(float(T))
        good.append(p % ell == 1)
        infos.append(dict(p=p, T=T, NC=NC, NSS=NSS, traces=traces,
                          hk_evals=hk_evals, vals=vals, full_r=full_r,
                          sym=(p % ell == ell - 1)))
    return ps, Ts, good, infos


def run(p: int, l: int, k: int, n: int, use_HCP=False, use_CN=False):
    primes = list(primerange(5, 50)) if p == -1 else [p]
    # primes = list(primerange(10**6, 10**6 + 100)) if p == -1 else [p]
    # 1000033
    # 10093
    # p=1091
    # 100043
    # 1000003
    # 10050013, 10050017, 10050023, 10050049, 10050059, 10050071, 10050083, 10050101, 10050133, 10050137, 10050167, 10050181, 10050191, 10050197, 10050203, 10050217, 10050223, 10050233, 10050253, 10050283, 10050317, 10050319, 10050331, 10050353, 10050367, 10050377, 10050389, 10050407, 10050413, 10050419, 10050427, 10050437, 10050463, 10050493
    p_powers = [n]
    dsize = len(primes)
    q_max = 10**20
    levels = [l] if l != -1 else list(primerange(2, 100))
    diffs = {}  # (ell, k) -> list of diffs per prime p

    good_p = []
    sym_p = []

    full_r_p = []

    for i in range(dsize):
        p = primes[i]

        # Clear t=0 claim cache to avoid cross-contamination between different p values
        nf = None
        if use_HCP or use_CN:
            NFC = NumberFieldsClassifier_Fq(p)
            nf = NFC.generate(p_powers, q_max=q_max)
        q = p**n
        if q > q_max:
            print(f"Skipping F_{q} due to size > {q_max}")
            continue
        CC = CurvesClassifier_Fq(p, n, NF=nf)
        CC.enumerate_curves(
            use_HCP=use_HCP, use_CN=use_CN, add_curves=True, add_SS=True
        )

        for ell in levels:

            if p % ell == 0:
                print(f"{Colors.WARNING}Skipping level ℓ={ell} for p={p} as it divides the field characteristic{Colors.ENDC}")
                continue
            # if p % ell == 1:
            #    print(f"{Colors.GREEN}CAUTION, we have p ≡ 1 (mod {ell}), p={p}{Colors.ENDC}")
            """print(
                f"\n{Colors.BLUE}=== Computing Hecke operator T_{ell} for weight {k} ==={Colors.ENDC}\n"
            )"""
            T, NC, NSS, traces, hk_evals, vals, full_r = CC.compute_hecke(
                k=args.k, level=ell, use_CN=use_CN
            )

            if full_r:
                full_r_p.append(p)

            sgn = 0
            if p % ell == 1:
                sgn = 1
                good_p.append(p)
            elif p % ell == ell - 1:
                sgn = -1
                sym_p.append(p)

            diff_comp = (ell-1)*(1+sgn**k) // 2
            EP = CC.count_EP(ell, use_CN=use_CN)

            trace_val = 0
            diff = 0

            if ell < 5 or k == 0:
                trace_val = CuspForms(Gamma1(ell), k + 2).hecke_operator(q).trace()
                diff = T - trace_val
                T = trace_val
            else:
                T -= diff_comp

            '''if diff_comp != diff and k > 0 and ell >= 5:
                print(
                    f"{Colors.FAIL}Discrepancy for p={p}, ell={ell}, k={k}: computed diff {diff_comp} does not match expected {diff}{Colors.ENDC}"
                )'''

            '''if ell < 5:
                T = trace_val
            else:
                T -= diff_comp'''

            diffs.setdefault((ell, k, (ell - 1) // 2), []).append(
                (q, diff, T, trace_val, NC, NSS, vals, traces, hk_evals)
            )

            print(
                f"{Colors.HEADER}p={p}, Total Hecke trace for level {ell} and weight {args.k}: {T}, sage trace: {trace_val}, difference: {T - trace_val}, EP: {EP}, NUM CURVES: {NC}, NUM SS: {NSS}, full_r: {full_r}{Colors.ENDC}"
            )
    print("\n" + "=" * 80)

    print(f"Good p (p ≡ 1 mod ell): {good_p}")
    print(f"Symmetric p (p ≡ -1 mod ell): {sym_p}")
    print(f"p with rank 2: {full_r_p}")
    print("=" * 80)
    for key, entries in sorted(diffs.items()):
        ell, k_, dim = key
        for (
            p_val,
            d,
            T_val,
            trace_val,
            NC_val,
            NSS_val,
            vals,
            traces_val,
            hk_evals_val,
        ) in entries:
            if p_val % ell == 1:
                color = Colors.GREEN
                label = " [q≡1 mod ell]"
            elif p_val % ell == ell - 1:
                color = Colors.WARNING if k_ % 2 == 1 else Colors.GREEN
                label = " [q≡-1 mod ell]"
            else:
                color = Colors.ENDC
                label = ""
            # print(
            #    f"    {color}q={p_val}: diff={d}, T={T_val}, sage_trace={trace_val}, vals={vals}, NC={NC_val}, NSS={NSS_val}, traces={traces_val}, hk_evals={hk_evals_val}, {label} [q≡{p_val % ell} mod {ell}]{Colors.ENDC}"
            # )
    return diffs, good_p, sym_p


if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    print("\n")
    print("=" * 80 + "")
    print(f"Using rank detection method: {Config.rank_method}")
    print("=" * 80 + "\n")
    start_hcp = time.time()
    diffs, good_p, sym_p = run(args.p, args.l, args.k, args.n, use_HCP=args.use_hcp, use_CN=args.use_cn)
    end_hcp = time.time()
    print(
        f"{Colors.HEADER}Hecke Trace computed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}"
    )

    if args.interactive:
        primes = [p for p in primerange(5, args.pmax)]
        print(f"Enumerating curves for {len(primes)} primes up to {args.pmax}...")
        cached = enumerate_only(primes, args.n, use_HCP=args.use_hcp, use_CN=args.use_cn)

        ell_list = list(primerange(2, 100))
        init_ell_idx = ell_list.index(args.l) if args.l in ell_list else 0
        init_k = args.k

        fig, ax = plt.subplots(figsize=(11, 6))
        plt.subplots_adjust(bottom=0.22)

        def get_colors(infos, good):
            colors = []
            for info, g in zip(infos, good):
                if info["NC"] == 0:
                    colors.append("red")
                elif g:
                    colors.append("green")
                elif info["sym"]:
                    colors.append("yellow")
                elif float(info["T"]) == 0 and info["NC"] != 0:
                    colors.append("turquoise")
                else:
                    colors.append("black")
            return colors

        def fmt_info(info):
            vals_fmt = [round(float(v), 3) for v in info["vals"]]
            return (
                f"p = {info['p']}\n"
                f"T = {info['T']},  NC = {info['NC']},  NSS = {info['NSS']},  full_r = {info['full_r']}\n"
                f"traces    = {info['traces']}\n"
                f"hk_evals  = {info['hk_evals']}\n"
                f"vals      = {vals_fmt}"
            )

        ps0, Ts0, good0, infos0 = compute_traces(cached, ell_list[init_ell_idx], init_k, args.use_cn)
        state = {"ps": ps0, "Ts": Ts0, "infos": infos0, "good": good0, "log": True}

        def draw_ax(ell, k, ps, Ts, good, infos):
            ax.cla()
            ax.scatter(ps, Ts, s=30, color=get_colors(infos, good), zorder=3)
            ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
            ax.set_xlabel("p", fontsize=13)
            ax.set_ylabel(f"tr $T_{{{ell}}}$", fontsize=13)
            ax.set_title(f"Hecke trace $T_{{{ell}}}$, weight $k={k}$", fontsize=14)
            ax.set_yscale("symlog" if state["log"] else "linear")
            ax.grid(True, alpha=0.3)
            a = ax.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
                            bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gray", alpha=0.95),
                            fontsize=8, fontfamily="monospace", visible=False)
            state["annot"] = a

        draw_ax(ell_list[init_ell_idx], init_k, ps0, Ts0, good0, infos0)

        ax_ell = plt.axes([0.15, 0.10, 0.65, 0.03])
        ax_k   = plt.axes([0.15, 0.05, 0.65, 0.03])
        ax_log = plt.axes([0.82, 0.10, 0.10, 0.04])
        from matplotlib.widgets import Button
        s_ell   = Slider(ax_ell, 'ℓ (index)', 0, len(ell_list) - 1, valinit=init_ell_idx, valstep=1)
        s_k     = Slider(ax_k,   'k',          0, 12,                valinit=init_k,       valstep=1)
        btn_log = Button(ax_log, 'log: ON', color="lightblue")
        ell_label = ax_ell.text(1.02, 0.5, f"ℓ={ell_list[init_ell_idx]}", transform=ax_ell.transAxes, va="center")

        def toggle_log(_):
            state["log"] = not state["log"]
            btn_log.label.set_text("log: ON" if state["log"] else "log: OFF")
            ell = ell_list[int(s_ell.val)]
            k   = int(s_k.val)
            draw_ax(ell, k, state["ps"], state["Ts"], state["good"], state["infos"])
            fig.canvas.draw_idle()

        btn_log.on_clicked(toggle_log)

        def update(_):
            ell = ell_list[int(s_ell.val)]
            k   = int(s_k.val)
            ps, Ts, good, infos = compute_traces(cached, ell, k, args.use_cn)
            state["ps"], state["Ts"], state["infos"], state["good"] = ps, Ts, infos, good
            draw_ax(ell, k, ps, Ts, good, infos)
            ell_label.set_text(f"ℓ={ell}")
            fig.canvas.draw_idle()

        def on_hover(event):
            annot = state.get("annot")
            if annot is None or event.inaxes != ax or not state["ps"]:
                if annot:
                    annot.set_visible(False)
                fig.canvas.draw_idle()
                return
            ps, Ts, infos = state["ps"], state["Ts"], state["infos"]
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            xs = (xlim[1] - xlim[0]) or 1
            ys = (ylim[1] - ylim[0]) or 1
            dists = [((event.xdata - px) / xs) ** 2 + ((event.ydata - ty) / ys) ** 2
                     for px, ty in zip(ps, Ts)]
            idx = min(range(len(dists)), key=lambda i: dists[i])
            if dists[idx] < 2e-4:
                annot.xy = (ps[idx], Ts[idx])
                annot.set_text(fmt_info(infos[idx]))
                annot.set_visible(True)
            else:
                annot.set_visible(False)
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", on_hover)
        s_ell.on_changed(update)
        s_k.on_changed(update)
        plt.show()
    else:
        for (ell, k_, dim), entries in sorted(diffs.items()):
            ps = [float(q) for (q, d, T, tv, NC, NSS, vals, tr, hk) in entries]
            Ts = [float(T) for (q, d, T, tv, NC, NSS, vals, tr, hk) in entries]
            colors = ["green" if p in good_p else "red" if t == 0 else "steelblue" for p, t in zip(ps, Ts)]
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.scatter(ps, Ts, s=30, color=colors, zorder=3)
            ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
            ax.set_xlabel("p", fontsize=13)
            ax.set_ylabel(f"tr $T_{{{ell}}}$", fontsize=13)
            ax.set_title(f"Hecke trace $T_{{{ell}}}$, weight $k={k_}$", fontsize=14)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()
