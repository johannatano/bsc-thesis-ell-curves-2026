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
    p.add_argument("--plist", type=int, nargs="*", default=None, help="List of specific primes to process (for debugging)")
    return p.parse_args()


def enumerate_only(primes, n, use_HCP=False, use_CN=False, q_max=10**20, plist=None, level=None):
    from lib.curves import reset_t0_cache
    cached = []
    for p in primes:
        if plist is not None and p not in plist:
            continue  # Filter to specific primes if plist is provided
        reset_t0_cache()
        nf = None
        if use_HCP or use_CN:
            NFC = NumberFieldsClassifier_Fq(p)
            nf = NFC.generate([n], q_max=q_max)
        q = p**n
        if q > q_max:
            continue

        print(
            f"{Colors.HEADER}--------------{level}------------------------Enumerate F_{q}, p={p}, n={n}{Colors.ENDC}"
        )

        '''if p % 3 == 2:
            print(f"{Colors.FAIL}j-invariant 0 is SS{Colors.ENDC}")

        print(f"gcd(6, q-1) = {gcd(6, q-1)}, p % 3 = {p % 3}, q % 3 = {q % 3}")

        if p % 4 == 3:
            print(f"{Colors.FAIL}j-invariant 1728 is SS{Colors.ENDC}")
        print(f"gcd(4, q-1) = {gcd(4, q-1)}")'''

        CC = CurvesClassifier_Fq(p, n, NF=nf)
        CC.enumerate_curves(
            use_HCP=use_HCP, use_CN=use_CN, add_curves=True, add_SS=True, special_only=False
        )
        cached.append((p, CC))

        CC.check_SS(level=level)  # Sanity check for supersingular curve counts
    return cached


def compute_traces(cached, ell, k, use_CN=False, plist=None):
    ps, Ts, sageTs, diffs, good, infos = [], [], [], [], [], []

    print("COMPUTE TRACES", ell)
    for p, CC in cached:
        if ell > 1 and p % ell == 0:
            continue

        if plist is not None and p not in plist:
            continue  # Filter to specific primes if plist is provided

        q = CC.field.q  # Get the actual field size q = p^n from the classifier
        T, NC, NSS, traces, hk_evals, vals, full_r = CC.compute_hecke(k=k, level=ell, use_CN=use_CN)

        sage_T = 0
        diff = 0

        # if ell < 5 or k == 0:
        sage_T = 0#CuspForms(Gamma1(ell), k + 2).hecke_operator(q).trace()

        diff = T - sage_T
        # T = sage_T
        # else:
        # should be able to compute true trace by diff
        sgn = 0
        if p % ell == 1:
            sgn = 1
        elif p % ell == ell - 1:
            sgn = -1

        # diff = (ell - 1) * (1 + sgn**k) // 2
        # T -= diff

        ps.append(float(q))  # Use q instead of p for plotting
        Ts.append(float(T))
        sageTs.append(float(sage_T))
        diffs.append(float(diff))

        good.append(p % ell == 1)
        infos.append(dict(p=p, q=q, T=T, sage_T=sage_T, diff=diff, NC=NC, NSS=NSS, traces=traces,
                          hk_evals=hk_evals, vals=vals, full_r=full_r,
                          sym=(p % ell == ell - 1)))
        
        print(f"p={p}, q={q}, q equiv ell = {q % ell}, T={T}, sage_T={sage_T}, diff={diff}, NC={NC}, NSS={NSS}, full_r={full_r}")
    return ps, Ts, sageTs, diffs, good, infos


if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    print("\n")
    print("=" * 80 + "")
    print(f"Using rank detection method: {Config.rank_method}")
    print("=" * 80 + "\n")

    primes = [p for p in primerange(5, args.pmax)]

    # Cache for extension degrees (computed on-demand)
    n_list = [1, 2, 3, 4, 5, 6, 7, 8]
    cached_all = {}

    ell_list = list(range(1, 101))
    init_ell_idx = ell_list.index(args.l) if args.l in ell_list else 0
    init_k = args.k

    level = ell_list[init_ell_idx]

    # Only enumerate for the initial n
    init_n = args.n if args.n in n_list else 1
    if args.plist:
        print(f"Enumerating curves for n={init_n}, primes {args.plist}...")
    else:
        print(f"Enumerating curves for n={init_n}, {len(primes)} primes up to {args.pmax}...")
    cached_all[init_n] = enumerate_only(primes, init_n, use_HCP=args.use_hcp, use_CN=args.use_cn, plist=args.plist, level=level)

    
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 10))
    plt.subplots_adjust(bottom=0.15, hspace=0.3)

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

    def fmt_info(info, show_sage=False, show_diff=False, ell=None, k=None, n=None):
        # Truncate long lists to avoid rendering overflow
        max_display = 10
        vals_fmt = [round(float(v), 3) for v in info["vals"]]
        traces = info['traces']
        hk_evals = info['hk_evals']

        traces_str = str(traces[:max_display]) + (f" ... ({len(traces)} total)" if len(traces) > max_display else "")
        hk_evals_str = str(hk_evals[:max_display]) + (f" ... ({len(hk_evals)} total)" if len(hk_evals) > max_display else "")
        vals_str = str(vals_fmt[:max_display]) + (f" ... ({len(vals_fmt)} total)" if len(vals_fmt) > max_display else "")

        if show_sage:
            result = (
                f"p = {info['p']},  q = {info['q']}\n"
                f"Sage T = {info['sage_T']}\n"
                # ,  NC = {info['NC']},  NSS = {info['NSS']},  full_r = {info['full_r']}
                # f"traces    = {traces_str}\n"
                # f"hk_evals  = {hk_evals_str}\n"
                # f"vals      = {vals_str}"
            )
        else:
            p = info['p']
            sgn = 0
            if p % ell == 1:
                sgn = 1
            elif p % ell == ell - 1:
                sgn = -1
            exp_diff = (
                (ell - 1) // 2 * (n - sgn**k*k * p ** (n + k - 1))
                if ell and k is not None and n
                else "N/A"
            )
            result = (
                f"p = {p},  q = {info['q']}\n"
                f"SUM = {info['T']}\n"
                f"EXP DIFF = {exp_diff}\n"
                # ,  NC = {info['NC']},  NSS = {info['NSS']},  full_r = {info['full_r']}
                # f"traces    = {traces_str}\n"
                # f"hk_evals  = {hk_evals_str}\n"
                # f"vals      = {vals_str}"
            )
        if show_diff:
            result += f"\ndiff = {info['diff']}"
        return result

    ps0, Ts0, sageTs0, diffs0, good0, infos0 = compute_traces(cached_all[init_n], ell_list[init_ell_idx], init_k, args.use_cn, plist=args.plist)
    
    print(f"Trace = {Ts0}")
    
    state = {"ps": ps0, "Ts": Ts0, "sageTs": sageTs0, "diffs": diffs0, "infos": infos0, "good": good0, "log": True, "n": init_n, "plist": args.plist}

    def draw_ax(ell, k, n, ps, Ts, sageTs, diffs, good, infos):
        # Top plot: Sage computed traces
        ax1.cla()
        ax1.scatter(ps, sageTs, s=30, color=get_colors(infos, good), zorder=3)
        ax1.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax1.set_xlabel("q", fontsize=13)
        ax1.set_ylabel(f"Sage tr $T_{{{ell}}}$", fontsize=13)
        ax1.set_title(f"Sage Hecke trace $T_{{{ell}}}$, weight $k={k}$, $n={n}$", fontsize=14)
        ax1.set_yscale("symlog" if state["log"] else "linear")
        ax1.grid(True, alpha=0.3)

        # Bottom plot: Computed T and T-diff
        ax2.cla()
        Ts_corrected = [T - d for T, d in zip(Ts, diffs)]
        ax2.scatter(ps, Ts, s=30, color="blue", label="T (computed)", alpha=0.7, zorder=3)
        ax2.scatter(ps, Ts_corrected, s=30, color="orange", label="T - diff", alpha=0.7, zorder=3)
        ax2.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax2.set_xlabel("q", fontsize=13)
        ax2.set_ylabel(f"Trace values", fontsize=13)
        ax2.set_title(f"Computed traces, weight $k={k}$, $n={n}$", fontsize=14)
        ax2.set_yscale("symlog" if state["log"] else "linear")
        ax2.legend(loc="upper right")
        ax2.grid(True, alpha=0.3)

        # Annotations for both plots
        a1 = ax1.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
                        bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gray", alpha=0.95),
                        fontsize=8, fontfamily="monospace", visible=False)
        a2 = ax2.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
                        bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gray", alpha=0.95),
                        fontsize=8, fontfamily="monospace", visible=False)
        state["annot1"] = a1
        state["annot2"] = a2

    draw_ax(ell_list[init_ell_idx], init_k, init_n, ps0, Ts0, sageTs0, diffs0, good0, infos0)

    ax_ell = plt.axes([0.15, 0.09, 0.65, 0.015])
    ax_k   = plt.axes([0.15, 0.06, 0.65, 0.015])
    ax_n   = plt.axes([0.15, 0.03, 0.65, 0.015])
    ax_log = plt.axes([0.82, 0.07, 0.10, 0.03])
    from matplotlib.widgets import Button
    s_ell   = Slider(ax_ell, 'ℓ (index)', 0, len(ell_list) - 1, valinit=init_ell_idx, valstep=1)
    s_k     = Slider(ax_k,   'k',          0, 12,                valinit=init_k,       valstep=1)
    s_n     = Slider(ax_n,   'n',          1, len(n_list),       valinit=init_n,       valstep=1)
    btn_log = Button(ax_log, 'log: ON', color="lightblue")
    ell_label = ax_ell.text(1.02, 0.5, f"ℓ={ell_list[init_ell_idx]}", transform=ax_ell.transAxes, va="center")

    def toggle_log(_):
        state["log"] = not state["log"]
        btn_log.label.set_text("log: ON" if state["log"] else "log: OFF")
        ell = ell_list[int(s_ell.val)]
        k   = int(s_k.val)
        n   = int(s_n.val)
        draw_ax(ell, k, n, state["ps"], state["Ts"], state["sageTs"], state["diffs"], state["good"], state["infos"])
        fig.canvas.draw_idle()

    btn_log.on_clicked(toggle_log)

    def update(_):
        ell = ell_list[int(s_ell.val)]
        k   = int(s_k.val)
        n   = int(s_n.val)

        # Enumerate on-demand if not cached
        if n not in cached_all:
            if args.plist:
                print(f"Enumerating curves for n={n}, primes {args.plist}...")
            else:
                print(f"Enumerating curves for n={n}, {len(primes)} primes up to {args.pmax}...")
            cached_all[n] = enumerate_only(primes, n, use_HCP=args.use_hcp, use_CN=args.use_cn, plist=args.plist)

        ps, Ts, sageTs, diffs, good, infos = compute_traces(cached_all[n], ell, k, args.use_cn, plist=state.get("plist"))
        
        
        state["ps"], state["Ts"], state["sageTs"], state["diffs"], state["infos"], state["good"], state["n"] = ps, Ts, sageTs, diffs, infos, good, n
        draw_ax(ell, k, n, ps, Ts, sageTs, diffs, good, infos)
        ell_label.set_text(f"ℓ={ell}")
        fig.canvas.draw_idle()

    def on_hover(event):
        # Handle hover for both axes
        annot1 = state.get("annot1")
        annot2 = state.get("annot2")

        if event.inaxes == ax1 and annot1 is not None and state["ps"]:
            ps, sageTs, infos = state["ps"], state["sageTs"], state["infos"]
            xlim, ylim = ax1.get_xlim(), ax1.get_ylim()
            xs = (xlim[1] - xlim[0]) or 1
            ys = (ylim[1] - ylim[0]) or 1
            dists = [((event.xdata - px) / xs) ** 2 + ((event.ydata - ty) / ys) ** 2
                     for px, ty in zip(ps, sageTs)]
            idx = min(range(len(dists)), key=lambda i: dists[i])
            if dists[idx] < 2e-4:
                annot1.xy = (ps[idx], sageTs[idx])
                ell = ell_list[int(s_ell.val)]
                k = int(s_k.val)
                n = int(s_n.val)
                annot1.set_text(fmt_info(infos[idx], show_sage=True, ell=ell, k=k, n=n))
                annot1.set_visible(True)
            else:
                annot1.set_visible(False)
            if annot2:
                annot2.set_visible(False)
        elif event.inaxes == ax2 and annot2 is not None and state["ps"]:
            ps, Ts, diffs, infos = state["ps"], state["Ts"], state["diffs"], state["infos"]
            xlim, ylim = ax2.get_xlim(), ax2.get_ylim()
            xs = (xlim[1] - xlim[0]) or 1
            ys = (ylim[1] - ylim[0]) or 1
            # Check both T and T-diff points
            Ts_corrected = [T - d for T, d in zip(Ts, diffs)]
            dists_T = [((event.xdata - px) / xs) ** 2 + ((event.ydata - ty) / ys) ** 2
                       for px, ty in zip(ps, Ts)]
            dists_Tcorr = [((event.xdata - px) / xs) ** 2 + ((event.ydata - ty) / ys) ** 2
                           for px, ty in zip(ps, Ts_corrected)]
            idx_T = min(range(len(dists_T)), key=lambda i: dists_T[i])
            idx_Tcorr = min(range(len(dists_Tcorr)), key=lambda i: dists_Tcorr[i])

            ell = ell_list[int(s_ell.val)]
            k = int(s_k.val)
            n = int(s_n.val)
            if dists_T[idx_T] < dists_Tcorr[idx_Tcorr] and dists_T[idx_T] < 2e-4:
                annot2.xy = (ps[idx_T], Ts[idx_T])
                annot2.set_text(f"T (computed)\n" + fmt_info(infos[idx_T], show_diff=True, ell=ell, k=k, n=n))
                annot2.set_visible(True)
            elif dists_Tcorr[idx_Tcorr] < 2e-4:
                annot2.xy = (ps[idx_Tcorr], Ts_corrected[idx_Tcorr])
                annot2.set_text(f"T - diff = {Ts_corrected[idx_Tcorr]} (diff = {diffs[idx_Tcorr]})\n" + fmt_info(infos[idx_Tcorr], show_diff=False, ell=ell, k=k, n=n))
                annot2.set_visible(True)
            else:
                annot2.set_visible(False)
            if annot1:
                annot1.set_visible(False)
        else:
            if annot1:
                annot1.set_visible(False)
            if annot2:
                annot2.set_visible(False)
        fig.canvas.draw_idle()
    
    fig.canvas.mpl_connect("motion_notify_event", on_hover)
    s_ell.on_changed(update)
    s_k.on_changed(update)
    s_n.on_changed(update)
    plt.show()
