#!/usr/bin/env python3
"""One diagram per species, showing only the levels that decay faster than a threshold (5 ns).

A species' full level list is a wall -- Ca I alone has 780 levels and 9% of them carry a
lifetime -- so the whole list is left in the CSV and the figure shows only what the project
can use: the fast end. Layout is two columns, even parity left and odd parity right, because
that is the structure the ladder rules use (an m=2 hop conserves parity, m=3 flips it).

Colour answers one question -- is this level inside the 1 - 1000 ps target window -- because
the lifetime is printed beside every level anyway. Levels of one fine-structure multiplet
share a label; the arrow gives the strongest decay channel, which is what decides whether the
state is a VUV/EUV emitter.

Usage:
    python3 scripts/diagram.py                 # every scraped species, threshold 5 ns
    python3 scripts/diagram.py "Li III"        # one species
    python3 scripts/diagram.py --tau-max 1ns   # a different threshold
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"

SUB_NS_LO, SUB_NS_HI = 1e-12, 1e-9    # the picosecond decade: a sanity floor at 1 ps, and a
                                      # ceiling at 1 ns where an ordinary ns laser would do

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8a85"
ACCENT = "#2a78d6"      # inside the target window
NEUTRAL = "#7d7d77"     # outside it, but still under the threshold
LIMIT = "#a8a8a2"
AUTO = "#f1efea"
GRID = "#e9e9e5"

SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
UNITS = {"ps": 1e-12, "ns": 1e-9, "us": 1e-6, "µs": 1e-6, "ms": 1e-3, "s": 1.0}


# --------------------------------------------------------------------------- formatting

def fmt_tau(t: float) -> str:
    for scale, unit in ((1e-12, "ps"), (1e-9, "ns"), (1e-6, "µs"), (1e-3, "ms"), (1.0, "s")):
        if t < scale * 1000:
            v = t / scale
            return f"{v:.0f} {unit}" if v >= 100 else (f"{v:.1f} {unit}" if v >= 10 else f"{v:.2f} {unit}")
    return f"{t:.0f} s"


def fmt_tau_range(taus) -> str:
    lo, hi = min(taus), max(taus)
    if lo <= 0 or (hi - lo) / lo < 0.08:
        return fmt_tau(sum(taus) / len(taus))
    a, b = fmt_tau(lo), fmt_tau(hi)
    return f"{a.split()[0]}–{b}" if a.split()[-1] == b.split()[-1] else f"{a} – {b}"


def fmt_j(j) -> str:
    if j is None:
        return ""
    return str(int(round(j))) if abs(j - round(j)) < 1e-6 else f"{int(round(j * 2))}/2"


def term_mathtext(term: str, js) -> str:
    t = term.strip().rstrip("?").replace("*", "")
    odd = "*" in term
    m = re.match(r"^\s*(?:[a-z]\s*)?(\d+)(\[[\d/]+\]|[A-Z])\s*$", t)
    jj = ",".join(fmt_j(j) for j in sorted({j for j in js if j is not None}))
    sub = f"_{{{jj}}}" if jj else ""
    if not m:
        return f"$\\mathrm{{{t}}}{sub}$"
    mult, letter = m.groups()
    core = letter if letter.startswith("[") else f"\\mathrm{{{letter}}}"
    deg = "^{\\circ}" if odd else ""
    return f"$^{{{mult}}}{core}{deg}{sub}$"


def compact_conf(conf: str) -> str:
    """'3s2.3p4.(3P).4s' -> '3s²3p⁴4s'; the parent-term parentheses identify core coupling,
    which is not what the reader needs at a glance."""
    conf = re.sub(r"\([^)]*\)", "", conf)
    parts = [p for p in conf.split(".") if p]
    return "".join(re.sub(r"(\d+)$", lambda m: m.group(1).translate(SUP), p) for p in parts)


def parse_tau_arg(s: str) -> float:
    m = re.match(r"^\s*([\d.]+)\s*(ps|ns|us|µs|ms|s)\s*$", s, re.I)
    if not m:
        raise argparse.ArgumentTypeError(f"cannot read a lifetime from {s!r} (try '5ns', '500ps')")
    return float(m.group(1)) * UNITS[m.group(2).lower().replace("µ", "u")]


def tau_slug(t: float) -> str:
    for scale, unit in ((1e-12, "ps"), (1e-9, "ns"), (1e-6, "us"), (1e-3, "ms")):
        if t < scale * 1000:
            v = t / scale
            return f"{v:g}{unit}"
    return f"{t:g}s"


# --------------------------------------------------------------------------- data

def load(d: Path, tau_max: float):
    rows = list(csv.DictReader((d / "levels_with_lifetimes.csv").open(encoding="utf-8")))
    keep = []
    for r in rows:
        if not r["tau_s"] or not r["energy_eV"]:
            continue
        tau = float(r["tau_s"])
        if tau >= tau_max:
            continue
        keep.append({
            "conf": r["configuration"], "term": r["term"],
            "J": float(r["J"]) if r["J"] else None,
            "E": float(r["energy_eV"]), "tau": tau,
            "parity": r["parity"] or "even",
            "decay": float(r["strongest_decay_eV"]) if r["strongest_decay_eV"] else None,
        })
    meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
    return keep, meta, len(rows)


def multiplets(levels):
    groups = defaultdict(list)
    for lv in levels:
        groups[(lv["conf"], lv["term"], lv["parity"])].append(lv)
    out = []
    for (conf, term, parity), members in groups.items():
        members.sort(key=lambda lv: lv["E"])
        taus = [lv["tau"] for lv in members]
        decays = [lv["decay"] for lv in members if lv["decay"] is not None]
        in_win = any(SUB_NS_LO < t < SUB_NS_HI for t in taus)
        txt = f"{compact_conf(conf)} {term_mathtext(term, [lv['J'] for lv in members])}  ·  {fmt_tau_range(taus)}"
        if decays:
            txt += f"  ↓{max(decays):.1f} eV"
        out.append({"E": sum(lv["E"] for lv in members) / len(members),
                    "parity": parity, "text": txt, "in_win": in_win,
                    "members": members})
    out.sort(key=lambda g: g["E"])
    return out


def find_break(values, lo, hi):
    """Split the axis across a band that holds nothing. Li III's fast levels sit at 91.8 eV and
    then 108-121 eV; on one linear axis the empty 17 eV in between eats half the figure."""
    vs = sorted(v for v in values if lo <= v <= hi)
    if len(vs) < 3:
        return None
    span = hi - lo
    g, a, b = max((vs[i + 1] - vs[i], vs[i], vs[i + 1]) for i in range(len(vs) - 1))
    return (a + 0.10 * g, b - 0.10 * g) if (g > 0.30 * span and g > 1.5) else None


def spread(anchors, lo, hi, gap):
    if not anchors:
        return []
    order = sorted(range(len(anchors)), key=lambda i: anchors[i])
    ys = [anchors[i] for i in order]
    for _ in range(300):
        moved = False
        for i in range(len(ys) - 1):
            d = ys[i + 1] - ys[i]
            if d < gap - 1e-12:
                s = (gap - d) / 2
                ys[i] -= s
                ys[i + 1] += s
                moved = True
        if ys[0] < lo:
            ys[0] = lo
            for i in range(1, len(ys)):
                ys[i] = max(ys[i], ys[i - 1] + gap)
        if ys[-1] > hi:
            ys[-1] = hi
            for i in range(len(ys) - 2, -1, -1):
                ys[i] = min(ys[i], ys[i + 1] - gap)
        if not moved:
            break
    out = [0.0] * len(anchors)
    for slot, i in enumerate(order):
        out[i] = ys[slot]
    return out


# --------------------------------------------------------------------------- figures

def empty_card(spectrum: str, d: Path, meta: dict, n_all: int, tau_max: float, out: Path):
    """A species with nothing under the threshold still gets a figure -- and the reason matters:
    'no fast level' and 'no data at all' are completely different answers."""
    fig = plt.figure(figsize=(9.6, 2.5))
    fig.patch.set_facecolor(SURFACE)
    if meta.get("n_lines_with_Aki", 0) == 0 and meta.get("n_tau_from_lidb", 0) == 0:
        why = (f"NIST ASD lists {meta.get('n_lines', 0)} lines for this spectrum but not one "
               "transition probability, and ExoMol LiDB holds no entry above charge 1.\n"
               "No lifetime can be built from either database — untested, not negative.")
    else:
        why = (f"{meta.get('n_levels_with_tau', 0)} of {n_all} levels carry a lifetime, "
               f"but none is shorter than {fmt_tau(tau_max)}.\n"
               "The ladder here is slower than the burst needs.")
    fig.text(0.045, 0.80, f"{spectrum}   —   no level under {fmt_tau(tau_max)}",
             ha="left", va="top", fontsize=15.5, color=INK)
    fig.text(0.045, 0.50, why, ha="left", va="top", fontsize=10, color=INK_2, linespacing=1.6)
    fig.savefig(out.with_suffix(".pdf"), facecolor=SURFACE)
    fig.savefig(out.with_suffix(".png"), dpi=300, facecolor=SURFACE)
    plt.close(fig)


def make_figure(spectrum: str, d: Path, tau_max: float) -> str:
    levels, meta, n_all = load(d, tau_max)
    out = d / f"levels_under_{tau_slug(tau_max)}"
    if not levels:
        empty_card(spectrum, d, meta, n_all, tau_max, out)
        return f"none under {fmt_tau(tau_max)}"

    groups = multiplets(levels)
    es = [lv["E"] for lv in levels]
    lo, hi = min(es), max(es)
    span = max(hi - lo, 0.6)
    ylo, yhi = lo - 0.14 * span, hi + 0.14 * span
    ion = meta.get("ionization_eV")
    show_ion = ion is not None and ylo < ion < hi + 0.5 * span
    if show_ion:
        yhi = max(yhi, ion + 0.06 * span)

    brk = find_break(es, ylo, yhi)
    ranges = [(brk[1], yhi), (ylo, brk[0])] if brk else [(ylo, yhi)]

    def side_max(a, b):
        return max(sum(1 for g in groups if g["parity"] == p and a <= g["E"] <= b)
                   for p in ("even", "odd"))

    per_side = side_max(ylo, yhi)
    fig_h = min(16.0, max(5.4, 2.7 + 0.46 * per_side))
    top_m, bot_m = 1.00, 0.78
    plot_h = fig_h - top_m - bot_m
    total = max(1, sum(side_max(a, b) for a, b in ranges))
    weights = [max(0.11, 0.80 * side_max(a, b) / total + 0.20 * (b - a) / (yhi - ylo))
               for a, b in ranges]
    weights = [w / sum(weights) for w in weights]

    fig, axes = plt.subplots(len(ranges), 1, figsize=(13.0, fig_h), sharex=True,
                             gridspec_kw={"height_ratios": weights, "hspace": 0.05})
    axes = list(axes) if len(ranges) > 1 else [axes]
    fig.patch.set_facecolor(SURFACE)

    for ax, (plo, phi), w in zip(axes, ranges, weights):
        ax.set_facecolor(SURFACE)
        ax.set_xlim(-4.6, 4.6)
        ax.set_ylim(plo, phi)
        ax.set_xticks([])
        for s in ("top", "right", "bottom"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_color("#6b6b66")
        ax.tick_params(axis="y", labelsize=9.5, colors=INK_2)
        ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
        pspan = phi - plo

        if show_ion and plo < ion < phi:
            ax.axhspan(ion, phi, color=AUTO, zorder=0)
            ax.axhline(ion, color=LIMIT, lw=1.0, ls=(0, (7, 4)), zorder=1)
            ax.text(4.5, ion - 0.012 * pspan, f"ionization limit  {ion:.4g} eV",
                    ha="right", va="top", fontsize=8.4, color=INK_3, zorder=6)

        for lv in levels:
            if not (plo <= lv["E"] <= phi):
                continue
            sign = -1.0 if lv["parity"] == "even" else 1.0
            win = SUB_NS_LO < lv["tau"] < SUB_NS_HI
            ax.plot([sign * 0.22, sign * 1.12], [lv["E"]] * 2,
                    color=ACCENT if win else NEUTRAL, lw=2.4 if win else 1.7,
                    solid_capstyle="butt", zorder=4)

        gap = 0.20 / max(plot_h * w, 0.6) * pspan
        for parity, sign, ha in (("even", -1.0, "right"), ("odd", 1.0, "left")):
            grp = [g for g in groups if g["parity"] == parity and plo <= g["E"] <= phi]
            ys = spread([g["E"] for g in grp], plo + 0.01 * pspan, phi - 0.01 * pspan, gap)
            for g, y in zip(grp, ys):
                col = ACCENT if g["in_win"] else INK_2
                ax.annotate(g["text"], xy=(sign * 1.12, g["E"]), xytext=(sign * 1.24, y),
                            ha=ha, va="center", fontsize=9.2, color=col, zorder=7,
                            arrowprops=dict(arrowstyle="-", color=col, lw=0.55, alpha=0.55,
                                            shrinkA=0, shrinkB=1.5))

    if brk:
        kw = dict(marker=[(-1, -0.55), (1, 0.55)], markersize=8, linestyle="none",
                  color="#6b6b66", mec="#6b6b66", mew=1.0, clip_on=False)
        axes[0].plot([0], [0], transform=axes[0].transAxes, **kw)
        axes[1].plot([0], [1], transform=axes[1].transAxes, **kw)

    for x, name in ((-0.67, "even parity"), (0.67, "odd parity")):
        axes[0].text(x, 1.006, name, transform=axes[0].get_xaxis_transform(),
                     ha="center", va="bottom", fontsize=10, color=INK_2, style="italic")

    n_win = sum(1 for lv in levels if SUB_NS_LO < lv["tau"] < SUB_NS_HI)
    fig.text(0.032, 1 - 0.34 / fig_h, f"{spectrum}   levels faster than {fmt_tau(tau_max)}",
             ha="left", va="top", fontsize=16, color=INK)
    fig.text(0.032, 1 - 0.62 / fig_h,
             f"{len(levels)} of {n_all} levels, in {len(groups)} fine-structure multiplets"
             f"   ·   {n_win} inside the 1 – 1000 ps window"
             f"   ·   energies above the ground state at 0 eV",
             ha="left", va="top", fontsize=9.4, color=INK_2)

    fig.legend(handles=[Line2D([], [], color=ACCENT, lw=2.4, label="1 – 1000 ps  (target window)"),
                        Line2D([], [], color=NEUTRAL, lw=1.7,
                               label=f"outside it, still under {fmt_tau(tau_max)}")],
               loc="lower left", bbox_to_anchor=(0.032, 0.30 / fig_h), ncol=2, frameon=False,
               fontsize=8.8, handlelength=2.4, columnspacing=2.4, labelcolor=INK_2)
    fig.text(0.968, 0.34 / fig_h,
             "τ = 1/ΣAₖᵢ over every known channel   ·   ↓ = strongest decay channel",
             ha="right", va="bottom", fontsize=8.2, color=INK_3)

    fig.subplots_adjust(left=0.075, right=0.925, top=1 - top_m / fig_h, bottom=bot_m / fig_h)
    fig.savefig(out.with_suffix(".pdf"), facecolor=SURFACE)
    fig.savefig(out.with_suffix(".png"), dpi=300, facecolor=SURFACE)
    plt.close(fig)
    return f"{len(levels)} levels / {len(groups)} multiplets, {n_win} in window"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species", nargs="*", help='e.g. "Li III" (default: every scraped folder)')
    ap.add_argument("--tau-max", type=parse_tau_arg, default=5e-9,
                    help="lifetime threshold, e.g. 5ns (default) or 500ps")
    args = ap.parse_args()

    targets = ([(s, OUT_ROOT / s.replace(" ", "_")) for s in args.species]
               or [(d.name.replace("_", " "), d) for d in sorted(OUT_ROOT.iterdir())
                   if d.is_dir() and (d / "levels_with_lifetimes.csv").exists()])
    rc = 0
    for spectrum, d in targets:
        try:
            print(f"  {spectrum:<8} {make_figure(spectrum, d, args.tau_max)}")
        except Exception as exc:
            print(f"  {spectrum:<8} FAILED: {exc}")
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
