#!/usr/bin/env python3
"""The candidate shortlist: every scraped level whose lifetime falls in 1 - 1000 ps.

One figure for all species, because most species have nothing in the window and a per-folder
diagram would be a blank page. The two window edges are the physics, not a styling choice:
1 ps is a sanity floor, 1 ns the ceiling above which an ordinary nanosecond laser would
already do the job.

Species identity rides on row grouping and text, not on hue -- five categorical colours in a
scatter cannot be told apart under colour-vision deficiency, and there is nothing here that
deserves emphasis yet. Colour therefore does one job: mark the data.

Usage:
    python3 scripts/subns_map.py
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diagram import SUB_NS_LO, SUB_NS_HI, SUP, fmt_tau, term_mathtext   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8a85"
ACCENT = "#2a78d6"
GRID = "#e6e6e2"
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]


def compact_conf(conf: str) -> str:
    """'3s2.3p4.(3P).4s' -> '3s²3p⁴4s'. Parent-term parentheses are dropped: they identify the
    core coupling, which is not what the reader needs at a glance."""
    conf = re.sub(r"\([^)]*\)", "", conf)
    parts = [p for p in conf.split(".") if p]
    return "".join(re.sub(r"(\d+)$", lambda m: m.group(1).translate(SUP), p) for p in parts)


def collect():
    hits, empty = defaultdict(list), []
    for d in sorted(OUT_ROOT.iterdir()):
        f = d / "levels_with_lifetimes.csv"
        if not d.is_dir() or not f.exists():
            continue
        spectrum = d.name.replace("_", " ")
        groups = defaultdict(list)
        for r in csv.DictReader(f.open(encoding="utf-8")):
            if not r["tau_s"]:
                continue
            tau = float(r["tau_s"])
            if SUB_NS_LO < tau < SUB_NS_HI:
                groups[(r["configuration"], r["term"])].append(r)
        if not groups:
            empty.append(spectrum)
            continue
        for (conf, term), members in groups.items():
            taus = [float(m["tau_s"]) for m in members]
            decays = [float(m["strongest_decay_eV"]) for m in members if m["strongest_decay_eV"]]
            hits[spectrum].append({
                "conf": conf, "term": term,
                "J": [float(m["J"]) for m in members if m["J"]],
                "E": sum(float(m["energy_eV"]) for m in members) / len(members),
                "tau_lo": min(taus), "tau_hi": max(taus),
                "decay": max(decays) if decays else None,
            })
    return hits, empty


def main() -> int:
    hits, empty = collect()
    if not hits:
        print("nothing in the window")
        return 1

    # Species with the most usable levels first -- that ordering is the screening result.
    order = sorted(hits, key=lambda s: (-len(hits[s]), s))
    rows = []
    for sp in order:
        rows.append({"kind": "head", "text": sp, "n": len(hits[sp])})
        for h in sorted(hits[sp], key=lambda h: -h["E"]):
            rows.append({"kind": "item", **h})

    n = len(rows)
    fig_h = 2.9 + 0.265 * n
    fig, ax = plt.subplots(figsize=(12.4, fig_h))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.set_xscale("log")
    ax.set_xlim(0.8, 1250.0)
    ax.set_ylim(n - 0.4, -1.4)
    ax.set_yticks([])

    ticks = [1, 3, 10, 30, 100, 300, 1000]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["1 ps", "3", "10", "30", "100", "300", "1 ns"])
    ax.tick_params(axis="x", length=0, labelsize=9.5, colors=INK_2, pad=7)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for s in ax.spines.values():
        s.set_visible(False)
    for t in ticks:
        ax.axvline(t, color=GRID, lw=0.8, zorder=0)

    # window edges
    for x, txt, ha in ((1, "1 ps", "left"),
                       (1000, "1 ns  ·  a ns laser would already reach it", "right")):
        ax.axvline(x, color="#b9b9b3", lw=1.2, zorder=1)
        ax.text(x, -1.28, txt, ha=ha, va="bottom", fontsize=8.4, color=INK_3)

    # Values live in right-hand columns rather than as labels on the marks: a number on every
    # point is noise, and here it also collided with the row below.
    x_lab, x_tau, x_e, x_d = -0.335, 1.135, 1.265, 1.415
    for i, r in enumerate(rows):
        if r["kind"] == "head":
            ax.text(x_lab, i, r["text"], transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=11, color=INK, weight="semibold")
            ax.text(x_lab + 0.092, i, f"{r['n']} multiplet" + ("s" if r["n"] > 1 else ""),
                    transform=ax.get_yaxis_transform(), ha="left", va="center",
                    fontsize=8.6, color=INK_3)
            continue

        ax.plot([0.8, 1250.0], [i, i], color=GRID, lw=0.6, ls=(0, (1, 3)), zorder=0)
        ax.text(x_lab + 0.01, i, f"{compact_conf(r['conf'])} {term_mathtext(r['term'], r['J'])}",
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=9.4, color=INK_2)

        lo, hi = r["tau_lo"] * 1e12, r["tau_hi"] * 1e12
        if hi / lo > 1.02:                       # fine-structure spread, drawn as a short bar
            ax.plot([lo, hi], [i, i], color=ACCENT, lw=3.2, solid_capstyle="round", zorder=3)
            txt = f"{lo:.0f}–{fmt_tau(r['tau_hi'])}"
        else:
            txt = fmt_tau(r["tau_lo"])
        ax.plot([(lo * hi) ** 0.5], [i], "o", ms=8.0, color=ACCENT,
                mec=SURFACE, mew=1.6, zorder=4)

        for x, val, col in ((x_tau, txt, INK), (x_e, f"{r['E']:.1f}", INK),
                            (x_d, f"{r['decay']:.1f}" if r["decay"] else "—", INK_2)):
            ax.text(x, i, val, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=9.4, color=col)

    for x, head in ((x_tau, "lifetime"), (x_e, "level\neV"), (x_d, "decay\nphoton eV")):
        ax.text(x, -1.12, head, transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=8.4, color=INK_3, linespacing=1.35)

    total = sum(len(v) for v in hits.values())
    fig.text(0.028, 1 - 0.40 / fig_h, "Levels with an intrinsic lifetime inside 1 – 1000 ps",
             ha="left", va="top", fontsize=16.5, color=INK)
    fig.text(0.028, 1 - 0.72 / fig_h,
             f"{total} fine-structure multiplets, in {len(hits)} of the "
             f"{len(hits) + len(empty)} species scraped from NIST ASD and ExoMol LiDB.",
             ha="left", va="top", fontsize=9.6, color=INK_2)
    fig.text(0.028, 1 - 0.94 / fig_h,
             "Lifetime is 1/ΣAₖᵢ summed over every known decay channel, J-resolved; "
             "the decay column is the strongest single channel.",
             ha="left", va="top", fontsize=9.0, color=INK_3)

    fig.text(0.028, 0.46 / fig_h, "Nothing in the window — " + ", ".join(empty) + ".",
             ha="left", va="bottom", fontsize=8.6, color=INK_3)
    fig.text(0.028, 0.20 / fig_h,
             "Xe III, Kr III, Xe IX and Kr IX carry no A-values in ASD at all and have no LiDB "
             "entry above charge 1 — untested rather than negative.",
             ha="left", va="bottom", fontsize=8.6, color=INK_3)

    fig.subplots_adjust(left=0.285, right=0.715, top=1 - 1.38 / fig_h, bottom=0.95 / fig_h)

    out = OUT_ROOT / "subns_candidates"
    fig.savefig(out.with_suffix(".pdf"), facecolor=SURFACE)
    fig.savefig(out.with_suffix(".png"), dpi=300, facecolor=SURFACE)
    plt.close(fig)
    print(f"  {total} multiplets, {len(hits)} species -> {out.with_suffix('.pdf').relative_to(ROOT)} (+ .png)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
