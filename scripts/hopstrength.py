#!/usr/bin/env python3
"""Two- and three-photon hop strength from the one-photon f-values already on disk.

ladder.py decides whether a hop is ALLOWED (energy, parity, J).  This decides whether it is
STRONG, which is a different question and the one that sets the required intensity.

The m-photon amplitude is a sum over real intermediate states,

    M(i->k) = sum_n <k|d|n><n|d|i> / (E_n - E_i - hv)          (m = 2)

and the one-photon f-values in transitions.csv give |<k|d|n>|^2 through the line strength

    S_ik = (3/2) g_i f_ik / dE[Ha]        (atomic units, summed over both levels' sublevels)

so every factor in M is available except the SIGNS of the matrix elements, which f cannot
carry.  Two bounds are therefore reported instead of one number:

    dominant   the single largest |term| -- a lower bound on |M|, and the term that would
               survive if one intermediate dominates (the usual case near resonance)
    aligned    sum of |term| -- the upper bound, reached only if every phase agreed

Both are quoted relative to the strongest hop in the comparison set: this is a ranking tool,
not a cross-section calculator.  Coverage is reported per hop, because a hop whose
intermediates simply have no published A-values scores zero for lack of data, not for lack of
strength.

Usage:
    python3 scripts/hopstrength.py "Cd I" "Hg II" "Xe I" ...
    python3 scripts/hopstrength.py --from-materials 12
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"
ANA = ROOT / "analysis"

HA = 27.211386245981     # eV per hartree
HC = 1239.841984         # eV nm
E_MATCH = 2e-4           # eV, level identification tolerance between the two ASD tables
DELTA_FLOOR = 0.005      # eV, below this a "virtual" intermediate is really a real step


def load_levels(d: Path):
    """(energy, J, g, is_autoionizing).  ASD lists no oscillator strengths INTO the
    autoionizing region, so a hop that lands there can never be scored from this data -- that
    number has to come from a photoabsorption measurement, not from a line list."""
    out = []
    for r in csv.DictReader((d / "levels_with_lifetimes.csv").open(encoding="utf-8")):
        if not r["energy_eV"]:
            continue
        try:
            out.append((float(r["energy_eV"]),
                        float(r["J"]) if r["J"] else None,
                        float(r["g"]) if r["g"] else None,
                        r["above_ionization"] == "yes"))
        except ValueError:
            continue
    out.sort(key=lambda x: x[0])
    return out


def strength_table(d: Path, levels):
    """(i, k) -> line strength in atomic units, i < k.  Levels are identified by energy and J,
    which come from the same ASD tables on both sides, so the match is exact to round-off."""
    import bisect
    E = [x[0] for x in levels]

    def find(e, j):
        k = bisect.bisect_left(E, e - E_MATCH)
        best, bd = None, E_MATCH
        while k < len(E) and E[k] <= e + E_MATCH:
            if (j is None or levels[k][1] is None or abs(levels[k][1] - j) < 1e-6) \
                    and abs(E[k] - e) <= bd:
                best, bd = k, abs(E[k] - e)
            k += 1
        return best

    tab, miss, tot = {}, 0, 0
    f = d / "transitions.csv"
    if not f.exists():
        return tab, 0, 0
    for r in csv.DictReader(f.open(encoding="utf-8")):
        tot += 1
        i = find(float(r["E_i_eV"]), float(r["J_i"]) if r["J_i"] else None)
        k = find(float(r["E_k_eV"]), float(r["J_k"]) if r["J_k"] else None)
        if i is None or k is None or i == k:
            miss += 1
            continue
        dE = abs(levels[k][0] - levels[i][0])
        if dE <= 0:
            continue
        s = 1.5 * float(r["g_i"]) * float(r["f_ik"]) / (dE / HA)
        tab[(min(i, k), max(i, k))] = max(s, tab.get((min(i, k), max(i, k)), 0.0))
    return tab, miss, tot


def partners(tab):
    """level -> the levels it has a published line strength with. The sums below run over these
    only; iterating the whole level list instead makes the three-photon case quadratic in a
    spectrum's size for no gain, since a missing strength contributes nothing."""
    adj = {}
    for a, b in tab:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    return adj


def hop_strength(levels, tab, adj, i, k, m, hv, memo):
    """Dominant and phase-aligned estimates of the m-photon amplitude, in atomic units.

    Memoized on the hop, not the ladder: neighbouring ladders share most of their rungs."""
    ck = (i, k, m, round(hv, 9))
    if ck in memo:
        return memo[ck]
    out = _hop_strength(levels, tab, adj, i, k, m, hv)
    memo[ck] = out
    return out


def _hop_strength(levels, tab, adj, i, k, m, hv):
    def S(a, b):
        return tab.get((min(a, b), max(a, b)))

    Ei = levels[i][0]
    if m == 1:
        s = S(i, k)
        return (s ** 0.5 if s else None), (s ** 0.5 if s else None), 1 if s else 0, None
    if m == 2:
        terms = []
        for n in adj.get(i, ()):
            if n in (i, k):
                continue
            a, b = S(i, n), S(n, k)
            if not a or not b:
                continue
            dn = levels[n][0] - Ei - hv
            terms.append(((a * b) ** 0.5 / max(abs(dn), DELTA_FLOOR), n, dn))
        if not terms:
            return None, None, 0, None
        terms.sort(reverse=True)
        return terms[0][0], sum(t[0] for t in terms), len(terms), (terms[0][1], terms[0][2])
    if m == 3:
        mid = [(n, levels[n][0] - Ei - hv) for n in adj.get(i, ()) if n not in (i, k)]
        terms = []
        for n, dn in mid:
            a = S(i, n)
            for p in adj.get(n, ()):
                if p in (i, k, n):
                    continue
                b, c = S(n, p), S(p, k)
                if not b or not c:
                    continue
                dp = levels[p][0] - Ei - 2 * hv
                terms.append(((a * b * c) ** 0.5
                              / (max(abs(dn), DELTA_FLOOR) * max(abs(dp), DELTA_FLOOR)),
                              n, dn))
        if not terms:
            return None, None, 0, None
        terms.sort(reverse=True)
        return terms[0][0], sum(t[0] for t in terms), len(terms), (terms[0][1], terms[0][2])
    return None, None, 0, None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species", nargs="*", help='e.g. "Cd I" "Hg II"')
    ap.add_argument("--from-materials", type=int, default=0,
                    help="take the top N species from analysis/materials.csv instead")
    ap.add_argument("--per-species", type=int, default=400,
                    help="how many of each species' best ladders to evaluate")
    ap.add_argument("--park-unknown", action="store_true",
                    help="pass through to ladder.py, so the ladders scored here are the same "
                         "ones the material table ranks")
    ap.add_argument("--out", default="hopstrength.csv")
    a = ap.parse_args()

    names = list(a.species)
    if a.from_materials:
        names = [r["spectrum"] for r in
                 csv.DictReader((ANA / "materials.csv").open(encoding="utf-8"))][:a.from_materials]
    if not names:
        print("nothing to do: name some species or pass --from-materials N")
        return 1

    rows = []
    for name in names:
        d = OUT_ROOT / name.replace(" ", "_")
        if not (d / "levels_with_lifetimes.csv").exists():
            print(f"  {name}: not scraped")
            continue
        # re-run the scan for this species alone, so the ladders carry level indices
        tmp = f"_hs_{d.name}.csv"
        subprocess.run([sys.executable, str(ROOT / "scripts" / "ladder.py"),
                        "--species", d.name, "--out", tmp]
                       + (["--park-unknown"] if a.park_unknown else []),
                       check=True, stdout=subprocess.DEVNULL)
        lad = list(csv.DictReader((ANA / tmp).open(encoding="utf-8")))
        (ANA / tmp).unlink()
        if not lad:
            continue
        levels = load_levels(d)
        tab, miss, tot = strength_table(d, levels)
        adj, memo = partners(tab), {}
        lad.sort(key=lambda r: ({"measured": 0, "autoionizing": 1, "none": 2}[r["subns_evidence"]],
                                {"A": 0, "B": 1, "C": 2}[r["detuning_tier"]],
                                -int(r["order"]), float(r["worst_detuning_meV"])))
        used = 0
        for r in lad[:a.per_species]:
            idx = [int(r["start_idx"])] + [int(x) for x in r["path_idx"].split()]
            ms = [int(x) for x in r["photons_per_hop"].split("+")]
            cs = [float(x) for x in r["colours_nm"].split("+")]
            dom, aln, cov, kind, weakest = [], [], [], [], None
            for h, (m, nm) in enumerate(zip(ms, cs)):
                ai = levels[idx[h + 1]][3]
                kind.append("AI" if ai else "bound")
                if ai:
                    dom.append(None)
                    aln.append(None)
                    cov.append(0)
                    continue
                dd, ss, n, best = hop_strength(levels, tab, adj, idx[h], idx[h + 1],
                                               m, HC / nm, memo)
                dom.append(dd)
                aln.append(ss)
                cov.append(n)
                if dd is not None and (weakest is None or dd < weakest):
                    weakest = dd
            bound = [x for x, k in zip(dom, kind) if k == "bound"]
            covered = bool(bound) and all(x is not None for x in bound)
            rows.append({
                "spectrum": name, "hops": r["hops"], "order": r["order"],
                "photons_per_hop": r["photons_per_hop"], "colours_nm": r["colours_nm"],
                "detuning_tier": r["detuning_tier"],
                "worst_detuning_meV": r["worst_detuning_meV"],
                "subns_evidence": r["subns_evidence"],
                "start": r["start"], "start_eV": r["start_eV"],
                "path": r["path"], "path_eV": r["path_eV"],
                "hop_dominant_au": " | ".join("-" if x is None else f"{x:.3g}" for x in dom),
                "hop_aligned_au": " | ".join("-" if x is None else f"{x:.3g}" for x in aln),
                "intermediates_used": " | ".join(str(x) for x in cov),
                "hop_lands_on": " | ".join(kind),
                "weakest_bound_hop_au": f"{weakest:.4g}" if covered and weakest else "",
                "bound_hops_scored": "yes" if covered else "no",
                "final_hop": ("into the autoionizing region: not scorable from a line list"
                              if kind[-1] == "AI" else
                              ("scored" if dom[-1] is not None else "no f data")),
            })
            used += covered
        print(f"  {name:9s} {len(tab):6d} strengths, {miss} of {tot} lines unmatched, "
              f"{used}/{min(len(lad), a.per_species)} ladders with every BOUND hop scored")

    scored = [r for r in rows if r["weakest_bound_hop_au"]]
    top = max((float(r["weakest_bound_hop_au"]) for r in scored), default=0.0)
    for r in rows:
        r["relative_strength"] = (f"{float(r['weakest_bound_hop_au']) / top:.3g}"
                                  if r["weakest_bound_hop_au"] and top else "")
    if scored:
        rows.sort(key=lambda r: -float(r["relative_strength"] or 0))
    ANA.mkdir(exist_ok=True)
    out = ANA / a.out
    cols = list(rows[0].keys()) if rows else []
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"  {len(rows)} ladders ({len(scored)} with every bound hop scored) -> "
          f"{out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
