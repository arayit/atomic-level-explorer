#!/usr/bin/env python3
"""Hop-grid scan: which cached species admit a real, sequential burst ladder?

The scheme (burst-ladder-project-handoff.md, s.2-3): split one impossibly high-order
multiphoton process into low-order hops through REAL intermediate states, one hop per pulse
of a GHz burst.  A ladder is emitted only if all of the following hold.

  1. Each hop is m photons of ONE available driver colour, m in {1,2,3} (floor set by --min-m).
  2. Hop energy matches a real level pair:  |dE - m*hv| <= --tol.
  3. Laporte: m photons flip parity m times, so parity_k == parity_i xor (m odd).
  4. |J_k - J_i| <= m, and J = 0 -> J = 0 is closed for odd m.
  5. Start state holds population for a whole burst: the ground state, or a level below the
     first allowed decay that is either >= 1 us or has no radiative channel at all.
  6. Every state the ladder PARKS population in -- all but the last -- has a measured
     lifetime longer than the intra-burst spacing, and sits below the ionization limit.
  7. >= 2 hops and total synthetic order sum(m) >= --min-order.

Criterion 8 of the handoff, "at least one rung intrinsically sub-ns", is NOT used as a filter
but recorded as an evidence class, because ASD cannot express it for half the candidates:

  measured      a populated state carries an ASD lifetime inside (spacing, 1 ns).
  autoionizing  the ladder lands above the ionization limit, so the state decays by
                autoionization -- intrinsically fs-ps, but its width lives in linewidth
                tables, not in ASD.  Needs one literature number to become quantitative.
  none          every populated state is ns or slower: a working ladder, but a ns laser
                could drive it too, so it carries no novelty.

tau from ASD is 1/sum(A_ki) over the channels NIST lists, so it is an UPPER BOUND.  That makes
"tau < 1 ns" safe (the truth is shorter) and "tau > spacing" only as good as the line list --
hence the completeness flag carried through to the output.

Usage:
    python3 scripts/ladder.py                   # 6 burst colours, m in {2,3}, +-150 meV
    python3 scripts/ladder.py --sh              # also allow each colour's second harmonic
    python3 scripts/ladder.py --species Ca_I --min-m 1
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"
ANA = ROOT / "analysis"
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]

HC = 1239.841984                       # eV nm
DRIVERS_NM = [1030, 1060, 1550, 1600, 1900, 2000]

SPACING_S = 20e-12                     # intra-burst spacing, 50 GHz (handoff s.1)
SUBNS_HI = 1e-9                        # top of the novelty window
METASTABLE_S = 1e-6                    # a start state must hold population for a whole burst
FAST_S = 1e-7                          # a level this fast marks the first allowed decay

DET_A, DET_B = 0.025, 0.060            # detuning quality bands, eV
MAX_HOPS = 4
BRANCH_PARK = 30                       # expansion hops kept per node, best detuning first
BRANCH_FINAL = 60                      # last-hop candidates kept per node
PATH_CAP = 250_000                     # per-species safety stop

GAS_AT_RT = {"H", "He", "N", "O", "F", "Ne", "Cl", "Ar", "Kr", "Xe"}
VAPOUR_1PA_K = {                       # CRC 84th ed.; metals: Alcock, Itkin & Horrigan (1984)
    "H": 10, "He": 1.3, "Li": 797, "Be": 1462, "B": 2348, "N": 37,
    "F": 38, "Ne": 12, "Na": 554, "Mg": 701, "Al": 1482, "Si": 1908, "P": 279, "S": 375,
    "Cl": 128, "Ar": 47, "K": 473, "Ca": 864, "Ti": 1982, "Fe": 1728, "Kr": 59, "Rb": 434,
    "Sr": 796, "Cd": 530, "Xe": 83, "Cs": 418, "Ba": 911, "Hg": 315,
    "V": 2101, "Cr": 1656, "Mn": 1228, "Co": 1790, "Ni": 1783, "Cu": 1509, "Zn": 610,
    "Ga": 1310, "Ge": 1644, "In": 1196, "Sn": 1497, "Tl": 882, "Pb": 978,
}
ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19,
    "Ca": 20, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Cd": 48,
    "In": 49, "Sn": 50, "Xe": 54, "Cs": 55, "Ba": 56, "Hg": 80, "Tl": 81, "Pb": 82,
}


def prep_cost(el: str, q: int) -> tuple[int, str]:
    """Two independent costs: getting the ELEMENT into the gas phase, and stripping q
    electrons off it.  Kept additive and coarse -- this ranks, it does not model."""
    if el in GAS_AT_RT:
        vap, note = 0, "gas at 20 C"
    else:
        t = VAPOUR_1PA_K.get(el)
        if t is None:
            vap, note = 4, "no vapour datum"
        elif t <= 500:
            vap, note = 1, f"1 Pa at {t - 273:.0f} C"
        elif t <= 900:
            vap, note = 2, f"1 Pa at {t - 273:.0f} C"
        else:
            vap, note = 3, f"1 Pa at {t - 273:.0f} C"
    if q == 1:
        return vap + 2, note + ", discharge"
    if q >= 2:
        return vap + 2 + 2 * q, note + f", {q}x ionized"
    return vap, note


class Level:
    __slots__ = ("i", "E", "J", "p", "tau", "floor", "conf", "term", "above", "complete", "decay")

    def __init__(self, i, E, J, p, tau, floor, conf, term, above, complete, decay):
        self.i, self.E, self.J, self.p = i, E, J, p
        self.tau, self.floor, self.conf, self.term = tau, floor, conf, term
        self.above, self.complete, self.decay = above, complete, decay

    def name(self):
        return f"{self.conf} {self.term}" + (f" J={self.J:g}" if self.J is not None else "")


def tau_floor(tau, channels: str):
    """ASD gives tau = 1/sum(A) over the channels it lists, so tau is an UPPER bound.  Assume no
    missing channel is stronger than the strongest listed one: with k of n channels known,
    sum(A) <= (n/k) sum(A_known), i.e. tau >= tau_listed * k/n.  Conservative, and enough to
    decide whether a level really outlives the pulse spacing."""
    if tau is None:
        return None
    m = re.match(r"\s*(\d+)\s*/\s*(\d+)\s*$", channels or "")
    if not m:
        return tau
    k, n = int(m.group(1)), int(m.group(2))
    return tau * k / n if n else tau


def load(d: Path) -> list[Level]:
    rows = []
    for r in csv.DictReader((d / "levels_with_lifetimes.csv").open(encoding="utf-8")):
        if not r["energy_eV"]:
            continue
        try:
            E = float(r["energy_eV"])
        except ValueError:
            continue
        tau = float(r["tau_s"]) if r["tau_s"] else None
        rows.append(Level(
            0, E,
            float(r["J"]) if r["J"] else None,
            1 if r["parity"] == "odd" else (0 if r["parity"] == "even" else None),
            tau, tau_floor(tau, r.get("asd_channels", "")),
            r["configuration"], r["term"],
            r["above_ionization"] == "yes",
            r["asd_tau_complete"] == "yes",
            float(r["max_decay_eV"]) if r["max_decay_eV"] else None,
        ))
    rows.sort(key=lambda x: x.E)
    for k, lv in enumerate(rows):
        lv.i = k
    return rows


def allowed(a: Level, b: Level, m: int) -> bool:
    if a.p is None or b.p is None:
        return False                                   # unknown parity is not a licence
    if (a.p ^ b.p) != (m & 1):
        return False
    if a.J is not None and b.J is not None:
        if abs(b.J - a.J) > m:
            return False
        if (m & 1) and a.J == 0 and b.J == 0:
            return False
    return True


def virtual_gap(lv, E, a: Level, m: int, hv: float) -> float:
    """How far the m-photon hop has to reach off resonance on its way up.

    An m-photon hop climbs through m-1 virtual points at E_i + s*hv.  Its rate carries a
    factor 1/prod(Delta_s), so the bottleneck is the LARGEST distance from a virtual point to
    a real level of the parity that s photons produce.  Returned in eV: small means the hop is
    stepwise-enhanced, large means it is genuinely virtual and therefore weak."""
    worst = 0.0
    for step in range(1, m):
        v = a.E + step * hv
        want = a.p ^ (step & 1)
        k = bisect.bisect_left(E, v)
        d = float("inf")
        for j in range(max(0, k - 12), min(len(lv), k + 12)):
            if lv[j].p == want:
                d = min(d, abs(lv[j].E - v))
        worst = max(worst, d)
    return worst


def start_states(lv, ion_eV):
    if not lv:
        return []
    out = [(lv[0], "ground")]
    fast = [x.E for x in lv if x.tau is not None and x.tau < FAST_S]
    ceiling = min(fast) if fast else (ion_eV if ion_eV is not None else 0.0)
    for x in lv[1:]:
        if x.E >= ceiling or x.above or (ion_eV is not None and x.E >= ion_eV):
            continue
        if x.tau is None or x.tau >= METASTABLE_S:
            out.append((x, "metastable"))
        if len(out) > 24:
            break
    return out


def build_hops(lv, grid, tol, parkable, measured):
    """node -> (expansion hops, final-hop candidates).  Expansion may only land on a level that
    can hold population; the last hop may land anywhere.

    Expansion candidates are ordered measured-lifetime-first, then by detuning.  Without that,
    switching --park-unknown on would push levels with a published lifetime out of the branch
    budget and the relaxed scan would return FEWER ladders than the strict one -- a pruning
    artefact rather than a result.  This ordering makes relaxed a true superset of strict."""
    E = [x.E for x in lv]
    n = len(lv)
    exp_, fin_ = {}, {}
    for a in lv:
        acc = []
        for nm, hv, m in grid:
            h = m * hv
            lo = bisect.bisect_left(E, a.E + h - tol)
            hi = bisect.bisect_right(E, a.E + h + tol)
            for k in range(lo, min(hi, n)):
                b = lv[k]
                if b.i == a.i or not allowed(a, b, m):
                    continue
                acc.append((b.i, m, nm, b.E - a.E - h))
        acc.sort(key=lambda t: abs(t[3]))
        fin_[a.i] = acc[:BRANCH_FINAL]
        ok = [t for t in acc if parkable[t[0]]]
        ok.sort(key=lambda t: (0 if measured[t[0]] else 1, abs(t[3])))
        exp_[a.i] = ok[:BRANCH_PARK]
    return exp_, fin_


def scan(lv, exp_, fin_, starts, max_hops, min_order):
    found, guard, cut = [], 0, False
    for s, kind in starts:
        stack = [(s.i, (), 0)]
        while stack:
            node, path, order = stack.pop()
            guard += 1
            if guard > PATH_CAP:
                cut = True
                break
            seen = {s.i} | {t for t, _, _, _ in path}
            for tgt, m, nm, delta in fin_.get(node, ()):
                if tgt in seen:
                    continue
                if len(path) + 1 >= 2 and order + m >= min_order:
                    found.append((s, kind, path + ((tgt, m, nm, delta),), order + m))
            if len(path) + 1 < max_hops:
                for tgt, m, nm, delta in exp_.get(node, ()):
                    if tgt not in seen:
                        stack.append((tgt, path + ((tgt, m, nm, delta),), order + m))
    return found, cut


def dedupe(found, lv):
    """The fine-structure components of one multiplet are the same physical ladder; keep the
    best-detuned representative of each (start term, landing terms, m-pattern, colour-pattern).

    Whether each rung's lifetime is published is part of the key. Without it, a J component
    ASD has no lifetime for can displace a measured sibling on a slightly better detuning, and
    the ladder silently stops being verifiable -- so --park-unknown would DROP verified ladders
    that the strict scan finds."""
    best = {}
    for s, kind, path, order in found:
        key = (s.conf, s.term,
               tuple((lv[t].conf, lv[t].term) for t, _, _, _ in path),
               tuple(m for _, m, _, _ in path),
               tuple(nm for _, _, nm, _ in path),
               tuple(lv[t].tau is None for t, _, _, _ in path[:-1]))
        worst = max(abs(d) for _, _, _, d in path)
        if key not in best or worst < best[key][0]:
            best[key] = (worst, (s, kind, path, order))
    return [v[1] for v in best.values()]


COLS = ["spectrum", "element", "Z", "charge", "detuning_tier", "subns_evidence",
        "prep_cost", "prep_note", "start", "start_kind", "start_eV", "start_idx", "path_idx", "hops", "order",
        "worst_detuning_meV", "colours_nm", "photons_per_hop", "path", "path_eV",
        "path_tau_s", "min_parked_tau_s", "min_parked_tau_floor_s", "final_eV", "final_eV_neutral_frame",
        "final_above_ionization", "final_max_decay_eV", "tau_complete", "fast_rung_at", "fast_rung_tau_s", "unverified_parked",
        "worst_virtual_gap_eV", "direct_shortcut", "species_truncated"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=0.150, help="hop energy tolerance, eV")
    ap.add_argument("--min-m", type=int, default=2, help="lowest photon order allowed per hop")
    ap.add_argument("--max-hops", type=int, default=MAX_HOPS)
    ap.add_argument("--min-order", type=int, default=4, help="lowest total synthetic order")
    ap.add_argument("--spacing", type=float, default=SPACING_S, help="intra-burst spacing, s")
    ap.add_argument("--sh", action="store_true", help="also allow each colour's 2nd harmonic")
    ap.add_argument("--keep-plain", action="store_true", help="keep ladders with no fast rung")
    ap.add_argument("--species", default="", help="restrict to one folder name, for checking")
    ap.add_argument("--max-charge", type=int, default=99, help="skip spectra above this charge")
    ap.add_argument("--park-unknown", action="store_true",
                    help="also park on levels ASD gives no lifetime for, and count them")
    ap.add_argument("--out", default="ladders.csv")
    a = ap.parse_args()

    nms = list(DRIVERS_NM) + ([n / 2 for n in DRIVERS_NM] if a.sh else [])
    grid = [(nm, HC / nm, m) for nm in nms for m in range(max(1, a.min_m), 4)]

    ions = {}
    for d in sorted(OUT_ROOT.iterdir()):
        if (d / "meta.json").exists():
            m = json.loads((d / "meta.json").read_text(encoding="utf-8"))
            ions[m["spectrum"]] = m.get("ionization_eV")

    rows, seen_species = [], set()
    for d in sorted(OUT_ROOT.iterdir()):
        f = d / "levels_with_lifetimes.csv"
        if not d.is_dir() or not f.exists():
            continue
        if a.species and d.name != a.species:
            continue
        if ROMAN.index(d.name.split("_")[1]) > a.max_charge:
            continue
        sp = d.name.replace("_", " ")
        el, roman = sp.split()
        q = ROMAN.index(roman)
        lv = load(d)
        if len(lv) < 3:
            continue
        ip = ions.get(sp)
        bound = [not x.above and (ip is None or x.E < ip) for x in lv]
        measured = [x.floor is not None and x.floor > a.spacing and b
                    for x, b in zip(lv, bound)]
        parkable = [m or (a.park_unknown and x.tau is None and b)
                    for x, m, b in zip(lv, measured, bound)]
        starts = start_states(lv, ip)
        exp_, fin_ = build_hops(lv, grid, a.tol, parkable, measured)
        raw, cut = scan(lv, exp_, fin_, starts, a.max_hops, a.min_order)
        got = dedupe(raw, lv)
        Evals = [x.E for x in lv]
        print(f"    {sp:9s} {len(lv):5d} levels, {sum(parkable):4d} parkable, "
              f"{len(starts):2d} starts -> {len(got):6d} ladders"
              + ("  [path cap hit]" if cut else ""), flush=True)
        if not got:
            continue

        off, ok = 0.0, True
        for k in range(q):
            prev = ions.get(f"{el} {ROMAN[k]}")
            if prev is None:
                ok = False
                break
            off += prev
        cost, note = prep_cost(el, q)

        for s, kind, path, order in got:
            states = [lv[t] for t, _, _, _ in path]
            fin = states[-1]
            parked = [s] + states[:-1]
            inter = states[:-1]
            # A rung that PARKS population has to outlive the pulse spacing, so its usable
            # window is (spacing, 1 ns).  The last rung parks nothing -- the faster it decays
            # the better, so the only bound there is 1 ns.
            fast_i = [x for x in inter if x.tau is not None and a.spacing < x.tau < SUBNS_HI]
            fast_f = fin.tau is not None and fin.tau < SUBNS_HI
            if fast_i or fast_f:
                ev = "measured"
            elif fin.above:
                ev = "autoionizing"
            else:
                ev = "none"
            where = ("both" if fast_i and fast_f else
                     ("intermediate" if fast_i else ("final" if fast_f else "")))
            ftau = min([x.tau for x in fast_i] + ([fin.tau] if fast_f else []), default=None)
            if ev == "none" and not a.keep_plain:
                continue
            worst = max(abs(dd) for _, _, _, dd in path)
            used = {nm for _, _, nm, _ in path}
            short = ""
            for nm, hv, m in grid:
                if nm in used and abs(fin.E - s.E - m * hv) <= a.tol and allowed(s, fin, m):
                    short = f"{m}x{nm:g}nm"
            vgap = 0.0
            prev = s
            for tgt, m, nm, _ in path:
                vgap = max(vgap, virtual_gap(lv, Evals, prev, m, HC / nm))
                prev = lv[tgt]
            rows.append({
                "spectrum": sp, "element": el, "Z": ELEMENT_Z.get(el, ""), "charge": q,
                "detuning_tier": "A" if worst <= DET_A else ("B" if worst <= DET_B else "C"),
                "subns_evidence": ev, "prep_cost": cost, "prep_note": note,
                "start": s.name(), "start_kind": kind, "start_eV": round(s.E, 4),
                "start_idx": s.i, "path_idx": " ".join(str(t) for t, _, _, _ in path),
                "hops": len(path), "order": order,
                "worst_detuning_meV": round(worst * 1000, 1),
                "colours_nm": "+".join(f"{nm:g}" for _, _, nm, _ in path),
                "photons_per_hop": "+".join(str(m) for _, m, _, _ in path),
                "path": " -> ".join(x.name() for x in states),
                "path_eV": " -> ".join(f"{x.E:.3f}" for x in states),
                "path_tau_s": " -> ".join("?" if x.tau is None else f"{x.tau:.2e}" for x in states),
                "min_parked_tau_s": min((x.tau for x in parked if x.tau), default=""),
                "min_parked_tau_floor_s": min((x.floor for x in parked if x.floor), default=""),
                "final_eV": round(fin.E, 4),
                "final_eV_neutral_frame": round(off + fin.E, 4) if ok else "",
                "final_above_ionization": "yes" if fin.above else "no",
                "final_max_decay_eV": round(fin.decay, 3) if fin.decay else "",
                "tau_complete": "yes" if all(x.complete for x in states[:-1]) else "no",
                "fast_rung_at": where,
                "fast_rung_tau_s": f"{ftau:.3e}" if ftau else "",
                # counts the RUNGS only: the start is metastable by construction and ASD
                # routinely leaves its lifetime blank, which is not the same as unverified
                "unverified_parked": sum(1 for x in inter if x.tau is None),
                "worst_virtual_gap_eV": round(vgap, 3) if vgap != float("inf") else "",
                "direct_shortcut": short,
                "species_truncated": "yes" if cut else "no",
            })
            seen_species.add(sp)

    rank = {"A": 0, "B": 1, "C": 2}
    ev_rank = {"measured": 0, "autoionizing": 1, "none": 2}
    rows.sort(key=lambda r: (r["prep_cost"], ev_rank[r["subns_evidence"]],
                             rank[r["detuning_tier"]], -r["order"], r["worst_detuning_meV"]))
    ANA.mkdir(exist_ok=True)
    out = ANA / a.out
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)

    print(f"  {len(seen_species)} species carry at least one admissible ladder")
    print(f"  {len(rows)} distinct ladders -> {out.relative_to(ROOT)}")
    for e in ("measured", "autoionizing", "none"):
        n = sum(1 for r in rows if r["subns_evidence"] == e)
        if n:
            print(f"    sub-ns evidence {e:13s}: {n}")
    for t in "ABC":
        print(f"    detuning tier {t}: {sum(1 for r in rows if r['detuning_tier'] == t)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
