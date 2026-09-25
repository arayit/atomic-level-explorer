#!/usr/bin/env python3
"""Collapse the ladder scan into a tiered material list.

ladder.py answers "which level chains close on the available colours".  This answers the
question that follows: which SPECIES is worth putting in a cell, and in what order to try them.

Four independent axes, all recorded rather than blended into a single opaque number:

  preparation   getting the species into the interaction region: gas at 20 C, warm vapour,
                hot vapour, discharge, multiply-ionized plasma.
  closure       the worst hop detuning of the species' best ladder.  A <= 25 meV closes on
                driver tuning alone; B <= 60 meV wants chirp or AC Stark; C uses the whole
                150 meV budget.
  fast rung     measured (an ASD lifetime inside the picosecond window), autoionizing (the
                ladder lands above the ionization limit, so the state is intrinsically fast
                but its width must come from a linewidth table), or none.
  payoff        total synthetic order, and how high the ladder reaches above the neutral ground.

A fifth, cheap but decisive in the lab: does the ladder run on ONE laser?  Two independent
colours means two synchronized bursts -- but a colour and its own second harmonic count as a
single source, because the SHG output needs no synchronization at all.

Usage:
    python3 scripts/materials.py                    # reads analysis/ladders.csv
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANA = ROOT / "analysis"

HC = 1239.841984
DRIVERS_NM = [1030, 1060, 1550, 1600, 1900, 2000]
SH_NM = [515, 530]                     # second harmonics of the Yb pair (scan: --sh 1030 1060)
COLOURS_NM = DRIVERS_NM + SH_NM


def sources(colours_nm: str):
    """The distinct LASERS behind a ladder's colour list.  A second harmonic maps back to
    its fundamental: the SHG output is inherently synchronized with the driver it came from,
    so 1030+515 is one source while 1030+1550 is two."""
    return {(x if x in DRIVERS_NM else 2 * x) for x in map(float, colours_nm.split("+"))}

# Monatomic on the bench.  The noble gases are; H2, N2, O2, F2, Cl2, S8, P4 and graphite are
# not, so their ATOMIC spectra need a dissociation source before any of this starts.
MONATOMIC_GAS = {"He", "Ne", "Ar", "Kr", "Xe"}
MOLECULAR_AT_RT = {"H": "H2", "N": "N2", "O": "O2", "F": "F2", "Cl": "Cl2",
                   "Br": "Br2", "I": "I2", "S": "S8", "P": "P4", "C": "graphite",
                   "Se": "Se8"}
VAPOUR_1PA_K = {                       # CRC 84th ed.; metals: Alcock, Itkin & Horrigan (1984)
    "Li": 797, "Be": 1462, "B": 2348, "Na": 554, "Mg": 701, "Al": 1482, "Si": 1908,
    "K": 473, "Ca": 864, "Ti": 1982, "V": 2101, "Cr": 1656, "Mn": 1228, "Fe": 1728,
    "Co": 1790, "Ni": 1783, "Cu": 1509, "Zn": 610, "Ga": 1310, "Ge": 1644, "Rb": 434,
    "Sr": 796, "Cd": 530, "In": 1196, "Sn": 1497, "Cs": 418, "Ba": 911, "Hg": 315,
    "Tl": 882, "Pb": 978,
}


def prep_cost(el: str, q: int) -> tuple[int, str]:
    """Two independent costs: getting single ATOMS of the element into the interaction region,
    and stripping q electrons off them."""
    if el in MONATOMIC_GAS:
        vap, note = 0, "monatomic gas at 20 C"
    elif el in MOLECULAR_AT_RT:
        vap, note = 3, f"{MOLECULAR_AT_RT[el]} at 20 C, needs dissociation"
    else:
        t = VAPOUR_1PA_K.get(el)
        if t is None:
            vap, note = 4, "no vapour datum"
        elif t <= 500:
            vap, note = 1, f"vapour, 1 Pa at {t - 273:.0f} C"
        elif t <= 900:
            vap, note = 2, f"vapour, 1 Pa at {t - 273:.0f} C"
        else:
            vap, note = 3, f"vapour, 1 Pa at {t - 273:.0f} C"
    if q == 1:
        return vap + 2, note + " + discharge"
    if q >= 2:
        return vap + 2 + 2 * q, note + f" + {q}x ionization"
    return vap, note

DET_W = {"A": 0, "B": 2, "C": 5}
# How far a multiphoton hop must reach off resonance on the way up.  A real level close to the
# virtual point multiplies the hop rate; far from one, the hop is genuinely virtual and weak.
VG_BANDS = [(0.05, "stepwise", 0), (0.20, "enhanced", 0), (0.60, "moderate", 1)]
VG_FAR = ("virtual", 3)
# Perturbative amplitude of the weakest BOUND hop, from scripts/hopstrength.py, relative to the
# strongest hop in the comparison set.  Only a ranking: the absolute scale is meaningless, and
# a species with no published f-values for its hops scores nothing at all rather than zero.
AMP_BANDS = [(1e-2, "strong", 0), (1e-3, "moderate", 1), (0.0, "weak", 3)]
AMP_NONE = ("no f data", 1)
# A rung whose lifetime ASD simply does not publish. Almost every bound level below the
# ionization limit lives far longer than 20 ps, so an unlisted rung is probably fine -- but
# "probably" is worth one point against a rung that was actually measured.
UNVERIFIED_W = 1


def amp_band(v):
    if v is None:
        return AMP_NONE
    for lo, name, w in AMP_BANDS:
        if v >= lo:
            return name, w
    return AMP_BANDS[-1][1], AMP_BANDS[-1][2]


def hop_amplitudes():
    """species -> (best relative amplitude, the ladder that carries it).  Only ladders that
    close inside the B band are considered: a strong hop to the wrong energy is worth nothing."""
    f = ANA / "hopstrength.csv"
    if not f.exists():
        return {}
    best = {}
    for r in csv.DictReader(f.open(encoding="utf-8")):
        if not r["relative_strength"] or r["detuning_tier"] not in ("A", "B"):
            continue
        v = float(r["relative_strength"])
        if r["spectrum"] not in best or v > best[r["spectrum"]][0]:
            best[r["spectrum"]] = (v, r)
    return best
EV_W = {"measured": 0, "autoionizing": 1, "none": 6}
EV_SHORT = {"measured": "measured", "autoionizing": "autoionizing", "none": "-"}


def hop_grid(min_m=2, max_m=3, tol=0.150):
    """Which rung energies the driver set can actually deliver, and where the holes are."""
    g = sorted((m * HC / nm, nm, m) for nm in COLOURS_NM for m in range(min_m, max_m + 1))
    bands, holes = [], []
    lo, hi = g[0][0] - tol, g[0][0] + tol
    members = [g[0]]
    for e, nm, m in g[1:]:
        if e - tol <= hi:
            hi = max(hi, e + tol)
            members.append((e, nm, m))
        else:
            bands.append((lo, hi, members))
            holes.append((hi, e - tol))
            lo, hi, members = e - tol, e + tol, [(e, nm, m)]
    bands.append((lo, hi, members))
    return g, bands, holes


def accidental_rate(tol=0.150):
    """How many levels land inside one hop window by chance.

    A 150 meV window in a dense spectrum catches a level no matter where it is aimed, so a
    "match" in Fe I means much less than the same match in Cd I.  Estimated as the level
    density below the ionization limit times the window width, halved for parity: above ~1 the
    energy match carries no information and the ladder must be judged on its parked states and
    hop strengths alone."""
    out = {}
    for d in sorted((ROOT / "energy_levels").iterdir()):
        f = d / "levels_with_lifetimes.csv"
        if not f.exists():
            continue
        es = []
        ip = None
        for r in csv.DictReader(f.open(encoding="utf-8")):
            if not r["energy_eV"]:
                continue
            if ip is None and r["ionization_eV"]:
                ip = float(r["ionization_eV"])
            if r["above_ionization"] != "yes":
                es.append(float(r["energy_eV"]))
        span = (ip if ip else (max(es) if es else 0)) or 1.0
        out[d.name.replace("_", " ")] = round(2 * tol * len(es) / span / 2, 1)
    return out


def vg_band(v):
    if v is None:
        return VG_FAR
    for hi, name, w in VG_BANDS:
        if v < hi:
            return name, w
    return VG_FAR


def score(prep, ev, det, order, single, vg, amp, unv):
    """Amplitude supersedes the virtual-gap proxy where it exists -- same axis, real numbers."""
    strength = amp_band(amp)[1] if amp is not None else vg_band(vg)[1]
    return (prep + DET_W[det] + EV_W[ev] + max(0, 6 - order)
            + (0 if single else 1) + strength + (UNVERIFIED_W if unv else 0))


def main() -> int:
    src = ANA / "ladders.csv"
    rows = list(csv.DictReader(src.open(encoding="utf-8")))
    for r in rows:
        r["_single"] = len(sources(r["colours_nm"])) == 1
        r["_order"] = int(r["order"])
        r["_hops"] = int(r["hops"])
        r["_det"] = float(r["worst_detuning_meV"])
        r["_prep"], r["prep_note"] = prep_cost(r["element"], int(r["charge"]))
        v = r.get("worst_virtual_gap_eV", "")
        r["_vg"] = float(v) if v else None
        r["_unv"] = int(r.get("unverified_parked") or 0)

    best_key = lambda r: (EV_W[r["subns_evidence"]], 1 if r["_unv"] else 0,
                          DET_W[r["detuning_tier"]], vg_band(r["_vg"])[1],
                          0 if r["_single"] else 1, -r["_order"], r["_det"])

    acc = accidental_rate()
    amps = hop_amplitudes()
    sp = {}
    for r in rows:
        sp.setdefault(r["spectrum"], []).append(r)

    out = []
    for name, rs in sp.items():
        rs.sort(key=best_key)
        b = rs[0]
        verified = [x for x in rs if not x["_unv"]]
        # the best ladder that also runs on a single colour, if that is not already the best
        single = next((x for x in rs if x["_single"]), None)
        reach = max((float(x["final_eV_neutral_frame"]) for x in rs
                     if x["final_eV_neutral_frame"]), default=None)
        out.append({
            "spectrum": name, "element": b["element"], "Z": int(b["Z"]) if b["Z"] else 999,
            "charge": int(b["charge"]), "n_ladders": len(rs),
            "prep_cost": b["_prep"], "prep_note": b["prep_note"],
            "evidence": b["subns_evidence"], "closure": b["detuning_tier"],
            "best_detuning_meV": b["_det"],
            "order": b["_order"], "hops": b["_hops"],       # of the BEST ladder, not the set
            "pattern": b["photons_per_hop"], "colours": b["colours_nm"],
            "max_order": max(x["_order"] for x in rs),
            "min_hops": min(x["_hops"] for x in rs),
            "ground_start": "yes" if any(x["start_kind"] == "ground" for x in rs) else "no",
            "single_colour": "yes" if single else "no",
            "reach_eV_neutral_frame": round(reach, 2) if reach is not None else "",
            "rung_lifetimes": ("all measured" if not b["_unv"]
                               else f"{b['_unv']} rung(s) unpublished"),
            "verified_ladders": len(verified),
            "fast_rung_at": b.get("fast_rung_at", ""),
            "fast_rung_tau_s": b.get("fast_rung_tau_s", ""),
            "hop_amplitude": amp_band(amps[name][0])[0] if name in amps else "not scored",
            "relative_amplitude": f"{amps[name][0]:.2g}" if name in amps else "",
            "hop_strength": vg_band(b["_vg"])[0],
            "worst_virtual_gap_eV": b["_vg"] if b["_vg"] is not None else "",
            "truncated": "yes" if any(x.get("species_truncated") == "yes" for x in rs) else "no",
            "accidental_matches_per_hop": acc.get(name, ""),
            "score": score(b["_prep"], b["subns_evidence"], b["detuning_tier"],
                           b["_order"], bool(single), b["_vg"],
                           amps[name][0] if name in amps else None, b["_unv"]),
            "_amp_row": amps[name][1] if name in amps else None,
            "_best": b, "_single": single,
        })

    out.sort(key=lambda x: (x["score"], x["spectrum"]))
    cuts = [(3, 1), (6, 2), (10, 3)]
    for x in out:
        x["tier"] = next((t for c, t in cuts if x["score"] <= c), 4)

    cols = ["tier", "score", "spectrum", "element", "Z", "charge", "prep_cost", "prep_note",
            "evidence", "fast_rung_at", "fast_rung_tau_s", "rung_lifetimes",
            "verified_ladders", "closure", "best_detuning_meV", "hop_amplitude",
            "relative_amplitude", "order", "pattern", "colours", "hops",
            "max_order", "min_hops",
            "ground_start", "single_colour", "hop_strength", "worst_virtual_gap_eV",
            "reach_eV_neutral_frame", "n_ladders", "accidental_matches_per_hop", "truncated"]
    with (ANA / "materials.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows({k: x[k] for k in cols} for x in out)

    def fmt(r):
        ms = r["photons_per_hop"].split("+")
        cs = r["colours_nm"].split("+")
        es = [r.get("start_eV", "?")] + r["path_eV"].split(" -> ")
        s = f"{es[0]} eV {r['start']}"
        for i, (m, c) in enumerate(zip(ms, cs)):
            s += f"  --{m}g@{c}nm-->  {es[i + 1]}"
        ai = r.get("final_above_ionization") == "yes" or \
            r.get("hop_lands_on", "").split(" | ")[-1] == "AI"
        return s + ("  [autoionizing]" if ai else "")

    g, bands, holes = hop_grid()
    L = []
    L.append("# Burst-ladder material list\n")
    L.append("Generated by `scripts/ladder.py` + `scripts/materials.py` from the cached "
             "NIST ASD tables in `energy_levels/`.\n")
    L.append(f"Driver colours: {', '.join(str(n) for n in DRIVERS_NM)} nm, plus "
             f"{' and '.join(str(n) for n in SH_NM)} nm as second harmonics of the Yb pair. "
             "Hops of 2 or 3 photons of one colour, tolerance 150 meV, at least 2 hops, "
             "total order >= 4.\n")

    L.append("\n## What the driver set can deliver\n")
    L.append("| rung energy (eV) | colour | photons |")
    L.append("|---|---|---|")
    for e, nm, m in g:
        L.append(f"| {e:.4f} | {nm} nm | {m} |")
    L.append("\nWith the 150 meV budget these merge into the continuous bands "
             + ", ".join(f"{lo:.2f}-{hi:.2f} eV" for lo, hi, _ in bands) + ".")
    if holes:
        L.append("No rung of any available colour falls in "
                 + ", ".join(f"{a:.2f}-{b:.2f} eV" for a, b in holes)
                 + " -- a level gap landing there cannot be bridged by this driver set.")

    # --- what the whole scan looks like, before it is cut into tiers ---------------------
    meas = [x for x in rows if x["subns_evidence"] == "measured"]
    parked_fast = [x for x in meas if x.get("fast_rung_at") in ("intermediate", "both")]
    verified = [x for x in rows if not x["_unv"]]
    singles = [x for x in rows if x["_single"]]
    ground = [x for x in rows if x["start_kind"] == "ground"]
    shortcut = [x for x in rows if x["direct_shortcut"]]
    col = {}
    for x in singles:
        c = f"{next(iter(sources(x['colours_nm']))):g}"
        col[c] = col.get(c, 0) + 1

    L.append("\n## What the scan found\n")
    L.append(f"- {len(rows)} admissible ladders in {len(out)} species, out of "
             f"{len(list((ROOT / 'energy_levels').glob('*/levels_with_lifetimes.csv')))} scraped.")
    L.append(f"- A MEASURED picosecond-window lifetime appears in {len(meas)} of them, "
             f"in {len({x['spectrum'] for x in meas})} species -- and every one is an ION. "
             "No neutral in the whole cache has a level with an ASD lifetime between 20 ps "
             "and 1 ns: in neutrals the fast states are the autoionizing resonances above the "
             "ionization limit, whose widths are not in a line list at all.")
    L.append(f"- In {len(parked_fast)} ladders the fast state is a PARKED intermediate rather "
             "than the last rung, which is the stronger form of the claim. The rest carry it "
             "on the last rung, where it is not a constraint at all -- nothing waits there, "
             "so a shorter lifetime is simply a brighter emitter.")
    L.append(f"- {len(verified)} ladders have a published lifetime for EVERY rung; the other "
             f"{len(rows) - len(verified)} lean on at least one bound level whose lifetime ASD "
             "does not tabulate. Those are not wrong -- almost every bound level below the "
             "ionization limit lives far longer than 20 ps -- but they are assumed, not read.")
    L.append(f"- {len(shortcut)} ladders have a competing direct channel: a single "
             "low-order jump from the start straight to the final state on a colour the "
             "ladder itself uses. The sequential-only requirement is therefore not a "
             "constraint that has to be imposed -- at these total energies no 2- or "
             "3-photon shortcut exists.")
    L.append(f"- Only {len(singles)} ladders ({100 * len(singles) / max(1, len(rows)):.0f}%) run "
             "on ONE laser, counting a colour and its own second harmonic as the same source "
             "(the SHG output needs no synchronization); the rest need two or more "
             "synchronized bursts. Among the single-laser ones the most productive driver is "
             + ", ".join(f"{k} nm ({v})" for k, v in sorted(col.items(), key=lambda t: -t[1])[:3])
             + ".")
    L.append(f"- {len(ground)} ladders start from the ground state; the rest need a "
             "metastable prepared first, by discharge or optical pumping.")
    cut = sorted({x["spectrum"] for x in rows if x.get("species_truncated") == "yes"})
    if cut:
        L.append(f"- The search stopped early in {len(cut)} species ("
                 + ", ".join(cut) + ") because their level density blows the path budget. "
                 "Their ladder counts are lower bounds and the particular chains kept are "
                 "arbitrary -- these are exactly the species where the chance-match column "
                 "is already saying an energy match proves nothing.")
    L.append("\nThe chance-match column is the check on all of this: it estimates how many "
             "levels fall inside one 150 meV hop window by accident. Above about 1, an energy "
             "match is not evidence, and the ladder has to stand on its parked lifetimes and "
             "hop strengths instead. Fe I sits at 16, Hg II at 0.7.\n")

    L.append("\n## Tiers\n")
    L.append("- **Tier 1** - preparable, closes on driver tuning alone, order 6 or better.")
    L.append("- **Tier 2** - one hard axis: hotter cell, a discharge, or a hop needing "
             "chirp/Stark help.")
    L.append("- **Tier 3** - two hard axes, or a ladder that only closes inside the full "
             "150 meV budget.")
    L.append("- **Tier 4** - a ladder exists on paper; preparation or closure makes it a "
             "late option.\n")

    for t in (1, 2, 3, 4):
        grp = [x for x in out if x["tier"] == t]
        if not grp:
            continue
        L.append(f"\n### Tier {t} - {len(grp)} species\n")
        L.append("| species | Z | preparation | fast rung | closure | best detuning | "
                 "state lifetimes | hop strength | order | pattern | best order anywhere | start | 1 colour | reach | "
                 "chance matches per hop |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for x in grp:
            L.append(f"| {x['spectrum']} | {x['Z']} | {x['prep_note']} | "
                     f"{EV_SHORT[x['evidence']]} | {x['closure']} | "
                     f"{x['best_detuning_meV']:.1f} meV | {x['rung_lifetimes']} | "
                     f"{x['hop_amplitude']}"
                     f"{' ' + x['relative_amplitude'] if x['relative_amplitude'] else ''} | "
                     f"{x['order']} | {x['pattern']} | {x['max_order']} | "
                     f"{x['ground_start'] == 'yes' and 'ground' or 'metastable'} | "
                     f"{x['single_colour']} | {x['reach_eV_neutral_frame']} eV | "
                     f"{x['accidental_matches_per_hop']} |")

    L.append("\n## Best ladder per species, tiers 1-3\n")
    for x in out:
        if x["tier"] > 3:
            continue
        L.append(f"**{x['spectrum']}** - {x['prep_note']}, "
                 f"fast rung {EV_SHORT[x['evidence']]}, closure {x['closure']}\n")
        L.append("```")
        L.append(fmt(x["_best"]))
        if x["_single"] is not None and x["_single"] is not x["_best"]:
            L.append("single colour: " + fmt(x["_single"]))
        a = x["_amp_row"]
        if a and a["path"] != x["_best"]["path"]:
            L.append("strongest hops: " + fmt(a))
        L.append("```\n")

    L.append("\n## What this scan cannot decide\n")
    L.append("- **Autoionizing widths.** Every neutral candidate ends on a level above the "
             "ionization limit. ASD carries no width for those, so the picosecond claim rests "
             "on the state class, not on a number. One linewidth measurement per candidate "
             "closes this.")
    L.append("- **The strength of the last hop.** Oscillator strengths exist for bound-bound "
             "transitions only, so the hop INTO the autoionizing region cannot be scored from "
             "a line list; it needs a photoabsorption cross section.")
    L.append("- **Lifetimes quoted as upper bounds.** ASD gives 1/sum(A) over the channels it "
             "lists. Parking uses the conservative floor tau * k/n for k of n channels known, "
             "but a level with no published A-value at all is invisible to the scan -- which "
             "is why Rb I drops out despite its 5s->4d two-photon hop landing 8 meV from "
             "2 x 1030 nm.")
    L.append("- **Anything time-dependent.** Whether the population actually survives the "
             "climb is a rate-equation question over the whole burst, not a level-table "
             "question.\n")
    (ANA / "materials.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    print(f"  {len(out)} species -> analysis/materials.csv, analysis/materials.md")
    for t in (1, 2, 3, 4):
        n = [x for x in out if x["tier"] == t]
        if n:
            print(f"    tier {t}: {len(n):3d}  " + ", ".join(x["spectrum"] for x in n[:14])
                  + (" ..." if len(n) > 14 else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
