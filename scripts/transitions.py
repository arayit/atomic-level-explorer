#!/usr/bin/env python3
"""Derive absorption strengths from the A-values already on disk.

No new download: the lifetime and the absorption cross section come from the same matrix
element, so the ASD Lines tables we cached for tau also carry everything needed for f.

    f_ik   = 1.4992e-16 * (g_k / g_i) * A_ki * lambda[A]^2        (lambda in vacuum)
    gf     = g_i * f_ik
    int sigma dnu = (pi e^2 / m_e c) * f_ik = 2.654e-2 * f_ik      cm^2 Hz

The integrated cross section is quoted rather than a peak value because the peak depends on the
line profile, and therefore on temperature, pressure and geometry -- assumptions that belong to
an experiment, not to a table. Divide by the linewidth in Hz to get cm^2.

IMPORTANT: these are ONE-photon quantities. The ladder's m = 2 and m = 3 hops need
sigma^(2) and sigma^(3), which are sums over intermediate states and are not tabulated anywhere.

Usage:
    python3 scripts/transitions.py            # every scraped species
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrape import unq, parse_energy, parse_j, norm_conf, norm_term   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"

HC = 1239.841984          # eV nm
F_CONST = 1.4992e-16      # f = F_CONST * (gk/gi) * A * lambda[A]^2
SIGMA_CONST = 2.654e-2    # int sigma dnu = SIGMA_CONST * f, in cm^2 Hz

COLS = ["conf_i", "term_i", "J_i", "g_i", "E_i_eV",
        "conf_k", "term_k", "J_k", "g_k", "E_k_eV",
        "dE_eV", "lambda_vac_nm", "Aki_s-1", "f_ik", "gf", "sigma_int_cm2Hz"]


def g_of(raw):
    s = unq(raw)
    try:
        v = float(s)
        return v if v > 0 else None
    except ValueError:
        return None


def build(d: Path) -> tuple[int, int]:
    src = d / "asd_lines.csv"
    if not src.exists():
        return 0, 0
    rows = list(csv.DictReader(src.open(encoding="utf-8")))
    out, skipped = [], 0
    for r in rows:
        a_raw = unq(r.get("Aki(s^-1)", ""))
        if not a_raw:
            continue
        try:
            a_ki = float(a_raw)
        except ValueError:
            continue
        e_i, e_k = parse_energy(r.get("Ei(eV)", "")), parse_energy(r.get("Ek(eV)", ""))
        g_i, g_k = g_of(r.get("g_i", "")), g_of(r.get("g_k", ""))
        if None in (e_i, e_k, g_i, g_k) or e_k <= e_i:
            skipped += 1          # degeneracies or level energies missing: f is not defined
            continue
        d_e = e_k - e_i
        lam_nm = HC / d_e
        f_ik = F_CONST * (g_k / g_i) * a_ki * (lam_nm * 10.0) ** 2
        out.append({
            "conf_i": unq(r.get("conf_i", "")), "term_i": unq(r.get("term_i", "")),
            "J_i": parse_j(r.get("J_i", "")), "g_i": g_i, "E_i_eV": round(e_i, 6),
            "conf_k": unq(r.get("conf_k", "")), "term_k": unq(r.get("term_k", "")),
            "J_k": parse_j(r.get("J_k", "")), "g_k": g_k, "E_k_eV": round(e_k, 6),
            "dE_eV": round(d_e, 6), "lambda_vac_nm": round(lam_nm, 5),
            "Aki_s-1": a_ki, "f_ik": float(f"{f_ik:.5g}"),
            "gf": float(f"{g_i * f_ik:.5g}"),
            "sigma_int_cm2Hz": float(f"{SIGMA_CONST * f_ik:.5g}"),
        })
    out.sort(key=lambda x: (x["E_i_eV"], x["dE_eV"]))
    with (d / "transitions.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(out)
    return len(out), skipped


def main() -> int:
    tot = skip = done = 0
    for d in sorted(OUT_ROOT.iterdir()):
        if not d.is_dir() or not (d / "asd_lines.csv").exists():
            continue
        n, s = build(d)
        tot += n
        skip += s
        done += 1
    print(f"  {done} species -> transitions.csv")
    print(f"  {tot} transitions with an oscillator strength")
    print(f"  {skip} lines carried an A-value but lacked g or level energies, so f was left out")
    return 0


if __name__ == "__main__":
    sys.exit(main())
