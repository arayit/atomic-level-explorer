#!/usr/bin/env python3
"""Pull energy levels + radiative lifetimes for a list of species into energy_levels/<Sym>_<Roman>/.

Three sources, in descending order of trustworthiness for our purpose:

  1. NIST ASD Levels  -> configuration, term, J, g, energy (eV), parity, ionization limits.
     Complete for essentially every spectrum, neutral or ionized. No lifetime column exists.
  2. NIST ASD Lines   -> Aki per transition. tau_k = 1 / sum(Aki over lines with k as upper level).
     J-resolved. Coverage is the binding constraint: many ions carry zero A-values.
  3. ExoMol LiDB      -> precomputed radiative tau. Neutrals + singly-charged only, and only
     for states with complete A-value data. J-lumped for most species. Used as cross-check
     and as fallback where ASD Lines has no Aki.

Autoionizing / predissociating lifetimes are in NO database -- they come from primary
linewidth papers (tau[ps] = 0.658 / Gamma[meV]). Levels above the first ionization limit are
flagged `above_ionization`, and their tau column is left empty on purpose, not filled with a
radiative value that does not apply.

Usage:
    python3 scripts/scrape.py                      # every species in scripts/species.txt
    python3 scripts/scrape.py "Xe III" "Ar II"     # just these
    python3 scripts/scrape.py --refresh            # ignore cache, re-download
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"
SPECIES_FILE = Path(__file__).resolve().parent / "species.txt"

# NIST returns 403 to the default Python-urllib User-Agent. This is the whole reason
# the manual Firefox export existed; with a UA set the CGI endpoint is fully scriptable
# and reproduces the manual export byte for byte.
USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
DELAY_S = 1.5  # between requests, to stay a polite client of two public services

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]


# --------------------------------------------------------------------------- fetching

def http_get(url: str, timeout: int = 120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if 400 <= exc.code < 500:      # a definitive answer; retrying cannot help
                raise
            if attempt == 2:
                raise
            print(f"      retry {attempt + 1}/2 after {exc}")
            time.sleep(3 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 2:
                raise
            print(f"      retry {attempt + 1}/2 after {exc}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("unreachable")


def cached_get(url: str, dest: Path, refresh: bool) -> str:
    """Download to dest unless it is already there. Returns the text either way."""
    if dest.exists() and not refresh:
        return dest.read_text(encoding="utf-8", errors="replace")
    text = http_get(url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text, encoding="utf-8")
    time.sleep(DELAY_S)
    return text


def asd_levels_url(spectrum: str) -> str:
    q = {
        "de": "0",
        "spectrum": spectrum,
        "units": "1",            # 1 = eV
        "format": "2",           # 2 = CSV
        "output": "0",           # 0 = entire table on one page
        "page_size": "15",
        "multiplet_ordered": "0",  # 0 = energy ordered (term ordering scrambles the ladder)
        "conf_out": "on", "term_out": "on", "level_out": "on", "unc_out": "on",
        "j_out": "on", "g_out": "on", "perc_out": "on",
        "biblio": "on",          # measured-vs-calculated provenance
        "temp": "",
        "submit": "Retrieve Data",
    }
    return "https://physics.nist.gov/cgi-bin/ASD/energy1.pl?" + urllib.parse.urlencode(q)


def asd_lines_url(spectrum: str) -> str:
    q = {
        "spectra": spectrum,
        "limits_type": "0",
        "low_w": "", "upp_w": "",  # no wavelength window at all -- the 200 nm default
        "unit": "1",               # silently drops VUV/EUV decay channels and inflates tau
        "de": "0",
        "format": "2",             # 2 = CSV
        "line_out": "0",           # all lines
        "en_unit": "1",            # eV
        "output": "0",
        "bibrefs": "1",
        "page_size": "15",
        "show_obs_wl": "1", "show_calc_wl": "1", "unc_out": "1", "order_out": "0",
        "show_av": "3",            # 3 = vacuum wavelengths throughout (2 still returns air)
        "A_out": "0",              # include the Aki column -- the only reason we pull Lines
        "intens_out": "on",
        "allowed_out": "1", "forbid_out": "1",  # forbidden lines matter for metastables
        "conf_out": "on", "term_out": "on", "enrg_out": "on", "J_out": "on", "g_out": "on",
        "submit": "Retrieve Data",
    }
    return "https://physics.nist.gov/cgi-bin/ASD/lines1.pl?" + urllib.parse.urlencode(q)


def lidb_url(key: str, category: str) -> str:
    # The published LiDB API page documents only molecules and was never updated for the
    # 2025 atomic addition, but `molecule` works as the generic species key ("Ca", "Xe+").
    return "https://www.exomol.com/lidb/api/?" + urllib.parse.urlencode(
        {"molecule": key, "category": category, "format": "csv"})


# --------------------------------------------------------------------------- parsing

def unq(v) -> str:
    """ASD wraps every field as ="value" (Excel formula quoting); g is left bare."""
    s = (v or "").strip()
    m = re.match(r'^="(.*)"$', s)
    if m:
        s = m.group(1)
    return s.strip('"').strip()


def parse_energy(raw: str):
    """Energies carry [] for interpolated, () for theoretical, +x/+y for unknown offsets, ? for doubtful."""
    s = unq(raw).replace("[", "").replace("]", "").replace("(", "").replace(")", "")
    s = re.sub(r"\+\s*[xyz]", "", s).replace("?", "").replace("a", "").strip()
    m = re.search(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s)
    return float(m.group(0)) if m else None


def parse_j(raw: str):
    s = unq(raw)
    if not s or s == "---":
        return None
    if "," in s or "or" in s:      # unresolved multi-J entries
        return None
    if "/" in s:
        try:
            num, den = s.split("/")
            return float(num) / float(den)
        except ValueError:
            return None
    try:
        return float(s)
    except ValueError:
        return None


def parity_of(term: str):
    """'*' in the term string marks odd parity; this holds for LS and Racah notation alike."""
    t = unq(term)
    if not t or t == "Limit":
        return ""
    return "odd" if "*" in t else "even"


def norm_term(term: str) -> str:
    """Join key. LiDB writes 3Po for ASD's 3P*, and prefixes some terms with a lowercase letter."""
    t = unq(term).strip()
    t = re.sub(r"^[a-z]\s*", "", t)          # 'w 1P*' / 'w1Po' -> '1P*'
    if t.endswith("o"):                       # LiDB odd-parity marker
        t = t[:-1] + "*"
    return t.replace(" ", "")


def norm_conf(conf: str) -> str:
    return unq(conf).replace(" ", "").lower()


def parse_asd_levels(text: str):
    """Returns (levels, limits). Limit rows are ionization-series limits, not levels."""
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows or "Level (eV)" not in (rows[0] or {}):
        return [], []
    levels, limits = [], []
    for r in rows:
        term = unq(r.get("Term", ""))
        energy = parse_energy(r.get("Level (eV)", ""))
        if energy is None:
            continue
        conf = unq(r.get("Configuration", ""))
        if term == "Limit":
            limits.append({"label": conf, "energy_eV": energy})
            continue
        levels.append({
            "configuration": conf,
            "term": term,
            "J": parse_j(r.get("J", "")),
            "g": unq(r.get("g", "")),
            "energy_eV": energy,
            "unc_eV": parse_energy(r.get("Uncertainty (eV)", "")),
            "parity": parity_of(term),
            "doubtful": "yes" if unq(r.get("Suffix", "")) == "?" else "",
            "ref": unq(r.get("Reference", "")),
        })
    return levels, limits


def parse_asd_lines(text: str):
    """tau_k = 1 / sum(Aki) over every line with k as the upper level, keyed J-resolved.

    Also records how many of the level's known decay channels actually carry an A-value:
    a level whose channels are only partly covered has an OVERESTIMATED tau (upper bound),
    because the missing channels would only add to sum(Aki).
    """
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows or "Aki(s^-1)" not in (rows[0] or {}):
        return {}, 0, 0, 0
    agg: dict[tuple, dict] = {}
    n_ident = n_with_a = 0
    for r in rows:
        conf_k, term_k, j_k = unq(r.get("conf_k", "")), unq(r.get("term_k", "")), unq(r.get("J_k", ""))
        if not (conf_k or term_k):
            continue  # unclassified line: observed wavelength only, no level assignment, never an Aki
        n_ident += 1
        key = (norm_conf(conf_k), norm_term(term_k), parse_j(j_k))
        slot = agg.setdefault(key, {"sum_A": 0.0, "n_channels": 0, "n_with_A": 0,
                                    "strongest_A": 0.0, "strongest_dE": None, "max_dE": None})
        slot["n_channels"] += 1
        e_i, e_k = parse_energy(r.get("Ei(eV)", "")), parse_energy(r.get("Ek(eV)", ""))
        d_e = (e_k - e_i) if (e_i is not None and e_k is not None) else None
        if d_e is not None and (slot["max_dE"] is None or d_e > slot["max_dE"]):
            slot["max_dE"] = d_e
        a_raw = unq(r.get("Aki(s^-1)", ""))
        if not a_raw:
            continue
        try:
            a_ki = float(a_raw)
        except ValueError:
            continue
        n_with_a += 1
        slot["n_with_A"] += 1
        slot["sum_A"] += a_ki
        if a_ki > slot["strongest_A"]:
            slot["strongest_A"] = a_ki
            slot["strongest_dE"] = d_e
    out = {}
    for key, s in agg.items():
        if s["sum_A"] > 0:
            s["tau_s"] = 1.0 / s["sum_A"]
            out[key] = s
    return out, len(rows), n_ident, n_with_a


def parse_lidb_states(text: str):
    """LiDB labels look like 'Ca 3Po;v=3p6.4s.4p'. Lifetimes are J-lumped for most species.

    A blank lifetime means no allowed radiative channel, i.e. metastable -- that is data,
    not a gap, so it is kept as an explicit 'metastable' flag rather than dropped.
    """
    rows = list(csv.DictReader(io.StringIO(text)))
    if not rows or "Lifetime /s" not in (rows[0] or {}):
        return {}
    out = {}
    for r in rows:
        label = (r.get("State") or "").strip()
        m = re.match(r"^(\S+)\s+(.*?);v=(.*)$", label)
        if not m:
            continue
        _sp, term, conf = m.groups()
        raw = (r.get("Lifetime /s") or "").strip()
        try:
            tau = float(raw) if raw else None
        except ValueError:
            tau = None
        out[(norm_conf(conf), norm_term(term))] = {
            "tau_s": tau,
            "metastable": raw == "",
            "energy_eV": parse_energy(r.get("Energy /eV", "")),
        }
    return out


# --------------------------------------------------------------------------- per species

def species_paths(spectrum: str) -> Path:
    sym, roman = spectrum.split()
    return OUT_ROOT / f"{sym}_{roman}"


def lidb_key(spectrum: str):
    """LiDB covers neutrals and singly-charged ions only; charge is encoded as a trailing '+'."""
    sym, roman = spectrum.split()
    charge = ROMAN.index(roman)
    if charge == 0:
        return sym
    if charge == 1:
        return sym + "+"
    return None


def build_species(spectrum: str, refresh: bool) -> dict:
    d = species_paths(spectrum)
    d.mkdir(parents=True, exist_ok=True)
    sym, roman = spectrum.split()
    print(f"  {spectrum}")

    levels_txt = cached_get(asd_levels_url(spectrum), d / "asd_levels.csv", refresh)
    levels, limits = parse_asd_levels(levels_txt)
    print(f"      ASD levels : {len(levels):5d} levels, {len(limits)} ionization limits")

    lines_txt = cached_get(asd_lines_url(spectrum), d / "asd_lines.csv", refresh)
    tau_asd, n_lines, n_ident, n_with_a = parse_asd_lines(lines_txt)
    print(f"      ASD lines  : {n_lines:6d} lines ({n_ident} classified, {n_with_a} with Aki)"
          f" -> tau for {len(tau_asd)} levels")

    key = lidb_key(spectrum)
    lidb = {}
    if key:
        try:
            st = cached_get(lidb_url(key, "states"), d / "lidb_states.csv", refresh)
        except Exception:
            st = ""                       # absent from LiDB: reported below, not fatal
        if st.lstrip().startswith("State,"):
            lidb = parse_lidb_states(st)
            cached_get(lidb_url(key, "transitions"), d / "lidb_transitions.csv", refresh)
            print(f"      LiDB       : {len(lidb)} states")
        else:
            (d / "lidb_states.csv").unlink(missing_ok=True)
            print("      LiDB       : species not in database")
    else:
        print("      LiDB       : n/a (charge >= 2)")

    ion_eV = min((l["energy_eV"] for l in limits), default=None)

    for lv in levels:
        ck, tk = norm_conf(lv["configuration"]), norm_term(lv["term"])
        a = tau_asd.get((ck, tk, lv["J"]))
        lv["tau_asd_s"] = a["tau_s"] if a else None
        lv["asd_channels"] = f"{a['n_with_A']}/{a['n_channels']}" if a else ""
        lv["asd_tau_complete"] = ("yes" if a["n_with_A"] == a["n_channels"] else "no") if a else ""
        lv["strongest_decay_eV"] = a["strongest_dE"] if a else None
        lv["max_decay_eV"] = a["max_dE"] if a else None

        b = lidb.get((ck, tk))
        lv["tau_lidb_s"] = b["tau_s"] if b else None
        lv["lidb_metastable"] = "yes" if (b and b["metastable"]) else ""

        lv["above_ionization"] = "yes" if (ion_eV is not None and lv["energy_eV"] > ion_eV) else "no"
        lv["ionization_eV"] = ion_eV

        # Radiative tau is meaningless above the first ionization limit: those states
        # autoionize, and no database carries that width. Leave tau empty and say why.
        if lv["above_ionization"] == "yes":
            lv["tau_s"], lv["tau_source"] = None, "autoionizing - not in any database"
        elif lv["tau_asd_s"] is not None:
            lv["tau_s"], lv["tau_source"] = lv["tau_asd_s"], "ASD 1/sum(Aki), J-resolved"
        elif lv["tau_lidb_s"] is not None:
            lv["tau_s"], lv["tau_source"] = lv["tau_lidb_s"], "LiDB (J-lumped)"
        elif b and b["metastable"]:
            lv["tau_s"], lv["tau_source"] = None, "metastable - no allowed radiative channel"
        else:
            lv["tau_s"], lv["tau_source"] = None, ""

    cols = ["spectrum", "element", "charge", "configuration", "term", "J", "g", "energy_eV",
            "unc_eV", "parity", "doubtful", "above_ionization", "ionization_eV",
            "tau_s", "tau_source", "tau_asd_s", "asd_channels", "asd_tau_complete",
            "strongest_decay_eV", "max_decay_eV", "tau_lidb_s", "lidb_metastable", "ref"]
    charge = ROMAN.index(roman)
    with (d / "levels_with_lifetimes.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for lv in levels:
            w.writerow({**{c: "" for c in cols}, **lv,
                        "spectrum": spectrum, "element": sym, "charge": charge})

    with (d / "limits.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "energy_eV"])
        w.writeheader()
        w.writerows(limits)

    n_tau = sum(1 for lv in levels if lv["tau_s"] is not None)
    n_sub_ns = sum(1 for lv in levels if lv["tau_s"] and 1e-12 < lv["tau_s"] < 1e-9)
    n_above = sum(1 for lv in levels if lv["above_ionization"] == "yes")
    summary = {
        "spectrum": spectrum, "element": sym, "charge": charge,
        "n_levels": len(levels), "ionization_eV": ion_eV,
        "max_level_eV": max((lv["energy_eV"] for lv in levels), default=None),
        "n_levels_above_ionization": n_above,
        "n_lines": n_lines, "n_lines_classified": n_ident, "n_lines_with_Aki": n_with_a,
        "n_levels_with_tau": n_tau,
        "n_tau_from_asd": sum(1 for lv in levels if lv["tau_source"].startswith("ASD")),
        "n_tau_from_lidb": sum(1 for lv in levels if lv["tau_source"].startswith("LiDB")),
        "n_metastable": sum(1 for lv in levels if lv["lidb_metastable"] == "yes"),
        "n_sub_ns": n_sub_ns,
        "tau_coverage": round(n_tau / len(levels), 3) if levels else 0.0,
    }
    (d / "meta.json").write_text(json.dumps(
        {**summary,
         "sources": {"asd_levels": asd_levels_url(spectrum), "asd_lines": asd_lines_url(spectrum),
                     "lidb_states": lidb_url(key, "states") if key else None}},
        indent=2), encoding="utf-8")
    print(f"      -> {n_tau}/{len(levels)} levels with tau, {n_sub_ns} in the 1 - 1000 ps window")
    return summary


# --------------------------------------------------------------------------- main

def read_species_file() -> list[str]:
    out = []
    for line in SPECIES_FILE.read_text(encoding="utf-8").splitlines():
        line = line.split("#")[0].strip()
        if line:
            out.append(line)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("species", nargs="*", help='e.g. "Xe III" (default: scripts/species.txt)')
    ap.add_argument("--refresh", action="store_true", help="ignore the on-disk cache")
    args = ap.parse_args()

    todo = args.species or read_species_file()
    print(f"Scraping {len(todo)} species into {OUT_ROOT}\n")

    summaries, failed = [], []
    for spectrum in todo:
        try:
            summaries.append(build_species(spectrum, args.refresh))
        except Exception as exc:                       # keep going; one dead species is not fatal
            print(f"      FAILED: {exc}")
            failed.append((spectrum, str(exc)))

    if summaries:
        cov = OUT_ROOT / "coverage.csv"
        existing = {}
        if cov.exists():
            existing = {r["spectrum"]: r for r in csv.DictReader(cov.open(encoding="utf-8"))}
        for s in summaries:
            existing[s["spectrum"]] = s
        rows = sorted(existing.values(), key=lambda r: (str(r["element"]), int(r["charge"])))
        with cov.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summaries[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nCoverage table -> {cov}")

    if failed:
        print("\nFailed:")
        for sp, err in failed:
            print(f"  {sp}: {err}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
