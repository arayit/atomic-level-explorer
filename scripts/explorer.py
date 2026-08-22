#!/usr/bin/env python3
"""Build the interactive level explorer from the scraped CSVs.

Writes two files from one template:
  index.html                           the site: a standalone page, also served by Pages
  <scratch>/explorer_artifact.html     the same page as a body fragment, for publishing

The point of the page is the question a static figure cannot answer: a level sits at 91.8 eV
above ITS OWN ion's ground state, but that ion had to be made first. The reference toggle adds
the cumulative ionization cost, so a level can be read in the frame the laser actually works
in -- above the neutral atom.

Usage:
    python3 scripts/explorer.py [output.html]
"""

from __future__ import annotations

import csv
import datetime
import glob
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scrape import norm_conf, norm_term, parse_j            # noqa: E402
import ladders_tab                                          # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT / "energy_levels"
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX"]

WIN_LO, WIN_HI = 1e-12, 1e-9          # the picosecond decade
LAMBDA_NM = 1030                       # Yb fibre driver, default



ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19,
    "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
    "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
    "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "W": 74, "Au": 79,
    "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83,
}
# Periodic-table category. "metalloid" is the semiconductor group (B, Si, Ge, As, Sb, Te).
ELEMENT_CAT = {
    **{e: "noble gas" for e in ("He", "Ne", "Ar", "Kr", "Xe")},
    **{e: "alkali metal" for e in ("Li", "Na", "K", "Rb", "Cs")},
    **{e: "alkaline earth" for e in ("Be", "Mg", "Ca", "Sr", "Ba")},
    **{e: "nonmetal" for e in ("H", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I")},
    **{e: "metalloid" for e in ("B", "Si", "Ge", "As", "Sb", "Te")},
    **{e: "transition metal" for e in ("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                                       "Zn", "Y", "Zr", "Nb", "Mo", "Ag", "Cd", "La", "W",
                                       "Au", "Hg")},
    **{e: "post-transition metal" for e in ("Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi")},
}


# Temperature at which the element reaches 1 Pa of vapour (about 1e14 cm^-3, a working
# vapour-cell density), in kelvin; None where the source tabulates no 1 Pa point.
# CRC Handbook of Chemistry and Physics, 84th ed.; metals after Alcock, Itkin & Horrigan (1984).
VAPOUR_1PA_K = {
    "H": 10, "He": 1.3, "Li": 797, "Be": 1462, "B": 2348, "C": None, "N": 37, "O": None,
    "F": 38, "Ne": 12, "Na": 554, "Mg": 701, "Al": 1482, "Si": 1908, "P": 279, "S": 375,
    "Cl": 128, "Ar": 47, "K": 473, "Ca": 864, "Ti": 1982, "Fe": 1728, "Kr": 59, "Rb": 434,
    "Sr": 796, "Cd": 530, "Xe": 83, "Cs": 418, "Ba": 911, "Hg": 315,
    "V": 2101, "Cr": 1656, "Mn": 1228, "Co": 1790, "Ni": 1783, "Cu": 1509, "Zn": 610,
    "Ga": 1310, "Ge": 1644, "In": 1196, "Sn": 1497, "Tl": 882, "Pb": 978,
    "Y": None, "Zr": None,          # source tabulates no measured point for these
}
# Carbon sublimes; the source tabulates no 1 Pa point, so its 100 Pa point is quoted instead.
VAPOUR_100PA_K = {"C": 2839}

# Boils below room temperature at 1 atm, i.e. already a gas on the bench.
GAS_AT_RT = {"H", "He", "N", "O", "F", "Ne", "Cl", "Ar", "Kr", "Xe"}


def vapour_note(el: str) -> str:
    """Condition for getting the PARENT ELEMENT into the gas phase -- a precondition for every
    charge state of it, not a property of the ion itself."""
    if el in GAS_AT_RT:
        return "already gas at 20 °C"
    t = VAPOUR_1PA_K.get(el)
    if t:
        return f"1 Pa at {t - 273:.0f} °C"
    t = VAPOUR_100PA_K.get(el)
    return f"100 Pa at {t - 273:.0f} °C" if t else ""


def sig(x, n=4):
    if x is None:
        return None
    return float(f"{x:.{n}g}")


def collect():
    species = {}
    for d in sorted(OUT_ROOT.iterdir()):
        f = d / "levels_with_lifetimes.csv"
        if not d.is_dir() or not f.exists():
            continue
        sp = d.name.replace("_", " ")
        el, roman = sp.split()
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        confs, terms = {}, {}
        confs_rev, terms_rev = [], []
        rows = []
        for r in csv.DictReader(f.open(encoding="utf-8")):
            if not r["energy_eV"]:
                continue
            if r["configuration"] not in confs:
                confs[r["configuration"]] = len(confs); confs_rev.append(r["configuration"])
            if r["term"] not in terms:
                terms[r["term"]] = len(terms); terms_rev.append(r["term"])
            ci, ti = confs[r["configuration"]], terms[r["term"]]
            tau = float(r["tau_s"]) if r["tau_s"] else None
            rows.append([
                ci, ti,
                float(r["J"]) if r["J"] else None,
                sig(float(r["energy_eV"]), 8),
                sig(tau, 3),
                sig(float(r["strongest_decay_eV"]), 4) if r["strongest_decay_eV"] else None,
                1 if r["parity"] == "odd" else 0,
                1 if r["above_ionization"] == "yes" else 0,
                1 if r["lidb_metastable"] == "yes" else 0,
            ])
        # Absorption channels: the same A-values that gave tau also give the oscillator
        # strength, so each level can carry the strongest transitions leading up out of it.
        index = {}
        for i, r in enumerate(rows):
            index[(norm_conf(confs_rev[r[0]]), norm_term(terms_rev[r[1]]), r[2])] = i
        chan = {}
        tf = d / "transitions.csv"
        if tf.exists():
            for t in csv.DictReader(tf.open(encoding="utf-8")):
                fv = float(t["f_ik"])
                if fv < 0.01:
                    continue
                lo = index.get((norm_conf(t["conf_i"]), norm_term(t["term_i"]),
                                float(t["J_i"]) if t["J_i"] else None))
                hi = index.get((norm_conf(t["conf_k"]), norm_term(t["term_k"]),
                                float(t["J_k"]) if t["J_k"] else None))
                if lo is None or hi is None or lo == hi:
                    continue
                chan.setdefault(lo, []).append([hi, float(f"{fv:.3g}")])
        chan = {k: sorted(v, key=lambda x: -x[1])[:6] for k, v in chan.items()}

        limits = []
        if (d / "limits.csv").exists():
            for r in csv.DictReader((d / "limits.csv").open(encoding="utf-8")):
                try:
                    limits.append([r["label"], sig(float(r["energy_eV"]), 7)])
                except (TypeError, ValueError):
                    pass
        species[sp] = {
            "sp": sp, "el": el, "q": ROMAN.index(roman),
            "ion": sig(meta.get("ionization_eV"), 7),
            "confs": list(confs), "terms": list(terms),
            "lv": rows, "limits": limits, "abs": chan,
            "nLines": meta.get("n_lines", 0), "nAki": meta.get("n_lines_with_Aki", 0),
            "vap": vapour_note(el),
        }

    # Absolute frame: a level of El^q+ sits this far above the NEUTRAL ground state.
    ions = {s: species[s]["ion"] for s in species}
    for sp, rec in species.items():
        total, ok = 0.0, True
        for k in range(rec["q"]):
            prev = ions.get(f"{rec['el']} {ROMAN[k]}")
            if prev is None:
                ok = False
                break
            total += prev
        rec["off"] = sig(total, 8) if ok else None
    return species


def element_table(species) -> dict:
    els = {r["el"] for r in species.values()}
    return {e: [ELEMENT_Z.get(e, 999), ELEMENT_CAT.get(e, "")] for e in els}


CSS = r"""
:root{
  color-scheme:light;
  --paper:#ffffff; --panel:#ffffff; --panel-2:#f2f2f0;
  --ink:#000000; --ink-2:#2e2e2e; --ink-3:#606060;
  --rule:#9a9a9a; --rule-2:#d6d6d6;
  --accent:#a4001e; --accent-ink:#a4001e; --accent-soft:#f6e9eb;
  --link:#00008b;
  --warm:#5a5a5a; --band:#efeeea;
  --plain:#b6b6b6; --known:#000000;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    color-scheme:dark;
    --paper:#0f0f0f; --panel:#141414; --panel-2:#1c1c1c;
    --ink:#ececec; --ink-2:#c4c4c4; --ink-3:#8c8c8c;
    --rule:#4a4a4a; --rule-2:#2c2c2c;
    --accent:#ff8080; --accent-ink:#ff9a9a; --accent-soft:#2a1618;
    --link:#8fb6ff;
    --warm:#9a9a9a; --band:#1b1a18;
    --plain:#484848; --known:#e4e4e4;
  }
}
:root[data-theme="dark"]{
  color-scheme:dark;
  --paper:#0f0f0f; --panel:#141414; --panel-2:#1c1c1c;
  --ink:#ececec; --ink-2:#c4c4c4; --ink-3:#8c8c8c;
  --rule:#4a4a4a; --rule-2:#2c2c2c;
  --accent:#ff8080; --accent-ink:#ff9a9a; --accent-soft:#2a1618;
  --link:#8fb6ff;
  --warm:#9a9a9a; --band:#1b1a18;
  --plain:#484848; --known:#e4e4e4;
}

*{box-sizing:border-box}
body{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:Georgia,"Times New Roman",Times,serif; font-size:14px; line-height:1.45;
}
.mono{font-family:ui-monospace,Menlo,Consolas,"Courier New",monospace; font-variant-numeric:tabular-nums}
a{color:var(--link)}
.app{display:flex; flex-direction:column; min-height:100vh}

/* ---------- head + controls ---------- */
.head{padding:14px 18px 10px}
.head h1{margin:0; font-size:19px; font-weight:normal; letter-spacing:.01em}
.head p{margin:4px 0 0; font-size:13px; color:var(--ink-2); max-width:96ch}
.controls{
  display:flex; flex-wrap:wrap; align-items:baseline; gap:6px 20px;
  padding:9px 18px; border-top:1px solid var(--rule); border-bottom:1px solid var(--rule);
  background:var(--panel-2); font-size:13px;
}
.controls label{color:var(--ink-2)}
.controls select,.controls input,.controls button{font:inherit; font-size:13px; color:var(--ink)}
.controls input[type=number]{width:74px; text-align:right; padding:1px 3px;
  font-family:ui-monospace,Menlo,Consolas,monospace}
.controls .derived{font-size:12px; color:var(--ink-3)}
.controls button{padding:1px 8px}

/* ---------- grid ---------- */
#view-levels{display:flex; flex-direction:column; flex:1; min-height:0}
/* an id selector outranks the user agent's [hidden] rule, so say it again here */
#view-levels[hidden]{display:none}
.grid{display:grid; grid-template-columns:196px minmax(0,1fr) 268px; flex:1; min-height:0}
.rail{border-right:1px solid var(--rule); overflow:auto; max-height:calc(100vh - 148px)}
.rail.right{border-right:0; border-left:1px solid var(--rule)}
.railhead{position:sticky; top:0; background:var(--paper); padding:8px 12px 5px;
  border-bottom:1px solid var(--rule-2); font-size:12px; color:var(--ink-3); z-index:2}
.elname{padding:8px 12px 3px; font-size:12px; color:var(--ink-3)}
.elname b{color:var(--ink); font-weight:bold; font-style:normal; font-size:13px}
.elname .z{margin-left:6px; font-family:ui-monospace,Menlo,Consolas,monospace; font-size:11px}
.elname .cat{display:block; font-style:italic; font-size:11px; margin-top:-1px}
.spbtn{display:flex; align-items:baseline; gap:8px; width:100%; text-align:left;
  border:0; background:transparent; color:var(--ink-2); font:inherit; font-size:13px;
  padding:2px 12px; cursor:pointer}
.spbtn:hover{background:var(--panel-2); color:var(--ink)}
.spbtn[aria-current="true"]{color:var(--ink); font-weight:bold; background:var(--panel-2)}
.spbtn .n{margin-left:auto; font-size:12px; color:var(--ink-3);
  font-family:ui-monospace,Menlo,Consolas,monospace}
.spbtn .n.hit{color:var(--accent-ink)}

/* ---------- stage ---------- */
.stage{position:relative; min-width:0; display:flex; flex-direction:column}
.canvaswrap{position:relative; flex:1; min-height:520px}
canvas{display:block; width:100%; height:100%; cursor:crosshair; touch-action:none}
.hint{padding:5px 16px; border-top:1px solid var(--rule); font-size:12px; color:var(--ink-3);
  display:flex; gap:16px; flex-wrap:wrap}
.tip{position:absolute; pointer-events:none; z-index:5; max-width:290px;
  background:var(--paper); border:1px solid var(--rule); padding:6px 9px; font-size:12.5px; opacity:0}
.tip.on{opacity:1}
.tip b{display:block; margin-bottom:3px; font-weight:bold}
.tip .row{display:flex; justify-content:space-between; gap:14px; color:var(--ink-2)}
.tip .row span:last-child{color:var(--ink)}

/* ---------- detail ---------- */
.detail{padding:10px 14px 20px}
.card{border-top:1px solid var(--rule-2); padding:9px 0}
.card:first-child{border-top:0}
.card h3{margin:0 0 5px; font-size:12px; color:var(--ink-3); font-weight:normal; font-style:italic}
.card h3.name{font-style:normal; font-weight:bold; font-size:15px; color:var(--ink)}
.kv{display:grid; grid-template-columns:auto 1fr; gap:1px 12px; font-size:12.5px}
.kv dt{color:var(--ink-3)}
.kv dd{margin:0; text-align:right; color:var(--ink)}
.empty{color:var(--ink-3); font-size:12.5px; font-style:italic}
.mark{color:var(--accent-ink); font-weight:bold; font-size:12.5px}
table.abs{border-collapse:collapse; width:100%; font-size:12px; margin-top:2px}
table.abs th{text-align:right; font-weight:normal; font-style:italic; color:var(--ink-3);
  padding:0 0 2px 8px; border-bottom:1px solid var(--rule-2)}
table.abs th:first-child,table.abs td:first-child{text-align:left; padding-left:0}
table.abs td{text-align:right; padding:1px 0 1px 8px; white-space:nowrap}
.card .note{margin:6px 0 0; font-size:11px; color:var(--ink-3); line-height:1.45}
.btn{font:inherit; font-size:12.5px; padding:1px 8px}

/* ---------- table ---------- */
.tablewrap{border-top:1px solid var(--rule)}
.tablehead{display:flex; align-items:baseline; gap:14px; padding:8px 18px; flex-wrap:wrap}
.tablehead h2{margin:0; font-size:14px; font-weight:normal; font-style:italic}
.tablehead .muted{font-size:12px; color:var(--ink-3)}
.scroll{overflow:auto; max-height:320px; border-top:1px solid var(--rule-2)}
table{border-collapse:collapse; width:100%; font-size:12.5px}
th,td{padding:2px 12px; text-align:right; white-space:nowrap}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}
thead th{position:sticky; top:0; background:var(--paper); z-index:1; color:var(--ink-2);
  font-weight:normal; font-style:italic; border-bottom:1px solid var(--rule); cursor:pointer}
tbody tr:hover{background:var(--panel-2)}
tbody tr.win td{color:var(--accent-ink)}
tbody tr.sel td{font-weight:bold}

/* ---------- sources ---------- */
.refs{border-top:1px solid var(--rule); padding:16px 18px 26px}
.refs h2{margin:0 0 8px; font-size:14px; font-weight:normal; font-style:italic}
.refs p{margin:0 0 6px; max-width:104ch; font-size:12.5px; color:var(--ink-2); line-height:1.5}
.refs p b{color:var(--ink); font-weight:bold}
.refs code{font-family:ui-monospace,Menlo,Consolas,monospace; font-size:11.5px}

@media (max-width:1080px){
  .grid{grid-template-columns:1fr}
  .rail{max-height:none; border-right:0; border-bottom:1px solid var(--rule)}
  .rail.right{border-left:0; border-top:1px solid var(--rule)}
  .canvaswrap{min-height:460px}
}
@media (prefers-reduced-motion:reduce){*{transition:none!important; animation:none!important}}
"""


BODY = r"""
<div class="app">
  <div class="head">
    <h1>Atomic level explorer</h1>
    <p>Energy levels and radiative lifetimes for 104 spectra, computed from NIST ASD transition
       probabilities. Levels whose lifetime falls between 1 and 1000 ps are marked throughout.
       Sources and method are listed at the foot of the page.</p>
  </div>

  <nav class="tabs" role="tablist">
    <button type="button" id="tab-levels" role="tab" aria-selected="true"
            aria-controls="view-levels">Levels</button>
    <button type="button" id="tab-ladders" role="tab" aria-selected="false"
            aria-controls="view-ladders">Ladders</button>
  </nav>

  <section id="view-levels">
  <form class="controls" onsubmit="return false">
    <label>Energy from
      <select id="frame">
        <option value="ion" selected>this ion&rsquo;s ground state</option>
        <option value="abs">the neutral atom&rsquo;s ground state</option>
      </select>
    </label>
    <label>Photon wavelength
      <input id="lam" type="number" min="1" max="2000" step="any" value="1030" list="lampresets"
             title="1 to 2000 nm" aria-label="photon wavelength in nanometres"> nm
    </label>
    <datalist id="lampresets">
      <option value="13.5" label="EUV lithography"></option>
      <option value="515" label="Yb second harmonic"></option>
      <option value="800" label="Ti:sapphire"></option>
      <option value="1030" label="Yb fibre"></option>
      <option value="1064" label="Nd:YAG"></option>
    </datalist>
    <span class="derived mono" id="lamout"></span>
    <label>Fine structure
      <select id="fine">
        <option value="0" selected>collapsed</option>
        <option value="1">resolved</option>
      </select>
    </label>
    <label>Show
      <select id="filt">
        <option value="all" selected>all levels</option>
        <option value="tau">levels with a lifetime</option>
        <option value="win">1&ndash;1000 ps only</option>
      </select>
    </label>
    <button type="button" id="reset">Reset view</button>
  </form>

  <div class="grid">
    <aside class="rail" id="rail">
      <div class="railhead">Spectrum &mdash; levels in 1&ndash;1000 ps</div>
      <div id="splist"></div>
    </aside>

    <main class="stage">
      <div class="canvaswrap" id="wrap">
        <canvas id="cv"></canvas>
        <div class="tip" id="tip"></div>
      </div>
      <div class="hint">
        <span>scroll to zoom</span><span>drag to pan</span>
        <span>click a level to measure every gap from it</span><span>esc clears</span>
        <span id="zoomnote"></span>
      </div>
    </main>

    <aside class="rail right detail" id="detail"></aside>
  </div>

  <section class="tablewrap">
    <div class="tablehead">
      <h2 id="tabtitle">Levels in view</h2>
      <span class="muted" id="tabnote"></span>
    </div>
    <div class="scroll"><table>
      <thead><tr>
        <th data-k="conf">Configuration</th><th data-k="term">Term</th><th data-k="J">J</th>
        <th data-k="E">Energy eV</th><th data-k="tau">Lifetime</th>
        <th data-k="decay">Strongest decay eV</th><th data-k="par">Parity</th>
      </tr></thead>
      <tbody id="tbody"></tbody>
    </table></div>
  </section>

  </section>

__LADDERS_BODY__

  <section class="refs">
    <h2>Sources</h2>
    <p><b>Energy levels and transition probabilities.</b> NIST Atomic Spectra Database v5.12,
      retrieved __DATE__. Kramida A., Ralchenko Yu., Reader J. and NIST ASD Team (2024),
      <a href="https://doi.org/10.18434/T4W30F" target="_blank" rel="noopener">doi:10.18434/T4W30F</a>.
      Cross-checked against ExoMol LiDB &mdash; Owens A. et al. (2025), <i>JQSRT</i> <b>330</b>,
      109242 &mdash; which is itself derived from ASD, so not an independent check.</p>
    <p><b>Lifetimes are not tabulated anywhere.</b> They are computed here as
      <span class="mono">&tau; = 1 / &Sigma; A<sub>ki</sub></span> over every listed decay channel,
      J-resolved, and are upper bounds where a level&rsquo;s channels are only partly covered.
      Levels above the first ionization limit autoionize and are left blank.
      Oscillator strengths come from the same A-values,
      <span class="mono">f = 1.4992&times;10<sup>&minus;16</sup> (g<sub>k</sub>/g<sub>i</sub>)
      A<sub>ki</sub> &lambda;<sup>2</sup></span> with &lambda; in &aring;ngstr&ouml;m, and are
      one-photon quantities only.
      The Ladders tab is built by <span class="mono">scripts/ladder.py</span> from these same
      tables: hops of 2 or 3 photons of one colour out of 1030, 1060, 1550, 1600, 1900 and
      2000&nbsp;nm, matched to a measured level gap within 150&nbsp;meV, parity and
      |&Delta;J|&nbsp;&le;&nbsp;m enforced per hop, and every rung but the last required to
      outlive a 20&nbsp;ps intra-burst spacing. Hop amplitudes are perturbative dominant-term
      estimates built from the same oscillator strengths, so they rank hops and do not measure
      them, and no oscillator strength exists for the hop into the autoionizing region.
      Vapour temperatures belong to the neutral element, not to the ion, and are a precondition
      for preparing any charge state of it: the 1 Pa points (about
      10<sup>14</sup> cm<sup>&minus;3</sup>) from the
      <i>CRC Handbook of Chemistry and Physics</i>, 84th ed.; metals after Alcock, Itkin &amp;
      Horrigan (1984).</p>
  </section>
</div>
"""


JS = r"""
const DATA = __DATA__;
const ELEM = __ELEM__;
const WIN_LO = __WIN_LO__, WIN_HI = __WIN_HI__;
const HC = 1239.841984;              // eV nm
const LAM_MIN = 1, LAM_MAX = 2000;
const nmLabel = v => v >= 100 ? String(Math.round(v)) : String(+v.toFixed(1));
const eVLabel = v => v >= 100 ? v.toFixed(1) : v >= 10 ? v.toFixed(2) : v.toFixed(3);
let lamNm = __LAMBDA__;
let PHOTON = HC / lamNm, PHOTON_SH = 2 * PHOTON;

// Every (n1030, n515) with m = n1030 + n515 <= 3, grouped by the energy it delivers.
// The same energy can arrive with two different photon counts -- that is the parity switch.
function buildHops(){
  const m = new Map();
  for (let a = 0; a <= 3; a++) for (let b = 0; b <= 3; b++){
    const n = a + b;
    if (!n || n > 3) continue;
    const e = a * PHOTON + b * PHOTON_SH, k = e.toFixed(4);
    if (!m.has(k)) m.set(k, {e, combos: []});
    m.get(k).combos.push([a, b, n]);
  }
  const out = [];
  for (const h of m.values()){
    const byParity = new Map();                       // parity flip -> cheapest combo giving it
    for (const c of h.combos.slice().sort((x, y) => x[2] - y[2])){
      if (!byParity.has(c[2] % 2)) byParity.set(c[2] % 2, c);
    }
    out.push({e: h.e, byParity: [...byParity.entries()]});
  }
  return out.sort((a, b) => a.e - b.e);
}
let HOPS = buildHops();
function comboLabel(c){
  const f = nmLabel(lamNm), h = nmLabel(lamNm / 2), p = [];
  if (c[0]) p.push((c[0] > 1 ? c[0] + '×' : '') + f);
  if (c[1]) p.push((c[1] > 1 ? c[1] + '×' : '') + h);
  return c[2] + 'γ ' + p.join('+');
}
// `from` names the control the edit came from, so writing the value back never fights typing
function setLambda(nm, from){
  const v = +nm;
  if (!isFinite(v) || v <= 0){
    if (from === 'blur') $('#lam').value = String(+lamNm.toFixed(2));
    return;
  }
  if (from === 'num' && (v < LAM_MIN || v > LAM_MAX)) return;   // mid-edit; settle on blur
  lamNm = Math.min(LAM_MAX, Math.max(LAM_MIN, v));
  PHOTON = HC / lamNm; PHOTON_SH = 2 * PHOTON;
  HOPS = buildHops();
  if (from !== 'num') $('#lam').value = String(+lamNm.toFixed(2));
  $('#lamout').textContent = 'ħω ' + eVLabel(PHOTON) + ' eV · SH ' + nmLabel(lamNm / 2)
                           + ' nm ' + eVLabel(PHOTON_SH) + ' eV';
  draw(); renderDetail();
}
// how many real levels sit on this guide: energy within tolerance, right column, |ΔJ| <= m
function countMatches(pool, pin, de, tgtPar, m){
  const TOL = 0.15;
  let n = 0;
  for (const l of pool){
    if (l.par !== tgtPar || l.id === pin.id) continue;
    if (Math.abs((l.E - pin.E) - de) > TOL) continue;
    const a = pin.Js.filter(v => v != null), b = l.Js.filter(v => v != null);
    if (a.length && b.length && !a.some(x => b.some(z => Math.abs(z - x) <= m))) continue;
    n++;
  }
  return n;
}

const $ = s => document.querySelector(s);
const cv = $('#cv'), wrap = $('#wrap'), tip = $('#tip');
const ctx = cv.getContext('2d');

let cur = null, frame = 'ion', filt = 'all', fine = false;
let view = {lo: 0, hi: 1};           // energy window, in the current frame
let pinned = null, hovered = null;
let hits = [];                       // screen geometry of what was drawn, for hit-testing
let sortKey = 'E', sortDir = 1;

const css = k => getComputedStyle(document.documentElement).getPropertyValue(k).trim();

// ---------- helpers ----------
function fmtTau(t){
  if (t == null) return '—';
  const u = [[1e-12,'ps'],[1e-9,'ns'],[1e-6,'µs'],[1e-3,'ms'],[1,'s']];
  for (const [s,n] of u) if (t < s*1000){
    const v = t/s;
    return (v>=100? v.toFixed(0) : v>=10? v.toFixed(1) : v.toFixed(2)) + ' ' + n;
  }
  return t.toExponential(2) + ' s';
}
function jList(js){
  const u = [...new Set((js||[]).filter(v => v != null))].sort((a,b) => a-b);
  return u.map(j => Number.isInteger(j) ? String(j) : (Math.round(j*2)+'/2')).join(',');
}
function termHTML(term, js){
  const t = (term||'').replace('*','');
  const odd = (term||'').includes('*');
  const m = t.match(/^\s*(?:[a-z]\s*)?(\d+)(\[[\d\/]+\]|[A-Z])\s*$/);
  const jj = jList(js);
  if (!m) return esc(t) + (jj? '<sub>'+jj+'</sub>' : '');
  return '<sup>'+m[1]+'</sup>'+esc(m[2])+(odd?'<sup>°</sup>':'')+(jj?'<sub>'+jj+'</sub>':'');
}
function termText(term, js){
  const t = (term||'').replace('*','°');
  const jj = jList(js);
  return t + (jj ? '(' + jj + ')' : '');
}
function tauLabelFull(l){          // detail panel: never hide a real spread
  if (l.tauLo == null) return l.meta ? 'metastable' : '—';
  if (l.tauHi / l.tauLo < 1.005) return fmtTau(l.tauLo);
  const a = fmtTau(l.tauLo), b = fmtTau(l.tauHi);
  return a.split(' ')[1] === b.split(' ')[1] ? a.split(' ')[0] + '–' + b : a + ' – ' + b;
}
function tauLabel(l){             // canvas + table: a multiplet within 8% reads as one number
  if (l.tauLo == null) return l.meta ? 'metastable' : '—';
  if (l.tauHi / l.tauLo < 1.08) return fmtTau(l.tau != null ? l.tau : (l.tauLo + l.tauHi) / 2);
  const a = fmtTau(l.tauLo), b = fmtTau(l.tauHi);
  return a.split(' ')[1] === b.split(' ')[1] ? a.split(' ')[0] + '–' + b : a + ' – ' + b;
}
function esc(s){ return String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c])); }
function inWin(t){ return t != null && t > WIN_LO && t < WIN_HI; }
function offset(rec){ return (frame === 'abs' && rec.off != null) ? rec.off : 0; }

let lvCache = {key: null, val: null};
function levels(rec){
  const key = rec.sp + '|' + frame + '|' + (fine ? 'r' : 'c');
  if (lvCache.key === key) return lvCache.val;
  const off = offset(rec);
  const raw = rec.lv.map((r,i) => ({
    id: 'r' + i, ri: i, ch: (rec.abs[i] || []).map(c => ['r' + c[0], c[1]]),
    conf: rec.confs[r[0]], term: rec.terms[r[1]], Js: [r[2]],
    E0: r[3], E: r[3] + off, tau: r[4], tauLo: r[4], tauHi: r[4], decay: r[5],
    par: r[6], above: r[7], meta: r[8], win: inWin(r[4]), n: 1, spread: 0
  }));
  let val = raw;
  if (!fine){
    const g = new Map();
    for (const l of raw){
      const k = l.conf + '|' + l.term;
      if (!g.has(k)) g.set(k, []);
      g.get(k).push(l);
    }
    val = [];
    const toGroup = new Map();        // raw index -> collapsed id, for remapping channels
    let c = 0;
    for (const arr of g.values()){
      arr.sort((a,b) => a.E - b.E);
      // one (config, term) label can be reused far up the ladder; that is not one multiplet
      const chunks = []; let ch = [arr[0]];
      for (let i = 1; i < arr.length; i++){
        if (arr[i].E - ch[ch.length-1].E < 0.3) ch.push(arr[i]);
        else { chunks.push(ch); ch = [arr[i]]; }
      }
      chunks.push(ch);
      for (const m of chunks){
        const gid = 'c' + c;
        for (const x of m) toGroup.set(x.ri, gid);
        const taus = m.map(x => x.tau).filter(t => t != null);
        const decs = m.map(x => x.decay).filter(t => t != null);
        val.push({
          id: 'c' + (c++), conf: m[0].conf, term: m[0].term,
          Js: m.map(x => x.Js[0]),
          E0: m.reduce((s,x) => s + x.E0, 0) / m.length,
          E:  m.reduce((s,x) => s + x.E,  0) / m.length,
          tau: taus.length ? taus.reduce((a,b) => a+b, 0) / taus.length : null,
          tauLo: taus.length ? Math.min(...taus) : null,
          tauHi: taus.length ? Math.max(...taus) : null,
          decay: decs.length ? Math.max(...decs) : null,
          par: m[0].par, above: m[0].above,
          meta: m.some(x => x.meta), win: m.some(x => x.win),
          n: m.length, spread: m[m.length-1].E - m[0].E,
          rawCh: m.flatMap(x => x.ch)
        });
      }
    }
    // a multiplet's channels are its members', remapped onto collapsed targets, strongest kept
    for (const l of val){
      const best = new Map();
      for (const [rid, f] of l.rawCh){
        const tgt = toGroup.get(+rid.slice(1));
        if (tgt && tgt !== l.id && (!best.has(tgt) || best.get(tgt) < f)) best.set(tgt, f);
      }
      l.ch = [...best.entries()].sort((a,b) => b[1] - a[1]).slice(0, 6);
      delete l.rawCh;
    }
    val.sort((a,b) => a.E - b.E);
  }
  lvCache = {key, val};
  return val;
}
function passes(l){
  if (filt === 'tau') return l.tau != null || l.meta;
  if (filt === 'win') return l.win;
  return true;
}

// ---------- species rail ----------
function buildRail(){
  const byEl = new Map();
  for (const sp of Object.keys(DATA)){
    const r = DATA[sp];
    if (!byEl.has(r.el)) byEl.set(r.el, []);
    byEl.get(r.el).push(r);
  }
  const order = [...byEl.keys()].sort((a,b) => (ELEM[a]?.[0] ?? 999) - (ELEM[b]?.[0] ?? 999));
  let html = '';
  for (const el of order){
    const list = byEl.get(el);
    list.sort((a,b) => a.q - b.q);
    const z = ELEM[el]?.[0], cat = ELEM[el]?.[1];
    html += '<div class="elgroup"><div class="elname"><b>' + esc(el) + '</b>'
          + (z && z < 999 ? '<span class="z">' + z + '</span>' : '')
          + (cat ? '<span class="cat">' + esc(cat) + '</span>' : '') + '</div>';
    for (const r of list){
      const n = r.lv.filter(x => inWin(x[4])).length;
      html += '<button class="spbtn" type="button" data-sp="' + esc(r.sp) + '">'
            + '<span>' + esc(r.sp) + '</span>'
            + '<span class="n' + (n ? ' hit' : '') + '">' + (n || '·') + '</span></button>';
    }
    html += '</div>';
  }
  $('#splist').innerHTML = html;
  $('#splist').addEventListener('click', e => {
    const b = e.target.closest('.spbtn');
    if (b) select(b.dataset.sp);
  });
}

function select(sp){
  cur = DATA[sp];
  pinned = null; hovered = null;
  document.querySelectorAll('.spbtn').forEach(b =>
    b.setAttribute('aria-current', b.dataset.sp === sp ? 'true' : 'false'));
  resetView();
}

function resetView(){
  if (!cur) return;
  let ls = levels(cur).filter(passes);
  if (!ls.length) ls = levels(cur);
  const off = offset(cur);
  let lo = 0, hi = 1;
  if (ls.length){
    lo = Math.min(...ls.map(l => l.E));
    hi = Math.max(...ls.map(l => l.E));
  }
  if (cur.ion != null) hi = Math.max(hi, cur.ion + off);
  if (frame === 'abs') lo = Math.min(lo, 0);
  const pad = Math.max((hi - lo) * 0.05, 0.4);
  view = {lo: lo - pad, hi: hi + pad};
  draw();
}

// ---------- drawing ----------
function layout(){
  const dpr = window.devicePixelRatio || 1;
  const w = wrap.clientWidth, h = wrap.clientHeight;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {w, h};
}
const PADL = 74, PADR = 22, PADT = 26, PADB = 18, RUG = 16;
const LABEL_FONT = '11.5px ui-monospace,SFMono-Regular,Menlo,monospace';
function labelText(l){
  return shortConf(l.conf) + ' ' + termText(l.term, l.Js) + '  ' + tauLabel(l);
}

function draw(){
  const {w, h} = layout();
  ctx.clearRect(0, 0, w, h);
  if (!cur) return;
  const off = offset(cur);
  const span = view.hi - view.lo;
  const pxPerEv = (h - PADT - PADB) / span;
  const y = e => PADT + (view.hi - e) * pxPerEv;
  const mid = PADL + RUG + (w - PADL - RUG - PADR) / 2;
  const colW = Math.min(150, (w - PADL - RUG - PADR) / 2 - 90);
  const inkA = css('--accent'), ink = css('--ink'), ink2 = css('--ink-2'), ink3 = css('--ink-3');
  const plain = css('--plain'), known = css('--known'), rule = css('--rule-2');

  // ionization band + limits
  if (cur.ion != null){
    const yi = y(cur.ion + off);
    ctx.fillStyle = css('--band');
    ctx.fillRect(PADL, PADT, w - PADL - PADR, Math.max(0, yi - PADT));
    ctx.strokeStyle = css('--warm'); ctx.lineWidth = 1; ctx.setLineDash([7, 4]);
    ctx.beginPath(); ctx.moveTo(PADL, yi + .5); ctx.lineTo(w - PADR, yi + .5); ctx.stroke();
    ctx.setLineDash([]);
    if (yi > PADT + 12 && yi < h - PADB){
      ctx.fillStyle = css('--warm'); ctx.font = '11px ui-sans-serif,system-ui,sans-serif';
      ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
      ctx.fillText('ionization limit ' + (cur.ion + off).toFixed(3) + ' eV', w - PADR - 3, yi - 3);
    }
  }
  // the ion's own ground state, only meaningful in the absolute frame
  if (frame === 'abs' && off > 0){
    const yg = y(off);
    ctx.strokeStyle = ink3; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(PADL, yg + .5); ctx.lineTo(w - PADR, yg + .5); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = ink3; ctx.font = '11px ui-sans-serif,system-ui,sans-serif';
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText(cur.sp + ' ground — ' + off.toFixed(2) + ' eV of ionization paid first',
                 PADL + 6, yg + 4);
  }

  // axis
  ctx.strokeStyle = css('--rule'); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(PADL + .5, PADT); ctx.lineTo(PADL + .5, h - PADB); ctx.stroke();
  const step = niceStep(span, (h - PADT - PADB) / 46);
  ctx.font = '11px ui-monospace,SFMono-Regular,Menlo,monospace';
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillStyle = ink3;
  for (let e = Math.ceil(view.lo / step) * step; e <= view.hi; e += step){
    const yy = y(e);
    ctx.strokeStyle = rule;
    ctx.beginPath(); ctx.moveTo(PADL, yy + .5); ctx.lineTo(w - PADR, yy + .5); ctx.stroke();
    ctx.fillText(fmtAxis(e, step), PADL - 7, yy);
  }

  const all = levels(cur);
  const shown = all.filter(passes);

  // density rug: nothing is ever hidden, even when the ladder is thinned
  ctx.strokeStyle = plain; ctx.lineWidth = 1; ctx.globalAlpha = pinned ? .3 : .75;
  for (const l of shown){
    const yy = y(l.E);
    if (yy < PADT || yy > h - PADB) continue;
    ctx.beginPath(); ctx.moveTo(PADL + 4, yy + .5); ctx.lineTo(PADL + RUG - 4, yy + .5); ctx.stroke();
  }
  ctx.globalAlpha = 1;

  // level of detail: plain levels only once the ladder is open enough to separate them
  const openEnough = pxPerEv > 9;
  const vis = shown.filter(l => {
    const yy = y(l.E);
    if (yy < PADT - 2 || yy > h - PADB + 2) return false;
    return l.win || l.tau != null || l.meta || openEnough;
  });


  hits = [];
  const labelRows = [];
  const FADE = 0.22;                 // once something is pinned, everything else recedes
  for (const l of vis){
    const yy = y(l.E);
    const sgn = l.par ? 1 : -1;
    const x0 = mid + sgn * 14, x1 = mid + sgn * (14 + colW);
    const isSel = pinned && pinned.id === l.id, isHov = hovered && hovered.id === l.id;
    ctx.globalAlpha = (pinned && !isSel && !isHov) ? FADE : 1;
    ctx.strokeStyle = l.win ? inkA : (l.tau != null || l.meta ? known : plain);
    ctx.lineWidth = l.win ? 2.4 : (l.tau != null ? 1.6 : 1);
    if (isSel){                      // the pinned level: full-width rule plus a heavy bar
      ctx.strokeStyle = ink; ctx.lineWidth = 0.8; ctx.globalAlpha = 0.5;
      ctx.beginPath(); ctx.moveTo(PADL + RUG, yy + .5); ctx.lineTo(w - PADR, yy + .5); ctx.stroke();
      ctx.globalAlpha = 1; ctx.strokeStyle = ink; ctx.lineWidth = 3.4;
    } else if (isHov){ ctx.strokeStyle = ink; ctx.lineWidth = 2.6; }
    ctx.beginPath(); ctx.moveTo(Math.min(x0,x1), yy + .5); ctx.lineTo(Math.max(x0,x1), yy + .5); ctx.stroke();
    ctx.globalAlpha = 1;
    hits.push({l, y: yy, par: l.par, mid});
    if (l.win || isSel || isHov) labelRows.push({l, yy, sgn, x: x1, sel: isSel || isHov});
  }

  // labels for the levels that matter, thinned so they never stack
  ctx.font = LABEL_FONT;
  ctx.textBaseline = 'middle';
  const taken = {'-1': [], '1': []};
  labelRows.sort((a, b) => a.yy - b.yy);
  for (const r of labelRows){
    const lane = taken[String(r.sgn)];
    if (lane.some(v => Math.abs(v - r.yy) < 13)) continue;
    lane.push(r.yy);
    ctx.textAlign = r.sgn > 0 ? 'left' : 'right';
    ctx.globalAlpha = (pinned && !r.sel) ? FADE : 1;
    ctx.fillStyle = r.sel ? ink : (r.l.win ? inkA : ink2);
    ctx.fillText(labelText(r.l), r.x + r.sgn * 7, r.yy);
    ctx.globalAlpha = 1;
  }

  // parity headings
  ctx.font = '11px ui-sans-serif,system-ui,sans-serif'; ctx.textBaseline = 'top';
  ctx.fillStyle = ink3;
  ctx.textAlign = 'right'; ctx.fillText('even parity', mid - 14, 6);
  ctx.textAlign = 'left';  ctx.fillText('odd parity',  mid + 14, 6);

  // Hop guides from the pinned level. Each photon flips parity, so an m-photon hop lands on
  // the same column when m is even and the other column when m is odd. A guide is therefore
  // drawn ONLY over the column it can legally land in -- a full-width line claimed both.
  if (pinned){
    ctx.font = '10.5px ui-sans-serif,system-ui,sans-serif';
    for (const hop of HOPS){
      const yy = y(pinned.E + hop.e);
      if (yy < PADT || yy > h - PADB) continue;
      for (const [flip, combo] of hop.byParity){
        const tgt = pinned.par ^ flip;                 // column this hop can reach
        const sgn = tgt ? 1 : -1;
        const xa = mid + sgn * 12, xb = mid + sgn * (18 + colW);
        const hit = countMatches(shown, pinned, hop.e, tgt, combo[2]);
        ctx.strokeStyle = inkA; ctx.globalAlpha = hit ? .85 : .32;
        ctx.lineWidth = hit ? 1.4 : 1;
        ctx.setLineDash(hit ? [5, 3] : [2, 4]);
        ctx.beginPath(); ctx.moveTo(Math.min(xa,xb), yy + .5); ctx.lineTo(Math.max(xa,xb), yy + .5); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = inkA; ctx.globalAlpha = hit ? 1 : .5;
        ctx.textAlign = sgn > 0 ? 'left' : 'right';
        ctx.textBaseline = 'bottom';
        ctx.fillText('+' + hop.e.toFixed(3) + ' eV · ' + comboLabel(combo)
                     + (hit ? '  → ' + hit + (hit > 1 ? ' levels' : ' level') : ''),
                     xb + sgn * 6, yy - 2);
        ctx.globalAlpha = 1;
      }
    }
  }

  $('#zoomnote').textContent = pxPerEv.toFixed(1) + ' px/eV · '
    + vis.length + ' of ' + shown.length + ' levels drawn'
    + (openEnough ? '' : ' · zoom in for the levels without lifetime data');
  scheduleTable(shown.filter(l => { const yy = y(l.E); return yy >= PADT && yy <= h - PADB; }));
}

function shortConf(c){
  return c.replace(/\([^)]*\)/g, '').replace(/\.+/g, '.').replace(/^\.|\.$/g, '');
}
function niceStep(span, targetTicks){
  const raw = span / Math.max(targetTicks, 2);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 2.5, 5, 10]) if (raw <= m * p) return m * p;
  return 10 * p;
}
function fmtAxis(e, step){
  const d = step >= 1 ? 0 : step >= 0.1 ? 1 : step >= 0.01 ? 2 : 3;
  return e.toFixed(d);
}

// ---------- detail panel ----------
function renderDetail(){
  const box = $('#detail');
  if (!cur){ box.innerHTML = ''; return; }
  const off = offset(cur);
  const l = hovered || pinned;
  let html = '<div class="card"><h3 class="name">' + esc(cur.sp) + '</h3><dl class="kv">'
    + row('Levels', cur.lv.length)
    + row('With a lifetime', cur.lv.filter(r => r[4] != null).length)
    + row('In 1–1000 ps', cur.lv.filter(r => inWin(r[4])).length)
    + row('Ionization limit', cur.ion != null ? (cur.ion).toFixed(4) + ' eV' : '—')
    + row('Ionization already paid', cur.q === 0 ? 'none — neutral'
        : cur.off != null ? cur.off.toFixed(2) + ' eV' : 'unknown')
    + row('ASD lines with A', cur.nAki + ' / ' + cur.nLines)
    + (cur.vap ? row(cur.el + ' vapour', cur.vap) : '')
    + '</dl></div>';

  if (l){
    html += '<div class="card"><h3>' + (pinned && (!hovered || hovered.id === pinned.id) ? 'Pinned level' : 'Level') + '</h3>'
      + '<div class="mono" style="font-size:13px;margin-bottom:7px">'
      + esc(shortConf(l.conf)) + ' ' + termHTML(l.term, l.Js) + '</div><dl class="kv">'
      + row('Energy, ion frame', l.E0.toFixed(4) + ' eV')
      + (cur.off ? row('Energy, neutral frame', (l.E0 + cur.off).toFixed(3) + ' eV') : '')
      + row('Lifetime', tauLabelFull(l))
      + (l.n > 1 ? row('Fine structure', l.n + ' levels · J = ' + jList(l.Js)
            + ' · spans ' + (l.spread*1000).toFixed(0) + ' meV') : '')
      + row('Strongest decay', l.decay != null ? l.decay.toFixed(3) + ' eV' : '—')
      + row('Parity', l.par ? 'odd' : 'even')
      + row('Above ionization', l.above ? 'yes — autoionizing' : 'no')
      + '</dl>'
      + (l.win ? '<div class="mark">inside 1–1000 ps</div>' : '')
      + '</div>';
    if (l.ch && l.ch.length){
      const map = new Map(levels(cur).map(x => [x.id, x]));
      html += '<div class="card"><h3>Linear absorption from here</h3>'
            + '<table class="abs"><thead><tr><th>upper level</th><th>&Delta;E</th>'
            + '<th>&lambda;</th><th>f</th></tr></thead><tbody>';
      for (const [id, f] of l.ch){
        const t = map.get(id);
        if (!t) continue;
        const dE = t.E - l.E;
        html += '<tr><td class="mono">' + esc(shortConf(t.conf)) + ' ' + termHTML(t.term, t.Js)
              + '</td><td class="mono">' + dE.toFixed(3) + '</td>'
              + '<td class="mono">' + (1239.841984 / dE).toFixed(0) + '</td>'
              + '<td class="mono">' + (f >= 0.1 ? f.toFixed(2) : f.toPrecision(2)) + '</td></tr>';
      }
      html += '</tbody></table><p class="note">One photon, electric dipole. '
            + '<i>f</i> is the oscillator strength; &int;&sigma;d&nu; = 2.654&times;10'
            + '<sup>&minus;2</sup> <i>f</i> cm<sup>2</sup>Hz. &Delta;E in eV, &lambda; in nm.</p></div>';
    }
  }

  if (pinned){
    const d = hovered && hovered.id !== pinned.id ? (hovered.E - pinned.E) : null;
    html += '<div class="card"><h3>Measuring from the pin</h3>'
      + '<div class="mono" style="font-size:12.5px;margin-bottom:6px">'
      + esc(shortConf(pinned.conf)) + ' ' + termHTML(pinned.term, pinned.Js) + '</div>';
    if (d != null){
      const n = d / PHOTON, nsh = d / PHOTON_SH;
      html += '<dl class="kv">'
        + row('Gap to hovered', (d >= 0 ? '+' : '') + d.toFixed(3) + ' eV')
        + row('at ' + nmLabel(lamNm) + ' nm', n.toFixed(2) + ' × ' + eVLabel(PHOTON) + ' eV')
        + row('at ' + nmLabel(lamNm / 2) + ' nm', nsh.toFixed(2) + ' × ' + eVLabel(PHOTON_SH) + ' eV')
        + row('Parity change', (pinned.par !== hovered.par) ? 'yes — needs odd m' : 'no — needs even m')
        + '</dl>';
    } else {
      html += '<p class="empty">Hover another level to read the gap.</p>';
    }
    html += '<div style="margin-top:9px"><button class="btn" id="unpin" type="button">Clear pin</button></div></div>';
  }
  box.innerHTML = html;
  const u = $('#unpin');
  if (u) u.addEventListener('click', () => { pinned = null; draw(); renderDetail(); });
}
function row(k, v){ return '<dt>' + esc(k) + '</dt><dd class="mono">' + esc(v) + '</dd>'; }

// ---------- table ----------
let tabTimer = null;
function scheduleTable(rows){
  if (tabTimer) clearTimeout(tabTimer);
  tabTimer = setTimeout(() => renderTable(rows), 80);
}
function renderTable(rows){
  const key = l => sortKey === 'J' ? (l.Js[0]) : l[sortKey];
  const s = rows.slice().sort((a, b) => {
    const va = key(a), vb = key(b);
    if (va == null) return 1;
    if (vb == null) return -1;
    return (va > vb ? 1 : va < vb ? -1 : 0) * sortDir;
  });
  const cap = 400;
  $('#tbody').innerHTML = s.slice(0, cap).map(l =>
    '<tr class="' + (l.win ? 'win ' : '') + (pinned && pinned.id === l.id ? 'sel' : '') + '" data-id="' + l.id + '">'
    + '<td class="mono">' + esc(shortConf(l.conf)) + '</td>'
    + '<td class="mono">' + termHTML(l.term, null) + '</td>'
    + '<td class="mono">' + (jList(l.Js) || '—') + (l.n > 1 ? ' <span style="opacity:.55">×' + l.n + '</span>' : '') + '</td>'
    + '<td class="mono">' + l.E.toFixed(4) + '</td>'
    + '<td class="mono">' + tauLabel(l) + '</td>'
    + '<td class="mono">' + (l.decay != null ? l.decay.toFixed(3) : '—') + '</td>'
    + '<td>' + (l.par ? 'odd' : 'even') + '</td></tr>').join('');
  $('#tabnote').textContent = rows.length + ' in view'
    + (rows.length > cap ? ' · first ' + cap + ' listed' : '')
    + ' · ' + rows.filter(l => l.win).length + ' inside 1–1000 ps';
}

// ---------- interaction ----------
function pick(ev){
  const r = cv.getBoundingClientRect();
  const mx = ev.clientX - r.left, my = ev.clientY - r.top;
  const side = hits.length && mx >= hits[0].mid ? 1 : 0;
  let best = null, bd = 7, bestAny = null, ba = 7;
  for (const hgt of hits){
    const d = Math.abs(hgt.y - my);
    if (d < ba){ ba = d; bestAny = hgt.l; }
    if (hgt.par === side && d < bd){ bd = d; best = hgt.l; }
  }
  return {l: best || bestAny, mx, my};
}
cv.addEventListener('mousemove', ev => {
  const {l, mx, my} = pick(ev);
  const changed = (l && (!hovered || hovered.id !== l.id)) || (!l && hovered);
  hovered = l;
  if (l){
    tip.classList.add('on');
    tip.innerHTML = '<b class="mono">' + esc(shortConf(l.conf)) + ' ' + termHTML(l.term, l.Js) + '</b>'
      + '<div class="row"><span>energy</span><span class="mono">' + l.E.toFixed(4) + ' eV</span></div>'
      + '<div class="row"><span>lifetime</span><span class="mono">' + tauLabel(l) + '</span></div>'
      + (l.n > 1 ? '<div class="row"><span>fine structure</span><span class="mono">' + l.n + ' levels, ' + (l.spread*1000).toFixed(0) + ' meV</span></div>' : '')
      + (l.decay != null ? '<div class="row"><span>strongest decay</span><span class="mono">' + l.decay.toFixed(2) + ' eV</span></div>' : '')
      + (pinned && pinned.id !== l.id ? '<div class="row"><span>from pin</span><span class="mono">'
          + ((l.E - pinned.E) >= 0 ? '+' : '') + (l.E - pinned.E).toFixed(3) + ' eV</span></div>' : '');
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    tip.style.left = Math.min(mx + 16, wrap.clientWidth - tw - 8) + 'px';
    tip.style.top = Math.max(6, Math.min(my - th / 2, wrap.clientHeight - th - 6)) + 'px';
  } else tip.classList.remove('on');
  if (changed){ draw(); renderDetail(); }
});
cv.addEventListener('mouseleave', () => { hovered = null; tip.classList.remove('on'); draw(); renderDetail(); });
cv.addEventListener('click', ev => {
  const {l} = pick(ev);
  if (l) pinned = (pinned && pinned.id === l.id) ? null : l;
  draw(); renderDetail();
});
cv.addEventListener('wheel', ev => {
  ev.preventDefault();
  const r = cv.getBoundingClientRect();
  const f = (ev.clientY - r.top - PADT) / Math.max(1, r.height - PADT - PADB);
  const at = view.hi - f * (view.hi - view.lo);
  const k = Math.exp(ev.deltaY * 0.0016);
  const lo = at - (at - view.lo) * k, hi = at + (view.hi - at) * k;
  if (hi - lo > 1e-4) { view = {lo, hi}; draw(); }
}, {passive: false});

let drag = null;
cv.addEventListener('pointerdown', ev => { drag = {y: ev.clientY, v: {...view}, moved: false}; cv.setPointerCapture(ev.pointerId); });
cv.addEventListener('pointermove', ev => {
  if (!drag) return;
  const dy = ev.clientY - drag.y;
  if (Math.abs(dy) > 2) drag.moved = true;
  const perPx = (drag.v.hi - drag.v.lo) / Math.max(1, cv.clientHeight - PADT - PADB);
  view = {lo: drag.v.lo + dy * perPx, hi: drag.v.hi + dy * perPx};
  draw();
});
cv.addEventListener('pointerup', ev => { drag = null; cv.releasePointerCapture(ev.pointerId); });

document.addEventListener('keydown', ev => {
  if (ev.key === 'Escape'){ pinned = null; draw(); renderDetail(); }
});
$('#tbody').addEventListener('mouseover', e => {
  const tr = e.target.closest('tr');
  if (!tr || !cur) return;
  const l = levels(cur).find(x => x.id === tr.dataset.id);
  if (l){ hovered = l; draw(); renderDetail(); }
});
$('#tbody').addEventListener('click', e => {
  const tr = e.target.closest('tr');
  if (!tr || !cur) return;
  const l = levels(cur).find(x => x.id === tr.dataset.id);
  if (l){ pinned = (pinned && pinned.id === l.id) ? null : l; draw(); renderDetail(); }
});
document.querySelectorAll('#view-levels thead th').forEach(th => th.addEventListener('click', () => {
  const k = th.dataset.k;
  if (sortKey === k) sortDir *= -1; else { sortKey = k; sortDir = 1; }
  draw();
}));
$('#frame').addEventListener('change', e => { frame = e.target.value; resetView(); renderDetail(); });
$('#fine').addEventListener('change', e => {
  fine = e.target.value === '1';
  pinned = null; hovered = null;      // ids differ between the two modes
  draw(); renderDetail();
});
$('#filt').addEventListener('change', e => { filt = e.target.value; resetView(); });
$('#lam').addEventListener('focus', e => e.target.select());
$('#lam').addEventListener('input', e => setLambda(e.target.value, 'num'));
$('#lam').addEventListener('blur', e => setLambda(e.target.value, 'blur'));
$('#reset').addEventListener('click', () => { pinned = null; resetView(); renderDetail(); });
window.addEventListener('resize', draw);
matchMedia('(prefers-color-scheme:dark)').addEventListener('change', draw);

$('#tab-levels').addEventListener('click', () => showTab('levels'));
$('#tab-ladders').addEventListener('click', () => showTab('ladders'));

buildRail();
select(Object.keys(DATA).find(s => DATA[s].lv.some(r => inWin(r[4]))) || Object.keys(DATA)[0]);
setLambda(lamNm, 'init');
initLadders();
"""


TITLE = "Atomic level explorer — burst ladder screening"


def build() -> tuple[str, str]:
    data = collect()
    newest = max(os.path.getmtime(f) for f in glob.glob(str(OUT_ROOT / "*" / "asd_levels.csv")))
    stamp = datetime.date.fromtimestamp(newest).strftime("%d %B %Y")
    mat, lad = ladders_tab.collect()
    js = ((JS + ladders_tab.JS)
            .replace("__MAT__", json.dumps(mat, separators=(",", ":")))
            .replace("__LAD__", json.dumps(lad, separators=(",", ":")))
            .replace("__ELEM__", json.dumps(element_table(data), separators=(",", ":")))
            .replace("__DATA__", json.dumps(data, separators=(",", ":")))
            .replace("__WIN_LO__", repr(WIN_LO)).replace("__WIN_HI__", repr(WIN_HI))
            .replace("__LAMBDA__", repr(LAMBDA_NM)))
    body = BODY.replace("__DATE__", stamp).replace("__LADDERS_BODY__", ladders_tab.BODY)
    css = CSS + ladders_tab.CSS
    fragment = (f"<title>{TITLE}</title>\n<style>{css}</style>\n{body}\n<script>{js}</script>\n")
    standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                  '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
                  f"<title>{TITLE}</title>\n<style>{css}</style>\n</head>\n<body>\n"
                  f"{body}\n<script>{js}</script>\n</body>\n</html>\n")
    return standalone, fragment


def main() -> int:
    standalone, fragment = build()
    out = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "index.html"
    out.write_text(standalone, encoding="utf-8")
    frag = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else ROOT / "explorer_fragment.html"
    frag.write_text(fragment, encoding="utf-8")
    print(f"  {out}  ({len(standalone)/1024:.0f} KB)")
    print(f"  {frag}  ({len(fragment)/1024:.0f} KB, body fragment for publishing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
