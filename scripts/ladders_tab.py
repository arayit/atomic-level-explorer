#!/usr/bin/env python3
"""The Ladders tab: the burst-ladder screen, rendered into the level explorer.

Reads what scripts/ladder.py and scripts/materials.py wrote to analysis/ and emits three
fragments -- CSS, markup, JS -- that scripts/explorer.py splices into the page.

Only a slice of the ladder set is carried into the page. The scan produces hundreds of
thousands of chains, most of them fine-structure neighbours of each other; the page keeps each
species' best few under two rankings at once -- the overall one, and the one restricted to
chains that end on a BOUND level, because a top state that autoionizes ends the ladder in an
ion rather than a photon.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANA = ROOT / "analysis"

PER_SPECIES = 26          # best overall
PER_SPECIES_BOUND = 12    # best that end on a bound level, kept even if outranked overall
PER_SPECIES_EMIT = 14     # best that end on a bound level which is itself fast, i.e. EMITS

EV_RANK = {"measured": 0, "autoionizing": 1, "none": 2}
DET_RANK = {"A": 0, "B": 1, "C": 2}


def _key(r):
    return (EV_RANK[r["subns_evidence"]], 1 if int(r["unverified_parked"] or 0) else 0,
            DET_RANK[r["detuning_tier"]],
            0 if len(set(r["colours_nm"].split("+"))) == 1 else 1,
            -int(r["order"]), float(r["worst_detuning_meV"]))


def _row(r):
    """One ladder, as a compact array. Field order is mirrored by unpack() in the page JS."""
    names = r["path"].split(" -> ")
    return [
        r["start"], float(r["start_eV"]), names,
        [float(x) for x in r["path_eV"].split(" -> ")],
        [None if t == "?" else float(f"{float(t):.3g}") for t in r["path_tau_s"].split(" -> ")],
        [int(x) for x in r["photons_per_hop"].split("+")],
        [float(x) for x in r["colours_nm"].split("+")],
        round(float(r["worst_detuning_meV"]), 1), int(r["order"]),
        r["subns_evidence"], r["fast_rung_at"],
        float(r["fast_rung_tau_s"]) if r["fast_rung_tau_s"] else None,
        1 if r["final_above_ionization"] == "yes" else 0,
        float(r["final_max_decay_eV"]) if r["final_max_decay_eV"] else None,
        float(r["final_eV_neutral_frame"]) if r["final_eV_neutral_frame"] else None,
        int(r["unverified_parked"] or 0),
        r["start_kind"],
    ]


def collect():
    mat = list(csv.DictReader((ANA / "materials.csv").open(encoding="utf-8")))
    by = {}
    for r in csv.DictReader((ANA / "ladders.csv").open(encoding="utf-8")):
        by.setdefault(r["spectrum"], []).append(r)

    lad = {}
    for sp, rs in by.items():
        rs.sort(key=_key)
        pick, seen = [], set()
        bound = [r for r in rs if r["final_above_ionization"] == "no"][:PER_SPECIES_BOUND]
        # The chains that answer "does this end in a photon?": bound top state, its own
        # lifetime inside the picosecond window, and a published decay energy to quote a
        # wavelength from. Ranked by photon energy, so the shortest wavelength leads.
        emit = [r for r in rs if r["final_above_ionization"] == "no"
                and r["final_max_decay_eV"] and r["fast_rung_at"] in ("final", "both")]
        emit.sort(key=lambda r: (-float(r["final_max_decay_eV"]),
                                 float(r["worst_detuning_meV"])))
        for r in rs[:PER_SPECIES] + bound + emit[:PER_SPECIES_EMIT]:
            k = (r["start"], r["path"], r["photons_per_hop"], r["colours_nm"])
            if k in seen:
                continue
            seen.add(k)
            pick.append(_row(r))
        lad[sp] = pick

    keep = ("tier", "spectrum", "element", "Z", "charge", "prep_note", "prep_cost", "evidence",
            "fast_rung_at", "rung_lifetimes", "closure", "best_detuning_meV", "hop_amplitude",
            "relative_amplitude", "order", "hops", "max_order", "ground_start",
            "single_colour", "reach_eV_neutral_frame", "n_ladders",
            "accidental_matches_per_hop", "truncated")
    return [{k: r.get(k, "") for k in keep} for r in mat], lad


CSS = r"""
/* ---------- tabs ---------- */
.tabs{display:flex; gap:0; padding:0 18px; border-bottom:1px solid var(--rule); background:var(--paper)}
.tabs button{
  font:inherit; font-size:13px; color:var(--ink-3); background:none; cursor:pointer;
  border:1px solid transparent; border-bottom:0; padding:5px 14px; margin-bottom:-1px;
}
.tabs button[aria-selected="true"]{
  color:var(--ink); background:var(--panel-2);
  border-color:var(--rule); border-bottom:1px solid var(--panel-2);
}
.tabs button:focus-visible{outline:2px solid var(--link); outline-offset:-2px}

/* ---------- ladders view ---------- */
.lad{padding:0 0 40px}
.lad .intro{padding:12px 18px 8px; font-size:13px; color:var(--ink-2); max-width:100ch}
.lad .intro p{margin:0 0 6px}
.lad h2{margin:22px 0 6px; padding:0 18px; font-size:15px; font-weight:normal}
.lad h2 .muted{font-size:12px; color:var(--ink-3)}
.lad .scroll{overflow-x:auto; padding:0 18px}
.lad table{border-collapse:collapse; font-size:12.5px; width:100%}
.lad th,.lad td{border-bottom:1px solid var(--rule-2); padding:4px 9px 4px 0; text-align:left;
  vertical-align:top; white-space:nowrap}
.lad th{border-bottom:1px solid var(--rule); font-weight:normal; color:var(--ink-3);
  position:sticky; top:0; background:var(--paper)}
.lad td.num{text-align:right; font-variant-numeric:tabular-nums}
.lad tr.t1 td:first-child{border-left:3px solid var(--accent); padding-left:6px}
.lad .pill{font-size:11px; color:var(--ink-3); border:1px solid var(--rule-2); padding:0 4px}
.lad .emit{color:var(--accent-ink)}
.lad .chain{font-family:ui-monospace,Menlo,Consolas,monospace; font-size:12px}
.lad .terms{color:var(--ink-3); font-size:11.5px; white-space:nowrap}
.lad .none{padding:14px 18px; color:var(--ink-3)}
.lad tbody tr:hover{background:var(--panel-2)}
"""

BODY = r"""
<section class="lad" id="view-ladders" hidden>
  <div class="intro">
    <p>Every chain below starts on a level that holds population, climbs in hops of 2 or 3
      photons of one colour, and lands on a real level. Each hop matches a measured level gap
      to within the stated detuning, obeys the parity rule for that photon number, and obeys
      |&Delta;J| &le; m. Every rung except the last outlives the intra-burst pulse spacing of
      20&nbsp;ps, so the next pulse arrives at an occupied level.</p>
    <p><b>Top state</b> is the axis to read first. A chain ending above the ionization limit
      lands on an autoionizing resonance: intrinsically fast, but it decays by throwing out an
      electron, so the ladder ends in an ion rather than a photon. A chain ending on a bound
      level with a picosecond lifetime ends in EUV light instead, and the emitted wavelength is
      given. Those exist only in ions &mdash; ionizing once roughly doubles the ionization
      limit while the level structure keeps its shape, so states that autoionize in the neutral
      sit safely below the limit one step up the isoelectronic sequence.</p>
  </div>

  <form class="controls" onsubmit="return false">
    <label>Top state
      <select id="lTop">
        <option value="any" selected>either</option>
        <option value="bound">bound &mdash; emits a photon</option>
        <option value="ai">autoionizing &mdash; ends in an ion</option>
      </select>
    </label>
    <label>Closure
      <select id="lDet">
        <option value="150" selected>within 150 meV</option>
        <option value="60">within 60 meV</option>
        <option value="25">within 25 meV &mdash; driver tuning alone</option>
      </select>
    </label>
    <label>Colours
      <select id="lCol"><option value="any" selected>any</option>
        <option value="one">one colour only</option></select>
    </label>
    <label>Rung lifetimes
      <select id="lVer"><option value="any" selected>any</option>
        <option value="all">all published</option></select>
    </label>
    <label>Start
      <select id="lStart"><option value="any" selected>either</option>
        <option value="ground">ground state</option>
        <option value="metastable">metastable</option></select>
    </label>
    <label>Species <select id="lSp"></select></label>
    <span class="derived mono" id="lCount"></span>
  </form>

  <h2>Materials <span class="muted" id="lMatNote"></span></h2>
  <div class="scroll"><table>
    <thead><tr>
      <th>Tier</th><th>Species</th><th>Z</th><th>Preparation</th><th>Fast rung</th>
      <th>Rung lifetimes</th><th>Closure</th><th>Hop amplitude</th><th>Order</th><th>Hops</th>
      <th>One colour</th><th>Reach</th><th>Chance matches</th>
    </tr></thead>
    <tbody id="lMat"></tbody>
  </table></div>

  <h2>Chains <span class="muted" id="lLadNote"></span></h2>
  <div class="scroll"><table>
    <thead><tr>
      <th>Species</th><th>Order</th><th>Hops</th><th>&Delta; worst</th><th>Chain</th>
      <th>Top state</th>
    </tr></thead>
    <tbody id="lLad"></tbody>
  </table></div>
</section>
"""

JS = r"""
/* ---------------- Ladders tab ---------------- */
const MAT = __MAT__, LAD = __LAD__;

function unpack(a){
  return {start:a[0], startE:a[1], names:a[2], Es:a[3], taus:a[4], ms:a[5], nms:a[6],
          det:a[7], order:a[8], ev:a[9], fastAt:a[10], fastTau:a[11], ai:a[12],
          decay:a[13], neutral:a[14], unver:a[15], kind:a[16]};
}

/* Number(x.toPrecision(3)) drops trailing zeros only after a decimal point. Stripping them
   with a regex instead turns 250 ps into 25 ps. */
const sig3 = v => String(+v.toPrecision(3));
function tauText(t){
  if (t == null) return "–";
  if (t < 1e-9) return sig3(t * 1e12) + " ps";
  if (t < 1e-6) return sig3(t * 1e9) + " ns";
  if (t < 1e-3) return sig3(t * 1e6) + " µs";
  return sig3(t) + " s";
}

function ladderRow(sp, L){
  const chain = [L.startE.toFixed(3)];
  const arrows = [];
  for (let i = 0; i < L.ms.length; i++){
    arrows.push(L.ms[i] + "γ " + nmLabel(L.nms[i]) + " nm");
    chain.push(L.Es[i].toFixed(3));
  }
  let line = esc(chain[0]);
  for (let i = 0; i < arrows.length; i++)
    line += ' <span class="muted">─' + esc(arrows[i]) + '→</span> ' + esc(chain[i + 1]);
  const terms = esc([L.start].concat(L.names).join("  →  "));

  let top;
  if (L.ai){
    top = '<span class="muted">autoionizing</span> ' + esc(L.Es[L.Es.length - 1].toFixed(2))
        + " eV";
  } else if (L.decay && L.fastTau != null && L.fastAt !== "intermediate"){
    top = '<span class="emit">emits ' + esc(nmLabel(HC / L.decay)) + " nm</span> · τ "
        + esc(tauText(L.fastTau));
  } else {
    top = '<span class="muted">bound</span> ' + esc(L.Es[L.Es.length - 1].toFixed(2)) + " eV"
        + (L.fastTau != null ? ' · fast rung τ ' + esc(tauText(L.fastTau)) : "");
  }
  const flag = L.unver ? ' <span class="pill">' + L.unver + ' rung(s) unpublished</span>' : "";
  return '<tr><td>' + esc(sp) + (L.kind === "ground" ? ' <span class="pill">ground</span>' : "")
       + '</td><td class="num">' + L.order + '</td><td class="num">' + L.ms.length
       + '</td><td class="num">' + L.det.toFixed(1) + ' meV</td>'
       + '<td><div class="chain">' + line + '</div><div class="terms">' + terms + '</div></td>'
       + '<td>' + top + flag + '</td></tr>';
}

function matRow(m){
  const cells = [m.tier, m.spectrum, m.Z, m.prep_note,
                 m.evidence + (m.fast_rung_at ? " (" + m.fast_rung_at + ")" : ""),
                 m.rung_lifetimes, m.closure + " · " + m.best_detuning_meV + " meV",
                 m.hop_amplitude + (m.relative_amplitude ? " " + m.relative_amplitude : ""),
                 m.order, m.hops, m.single_colour,
                 m.reach_eV_neutral_frame ? m.reach_eV_neutral_frame + " eV" : "–",
                 m.accidental_matches_per_hop];
  return '<tr class="t' + m.tier + '">'
       + cells.map((c, i) => '<td' + ([2, 8, 9, 12].includes(i) ? ' class="num"' : "") + '>'
                             + esc(String(c)) + '</td>').join("") + '</tr>';
}

function ladFilters(){
  return {top: $('#lTop').value, det: +$('#lDet').value, col: $('#lCol').value,
          ver: $('#lVer').value, kind: $('#lStart').value, sp: $('#lSp').value};
}

function ladPass(L, f){
  if (f.top === "bound" && L.ai) return false;
  if (f.top === "ai" && !L.ai) return false;
  if (L.det > f.det) return false;
  if (f.col === "one" && new Set(L.nms).size !== 1) return false;
  if (f.ver === "all" && L.unver) return false;
  if (f.kind !== "any" && L.kind !== f.kind) return false;
  return true;
}

function drawLadders(){
  const f = ladFilters();
  const names = f.sp === "any" ? Object.keys(LAD) : [f.sp];
  const rows = [], live = new Set();
  for (const sp of names)
    for (const a of (LAD[sp] || [])){
      const L = unpack(a);
      if (ladPass(L, f)){ rows.push([sp, L]); live.add(sp); }
    }
  const emits = L => (!L.ai && L.decay && L.fastTau != null && L.fastAt !== "intermediate");
  rows.sort((x, y) => (emits(y[1]) - emits(x[1]))
                   || (emits(x[1]) ? (y[1].decay - x[1].decay) : 0)
                   || (x[1].ai - y[1].ai) || (y[1].order - x[1].order) || (x[1].det - y[1].det));

  const mat = MAT.filter(m => live.has(m.spectrum));
  $('#lMat').innerHTML = mat.length ? mat.map(matRow).join("")
    : '<tr><td colspan="13" class="none">No species has a chain that passes these filters.</td></tr>';
  $('#lLad').innerHTML = rows.length ? rows.map(r => ladderRow(r[0], r[1])).join("")
    : '<tr><td colspan="6" class="none">No chain passes these filters.</td></tr>';
  $('#lMatNote').textContent = mat.length + " of " + MAT.length + " species";
  $('#lLadNote').textContent = rows.length + " chains shown, best first";
  $('#lCount').textContent = rows.length + " chains · " + live.size + " species";
}

function initLadders(){
  const sel = $('#lSp');
  sel.innerHTML = '<option value="any">all</option>'
    + MAT.filter(m => LAD[m.spectrum]).map(m =>
        '<option value="' + esc(m.spectrum) + '">' + esc(m.spectrum) + "</option>").join("");
  ["lTop", "lDet", "lCol", "lVer", "lStart", "lSp"].forEach(id =>
    $(id).addEventListener("change", drawLadders));
  drawLadders();
}

function showTab(which){
  const lv = which === "levels";
  $('#view-levels').hidden = !lv;
  $('#view-ladders').hidden = lv;
  $('#tab-levels').setAttribute("aria-selected", String(lv));
  $('#tab-ladders').setAttribute("aria-selected", String(!lv));
  if (lv) requestAnimationFrame(draw);
}
"""
