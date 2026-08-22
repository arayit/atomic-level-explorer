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
            "fast_rung_at", "fast_rung_tau_s", "rung_lifetimes", "closure",
            "best_detuning_meV", "hop_amplitude",
            "relative_amplitude", "order", "pattern", "colours", "hops",
            "max_order", "ground_start",
            "single_colour", "reach_eV_neutral_frame", "n_ladders", "truncated")
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
.lad .pill.warn{border-style:dashed}
.lad .emit{color:var(--accent-ink)}
.lad .chain{display:flex; align-items:flex-start; font-family:ui-monospace,Menlo,Consolas,monospace}
.lad .lv,.lad .ar{display:inline-flex; flex-direction:column; line-height:1.25}
.lad .lv{padding:0 1px}
.lad .lv b{font-weight:normal; font-size:12.5px}
.lad .ar{align-items:center; padding:0 7px; color:var(--ink-3)}
.lad .ar b{font-weight:normal; font-size:11px}
.lad .ar b::after{content:" →"}
.lad .lv i,.lad .ar i{font-style:normal; font-size:10.5px; color:var(--ink-3)}
.lad .terms{color:var(--ink-3); font-size:11.5px; white-space:nowrap}
.lad .none{padding:14px 18px; color:var(--ink-3)}
.lad tbody tr[data-i]{cursor:pointer}
.lad tbody tr:hover{background:var(--panel-2)}
.lad tbody tr.on{background:var(--panel-2)}
.lad tbody tr.on td:first-child{border-left:3px solid var(--ink); padding-left:6px}
.lad tr.diag td{padding:4px 0 12px}
.lad .dwrap{overflow-x:auto; max-width:100%}
.lad .cdiag{display:block}
.lad .cdiag .axis{stroke:var(--rule); stroke-width:1}
.lad .cdiag .lev{stroke:var(--ink); stroke-width:2}
.lad .cdiag .lev.fast{stroke:var(--accent)}
.lad .cdiag .lev.ai{stroke-dasharray:4 3}
.lad .cdiag .ip{stroke:var(--rule); stroke-width:1; stroke-dasharray:2 4}
.lad .cdiag .hop{stroke:var(--ink-3); stroke-width:1}
.lad .cdiag .head{fill:var(--ink-3)}
.lad .cdiag text{font-family:Georgia,"Times New Roman",serif}
.lad .cdiag .ax{font-size:10.5px; fill:var(--ink-3);
  font-family:ui-monospace,Menlo,Consolas,monospace}
.lad .cdiag .term{font-size:11.5px; fill:var(--ink)}
.lad .cdiag .tau,.lad .cdiag .hoplab,.lad .cdiag .hopdet,.lad .cdiag .iplab,.lad .cdiag .cap{
  font-size:10.5px; fill:var(--ink-3)}
"""

BODY = r"""
<section class="lad" id="view-ladders" hidden>
  <div class="intro">
    <p>Chains climb from a level that holds population in hops of 2 or 3 photons of one
      colour, each hop matching a measured level gap and obeying the parity and
      |&Delta;J|&nbsp;&le;&nbsp;m rules for that photon number, with every rung but the last
      outliving the 20&nbsp;ps intra-burst spacing.</p>
    <p>A chain ending above the ionization limit ends in an ion rather than a photon. One
      ending on a bound picosecond level emits instead, and the wavelength is given. Where the
      rung below the top sits less than two photons under the ionization limit, that is noted
      too: the field driving the last hop can ionize that rung instead of climbing it.</p>
  </div>

  <form class="controls" onsubmit="return false">
    <label>Top state
      <select id="lTop">
        <option value="any" selected>either</option>
        <option value="bound">bound &mdash; emits a photon</option>
        <option value="ai">above ionization &mdash; ends in an ion</option>
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
      <th>Tier</th><th>Species</th><th>Z</th><th>Preparation</th><th>Picosecond state</th>
      <th>Lifetimes</th><th>Energy match</th><th>Order</th>
      <th>Pattern</th><th>One colour</th><th>Reach</th>
    </tr></thead>
    <tbody id="lMat"></tbody>
  </table></div>

  <h2>Chains <span class="muted" id="lLadNote"></span></h2>
  <div class="scroll"><table>
    <thead><tr>
      <th>Species</th><th>Order</th><th>Pattern</th><th>Worst match</th><th>Chain</th>
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

/* The ladder CSV rounds level energies to 3 decimals, which is coarse enough to shift a
   per-hop detuning by a millielectronvolt, and it carries no lifetime for the state the
   ladder starts from. Both are already on the page in full precision, in the level data the
   other tab draws, so look them up there and fall back only if the name does not resolve. */
const _lvIdx = {};
function levelLookup(sp, name, approxE){
  const D = DATA[sp];
  if (!D) return null;
  let idx = _lvIdx[sp];
  if (!idx){
    idx = _lvIdx[sp] = {};
    for (const r of D.lv){
      const n = D.confs[r[0]] + " " + D.terms[r[1]] + (r[2] == null ? "" : " J=" + r[2]);
      (idx[n] || (idx[n] = [])).push([r[3], r[4]]);
    }
  }
  const c = idx[name];
  if (!c) return null;
  let best = c[0];
  for (const x of c) if (Math.abs(x[0] - approxE) < Math.abs(best[0] - approxE)) best = x;
  return best;                                   // [energy, lifetime or null]
}

function stateNote(L, i, tau){
  if (i === L.Es.length && L.ai) return "above ionization";
  if (tau != null) return tauText(tau);
  if (i === 0) return L.kind === "ground" ? "ground state" : "metastable";
  return "unpublished";
}

/* Text metrics for the diagram. A canvas 2D context measures the same string in the same
   font the SVG will render it in, so the boxes below are the real ones rather than a guess at
   average glyph width. Falls back to an estimate where no canvas exists. */
const _mc = (() => { try { return document.createElement("canvas").getContext("2d"); }
                     catch (e) { return null; } })();
const SERIF = 'Georgia,"Times New Roman",Times,serif';
const MONO = 'ui-monospace,Menlo,Consolas,monospace';
function textW(t, px, mono){
  if (_mc){ _mc.font = px + "px " + (mono ? MONO : SERIF); return _mc.measureText(t).width; }
  return String(t).length * px * (mono ? 0.6 : 0.52);
}

const overlap = (a, b) => a.x < b.x + b.w && b.x < a.x + a.w &&
                          a.y < b.y + b.h && b.y < a.y + a.h;

/* Segment against axis-aligned box: true if either end is inside, or the segment crosses any
   edge. Used to keep a hop's label off the arrow that hop draws. */
function segHitsBox(s, b){
  const inside = (px, py) => px >= b.x && px <= b.x + b.w && py >= b.y && py <= b.y + b.h;
  if (inside(s.x1, s.y1) || inside(s.x2, s.y2)) return true;
  const cross = (ax, ay, bx, by, cx, cy, dx, dy) => {
    const d = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);
    if (!d) return false;
    const t = ((cx - ax) * (dy - cy) - (cy - ay) * (dx - cx)) / d;
    const u = ((cx - ax) * (by - ay) - (cy - ay) * (bx - ax)) / d;
    return t >= 0 && t <= 1 && u >= 0 && u <= 1;
  };
  const e = [[b.x, b.y, b.x + b.w, b.y], [b.x + b.w, b.y, b.x + b.w, b.y + b.h],
             [b.x + b.w, b.y + b.h, b.x, b.y + b.h], [b.x, b.y + b.h, b.x, b.y]];
  return e.some(q => cross(s.x1, s.y1, s.x2, s.y2, q[0], q[1], q[2], q[3]));
}

/* A Grotrian-style sketch of one chain: the levels it stops on, at their real energies, and
   the hops between them. Inline SVG rather than canvas because it is small, static and
   text-heavy, and because inline SVG inherits the page's theme tokens.

   Level labels are placed first and treated as fixed -- each hangs under its own line and grows
   rightwards inside its own column. Each hop's label is then placed by trying candidate
   positions around its arrow until one collides with nothing: not the level labels, not an
   earlier hop's label, and not any arrow. Hand-tuned offsets only work until a chain arrives
   with a gap the rule did not anticipate. */
function chainDiagram(sp, L){
  const names = [L.start].concat(L.names);
  const Es = [L.startE].concat(L.Es);
  const taus = [null].concat(L.taus);
  for (let i = 0; i < names.length; i++){
    const hit = levelLookup(sp, names[i], Es[i]);
    if (hit){ Es[i] = hit[0]; if (taus[i] == null) taus[i] = hit[1]; }
  }
  const n = Es.length, ip = DATA[sp] && DATA[sp].ion;
  const PITCH = 196, LW = 76;
  const padL = 16, padR = 170, padT = 34, padB = 54;
  const W = padL + (n - 1) * PITCH + LW + padR, H = 330;
  let lo = Math.min.apply(null, Es), hi = Math.max.apply(null, Es);
  const raw = hi - lo || 1;
  /* The limit belongs on the axis only when it is near the climb. For a doubly charged ion it
     can sit 20 eV above the top rung, and drawing it would flatten the ladder into a line. */
  const showIp = ip != null && ip > lo && ip < hi + 0.6 * raw;
  if (showIp) hi = Math.max(hi, ip);
  const pad = (hi - lo || 1) * 0.16;
  lo -= pad; hi += pad;
  const y = e => padT + (H - padT - padB) * (hi - e) / (hi - lo);
  const x = i => padL + i * PITCH;

  const taken = [], segs = [];
  let g = "";
  /* px is the baseline x, py the baseline y; anchor "middle" and "end" shift the box. */
  function put(px, py, t, cls, size, mono, anchor){
    const w = textW(t, size, mono);
    const bx = anchor === "middle" ? px - w / 2 : (anchor === "end" ? px - w : px);
    taken.push({x: bx - 2, y: py - size * 0.82 - 1, w: w + 4, h: size * 1.12 + 2});
    g += '<text x="' + px.toFixed(1) + '" y="' + py.toFixed(1) + '" class="' + cls + '"'
       + (anchor ? ' text-anchor="' + anchor + '"' : "") + ">" + esc(t) + "</text>";
  }
  const short = t => t.length > 26 ? t.slice(0, 25) + "…" : t;

  let head = '<svg class="cdiag" viewBox="0 0 ' + W + " " + H + '" width="' + W + '" height="'
        + H + '" role="img" aria-label="level diagram for one chain">'
        + '<defs><marker id="ah" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7"'
        + ' markerHeight="7" orient="auto"><path d="M0,0 L8,4 L0,8 z" class="head"/></marker>'
        + "</defs>";

  if (showIp){
    head += '<line class="ip" x1="' + padL + '" y1="' + y(ip).toFixed(1) + '" x2="' + (W - 8)
          + '" y2="' + y(ip).toFixed(1) + '"/>';
    put(padL, y(ip) - 5, "ionization limit " + ip.toFixed(2) + " eV", "iplab", 10.5, false);
  }

  for (let i = 0; i < n; i++){
    const yi = y(Es[i]), xi = x(i);
    const ai = i === n - 1 && L.ai;
    const fast = taus[i] != null && taus[i] < 1e-9;
    head += '<line class="lev' + (ai ? " ai" : "") + (fast ? " fast" : "") + '" x1="' + xi
          + '" y1="' + yi.toFixed(1) + '" x2="' + (xi + LW) + '" y2="' + yi.toFixed(1) + '"/>';
    put(xi, yi - 8, Es[i].toFixed(3) + " eV", "ax", 10.5, true);
    put(xi, yi + 16, short(names[i].split(" ").slice(1).join(" ") || names[i]), "term", 11.5);
    put(xi, yi + 29, ai ? "above ionization" : stateNote(L, i, taus[i]), "tau", 10.5);
    if (i) segs.push({x1: x(i - 1) + LW, y1: y(Es[i - 1]), x2: xi, y2: yi});
  }
  for (const s of segs)
    head += '<line class="hop" x1="' + s.x1 + '" y1="' + s.y1.toFixed(1) + '" x2="' + s.x2
          + '" y2="' + s.y2.toFixed(1) + '" marker-end="url(#ah)"/>';

  put(padL, H - 14, sp + " · order " + L.order + " · " + L.ms.join("+") + " photons",
      "cap", 10.5);

  for (let i = 1; i < n; i++){
    const s = segs[i - 1];
    const lab = [L.ms[i - 1] + "γ " + nmLabel(L.nms[i - 1]) + " nm",
                 ((Es[i] - Es[i - 1] - L.ms[i - 1] * HC / L.nms[i - 1]) * 1000 >= 0 ? "+" : "−")
                 + Math.abs((Es[i] - Es[i - 1] - L.ms[i - 1] * HC / L.nms[i - 1]) * 1000)
                       .toFixed(1) + " meV"];
    const w = Math.max(textW(lab[0], 10.5), textW(lab[1], 10.5)), h = 25;
    const mx = (s.x1 + s.x2) / 2, my = (s.y1 + s.y2) / 2;
    const top = Math.min(s.y1, s.y2), bot = Math.max(s.y1, s.y2);
    const len = Math.hypot(s.x2 - s.x1, s.y2 - s.y1) || 1;
    const nx = (s.y2 - s.y1) / len, ny = -(s.x2 - s.x1) / len;   // normal, pointing up-left
    const cand = [[mx, top - 18], [mx, top - 34],
                  [mx + nx * 26, my + ny * 26], [mx + nx * 42, my + ny * 42],
                  [mx - nx * 26, my - ny * 26], [mx - nx * 42, my - ny * 42],
                  [mx, bot + 30], [mx, top - 50], [mx, bot + 46]];
    let pick = cand[0];
    for (const c of cand){
      const box = {x: c[0] - w / 2 - 3, y: c[1] - h / 2 - 3, w: w + 6, h: h + 6};
      if (box.x < 2 || box.x + box.w > W - 2 || box.y < 2 || box.y + box.h > H - 2) continue;
      if (taken.some(t => overlap(box, t))) continue;
      if (segs.some(q => segHitsBox(q, box))) continue;
      pick = c; break;
    }
    put(pick[0], pick[1] - 1, lab[0], "hoplab", 10.5, false, "middle");
    put(pick[0], pick[1] + 11, lab[1], "hopdet", 10.5, false, "middle");
  }
  return head + g + "</svg>";
}

function ladderRow(sp, L, IDX){
  const names = [L.start].concat(L.names);
  const Es = [L.startE].concat(L.Es);
  const taus = [null].concat(L.taus);
  for (let i = 0; i < names.length; i++){
    const hit = levelLookup(sp, names[i], Es[i]);
    if (hit){ Es[i] = hit[0]; if (taus[i] == null) taus[i] = hit[1]; }
  }

  let line = "";
  for (let i = 0; i < Es.length; i++){
    if (i){
      const d = (Es[i] - Es[i - 1] - L.ms[i - 1] * HC / L.nms[i - 1]) * 1000;
      line += '<span class="ar"><b>' + L.ms[i - 1] + "γ " + esc(nmLabel(L.nms[i - 1]))
            + ' nm</b><i>' + (d >= 0 ? "+" : "−") + Math.abs(d).toFixed(1) + " meV</i></span>";
    }
    line += '<span class="lv"><b>' + Es[i].toFixed(3) + '</b><i>'
          + esc(stateNote(L, i, taus[i])) + "</i></span>";
  }
  const terms = esc(names.join("  →  "));

  let top;
  if (L.ai){
    top = '<span class="muted">above ionization</span> '
        + esc(L.Es[L.Es.length - 1].toFixed(2)) + " eV";
  } else if (L.decay && L.fastTau != null && L.fastAt !== "intermediate"){
    top = '<span class="emit">emits ' + esc(nmLabel(HC / L.decay)) + " nm</span> · τ "
        + esc(tauText(L.fastTau));
  } else {
    top = '<span class="muted">bound</span> ' + esc(L.Es[L.Es.length - 1].toFixed(2)) + " eV"
        + (L.fastTau != null ? ' · fast rung τ ' + esc(tauText(L.fastTau)) : "");
  }
  const flag = L.unver ? ' <span class="pill">' + L.unver + ' rung(s) unpublished</span>' : "";
  /* The field that drives the last hop also ionizes the rung it starts from. Warn when that
     rung sits less than two photons below the ionization limit, because then a lower-order
     ionization competes with the hop and the ladder leaks before it reaches the top. */
  let leak = "";
  const ip = DATA[sp] && DATA[sp].ion, parked = Es[Es.length - 2];
  if (ip && parked != null){
    const n = (ip - parked) / (HC / L.nms[L.nms.length - 1]);
    if (n < 2) leak = ' <span class="pill warn">last rung ' + n.toFixed(1)
                    + 'γ below ionization</span>';
  }
  return '<tr data-i="' + IDX + '"><td>' + esc(sp) + (L.kind === "ground" ? ' <span class="pill">ground</span>' : "")
       + '</td><td class="num">' + L.order + '</td><td class="num">' + L.ms.join("+")
       + '</td><td class="num">' + L.det.toFixed(1) + ' meV</td>'
       + '<td><div class="chain">' + line + '</div><div class="terms">' + terms + '</div></td>'
       + '<td>' + top + flag + leak + '</td></tr>';
}

/* Where the ladder's picosecond state sits, and its lifetime where one was measured. Above
   the ionization limit no line list carries a width, so the column reports the fact -- the
   level sits above the limit -- rather than the inference drawn from it. */
const PS_WHERE = {intermediate: "mid-ladder", final: "top", both: "mid-ladder and top"};

function psState(m){
  if (m.evidence === "measured")
    return (PS_WHERE[m.fast_rung_at] || "–") + " · "
         + (m.fast_rung_tau_s ? tauText(+m.fast_rung_tau_s) : "–");
  if (m.evidence === "autoionizing") return "top · above ionization";
  return "–";
}

function matRow(m){
  const cells = [m.tier, m.spectrum, m.Z, m.prep_note, psState(m),
                 m.rung_lifetimes, m.closure + " · " + m.best_detuning_meV + " meV",
                 m.order, m.pattern, m.single_colour,
                 m.reach_eV_neutral_frame ? m.reach_eV_neutral_frame + " eV" : "–"];
  return '<tr class="t' + m.tier + '">'
       + cells.map((c, i) => '<td' + ([2, 8, 9].includes(i) ? ' class="num"' : "") + '>'
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

let LROWS = [], LOPEN = -1;

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
    : '<tr><td colspan="11" class="none">No species has a chain that passes these filters.</td></tr>';
  LROWS = rows; LOPEN = -1;
  $('#lLad').innerHTML = rows.length ? rows.map((r, i) => ladderRow(r[0], r[1], i)).join("")
    : '<tr><td colspan="6" class="none">No chain passes these filters.</td></tr>';
  $('#lMatNote').textContent = mat.length + " of " + MAT.length + " species";
  $('#lLadNote').textContent = rows.length
    + " chains shown, best first — click one to draw it";
  $('#lCount').textContent = rows.length + " chains · " + live.size + " species";
}

function initLadders(){
  const sel = $('#lSp');
  sel.innerHTML = '<option value="any">all</option>'
    + MAT.filter(m => LAD[m.spectrum])
         .slice().sort((a, b) => a.spectrum.localeCompare(b.spectrum)).map(m =>
        '<option value="' + esc(m.spectrum) + '">' + esc(m.spectrum) + "</option>").join("");
  ["#lTop", "#lDet", "#lCol", "#lVer", "#lStart", "#lSp"].forEach(sel =>
    $(sel).addEventListener("change", drawLadders));
  $('#lLad').addEventListener("click", e => {
    const tr = e.target.closest("tr[data-i]");
    if (!tr) return;
    const i = +tr.dataset.i, again = LOPEN === i;
    const open = $('#lLad').querySelector("tr.diag");
    if (open) open.remove();
    const lit = $('#lLad').querySelector("tr.on");
    if (lit) lit.classList.remove("on");
    if (again){ LOPEN = -1; return; }
    LOPEN = i;
    tr.classList.add("on");
    tr.insertAdjacentHTML("afterend", '<tr class="diag"><td colspan="6"><div class="dwrap">'
      + chainDiagram(LROWS[i][0], LROWS[i][1]) + "</div></td></tr>");
  });
  drawLadders();
}

/* Boot lives at the very end of the combined script. MAT and LAD are const, so the level
   view's boot block -- which runs earlier in the same script -- cannot touch them: reading a
   const before its declaration executes is a ReferenceError, not a hoisted undefined. */
function showTab(which){
  const lv = which === "levels";
  $('#view-levels').hidden = !lv;
  $('#view-ladders').hidden = lv;
  $('#tab-levels').setAttribute("aria-selected", String(lv));
  $('#tab-ladders').setAttribute("aria-selected", String(!lv));
  if (lv) requestAnimationFrame(draw);
}

$('#tab-levels').addEventListener('click', () => showTab('levels'));
$('#tab-ladders').addEventListener('click', () => showTab('ladders'));
initLadders();
"""
