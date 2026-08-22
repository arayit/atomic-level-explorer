# Atomic level explorer

Energy levels, radiative lifetimes and linear absorption strengths for **238 atomic spectra**
across **45 elements**, assembled from public databases and served as a single self-contained page.

**→ [Open the explorer](https://arayit.github.io/atomic-level-explorer/)**

Built for a screening question: which atomic states decay on a **picosecond** timescale — fast
enough to matter between the pulses of a GHz burst, slow enough that an ordinary nanosecond laser
could not reach them — and what would it cost to prepare the species that carries them.

| | |
|---|---|
| spectra | 238 |
| energy levels | 44,654 |
| levels with a lifetime in 1–1000 ps | 4,797 |
| transitions with an oscillator strength | 82,407 |

## What the page does

- Two columns, **even parity** left and **odd parity** right. An *m*-photon hop changes parity by
  (−1)^*m*, so a two-photon step stays in its column and a one- or three-photon step crosses.
- **Semantic zoom.** Zoomed out, only levels carrying a lifetime are drawn; zoom in and the rest
  appear. A density strip at the axis always shows every level, so nothing is hidden.
- **Energy reference switch.** A level of Al II sits 7.42 eV above *its own ion's* ground state,
  but 13.4 eV above neutral aluminium — the ionization had to be paid first. The switch shows both.
- **Click a level** to measure every gap from it, with photon-hop guides drawn only over the
  column each hop can legally reach, and marked where a real level actually sits there.
- **Photon wavelength** is typed in, 1–2000 nm; the hop grid recomputes.

## Lifetimes are computed, not looked up

No database tabulates them. They are built here as

    τ_k = 1 / Σ A_ki

summed over every line in NIST ASD with *k* as the upper level, J-resolved. Where a level's decay
channels are only partly covered by published A-values, its τ is an **upper bound**. Levels above
the first ionization limit autoionize; their widths are in no database, so they are left blank
rather than filled with a radiative value that does not apply.

Oscillator strengths come from the same A-values,
`f = 1.4992e-16 (g_k/g_i) A_ki λ²` with λ in ångström, and are **one-photon** quantities.

## Layout

    index.html            the site, self-contained, no external requests
    scripts/scrape.py     pulls NIST ASD levels + lines and ExoMol LiDB, caching raw responses
    scripts/transitions.py derives f, gf and ∫σdν from the cached A-values
    scripts/diagram.py    one PDF per spectrum, levels faster than a threshold
    scripts/explorer.py   builds index.html from the CSVs
    scripts/species.txt   the spectrum list — add a line to extend it
    energy_levels/<Sp>/   raw responses, derived tables, limits, meta.json, diagram

## Reproducing

    python3 scripts/scrape.py        # only fetches what is not already cached
    python3 scripts/transitions.py
    python3 scripts/diagram.py
    python3 scripts/explorer.py

## Sources

Kramida A., Ralchenko Yu., Reader J. and NIST ASD Team (2024). *NIST Atomic Spectra Database*
(version 5.12). National Institute of Standards and Technology, Gaithersburg MD.
[doi:10.18434/T4W30F](https://doi.org/10.18434/T4W30F)

Owens A., Chen T.-T., Hill C., Mohr S. and Tennyson J. (2025). LiDB: Database of atomic radiative
lifetimes for plasma processes. *J. Quant. Spectrosc. Radiat. Transfer* **330**, 109242.
LiDB is not an independent check — its atomic lifetimes are themselves derived from NIST ASD.

Vapour temperatures: *CRC Handbook of Chemistry and Physics*, 84th ed.; metals after Alcock,
Itkin & Horrigan (1984).

See [LICENSE-NOTES.md](LICENSE-NOTES.md).
