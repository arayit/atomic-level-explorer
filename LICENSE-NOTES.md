# Licensing

**Code** (`scripts/`) — MIT.

**Data** (`energy_levels/`) is derived from two public services and carries their terms:

- **NIST Atomic Spectra Database** (levels, transition probabilities). A work of the US
  Government; free to use, with the citation given in `README.md`.
- **ExoMol LiDB** (lifetime cross-check for neutrals and singly-charged ions), released under
  **CC BY-SA 4.0**. Files under `energy_levels/*/lidb_*.csv`, and the small number of lifetimes
  in `levels_with_lifetimes.csv` whose `tau_source` column reads `LiDB (J-lumped)`, inherit
  that licence.

**Vapour-pressure temperatures** quoted in the page come from the CRC Handbook of Chemistry and
Physics, 84th ed. (metals after Alcock, Itkin & Horrigan 1984) and are cited, not redistributed.
