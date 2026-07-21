# MF6 I/O Reference Notes

This note captures the `myflopy` assumptions that come directly from the
local MF6 I/O reference PDF:

- Source PDF: `C:\Users\lukem\OneDrive - Associated Earth Sciences Inc\References\Software\MODFLOW\mf6io.pdf`
- Reviewed on: `2026-05-11`

It is intentionally a summary, not a copy of the PDF.

## Core Budget File Structure

Relevant section:
- `mf6io.pdf` pages 365-369

Key points:
- MF6 groundwater-model budget files use `IMETH=6` for stress-package style
  list records.
- For `IMETH=6` records, each row contains:
  - `ID1`
  - `ID2`
  - one budget value
  - zero or more auxiliary values
- For GWF model budget files:
  - `ID1` is the GWF cell/node number
  - `ID2` is usually the package bound number, unless MF6 documents a different
    meaning for that flow type

## Advanced Package Output Budgets

Relevant section:
- `mf6io.pdf` pages 370-376

These package-output budget files are not shaped the same way as the main GWF
budget file. The package-specific tables below are the source of truth.

### LAK Package Output

Relevant section:
- `mf6io.pdf` pages 370-371, Table 51

For LAK package binary output:
- `GWF` means calculated flow from lake (`ID1`) to GWF cell (`ID2`)
- auxiliary `FLOW-AREA` is written for this term
- the package-output `GWF` rows do **not** include the LAK input-file
  `iconn` value

Implications for `myflopy`:
- LAK exchange sign convention (the **feature** reference frame). LAK's `GWF`
  term comes from the LAK *package budget file*, which — like every MF6
  advanced-package budget — is written from the feature's own balance, inflows
  positive:
  - positive `q` = groundwater entering the lake = **lake GAINS**
  - negative `q` = **lake LOSES** to groundwater
  - Verified on the canonical model: the perched lake leaks downward (loses) and
    reports `GWF` in roughly `-5888..0`.
  - myflopy keeps this RAW sign and names the column **`q_lake`** so the frame is
    explicit; the accessor stays `results.q`. See `ResultSpec.reference_frame`.
- LAK package-output rows must be matched back to `lak.connectiondata`
  using:
  - lake id
  - GWF cell
  - within-cell connection order
  - and that order must reset each stress period / time step
- `FLOW-AREA` is the authoritative exchange area coming from MF6 output

### SFR Package Output

Relevant section:
- `mf6io.pdf` pages 374-375, Table 53

For SFR package binary output:
- `GWF` means calculated flow from reach (`ID1`) to GWF cell (`ID2`)
- auxiliary `FLOW-AREA` is written for this term

Implications for `myflopy`:
- SFR exchange sign convention (the **gwf** reference frame). Unlike LAK, SFR's
  exchange here is the aquifer's cell-by-cell (`.cbc`) `SFR` record, written as
  flow FROM the reach TO the GWF cell:
  - positive `q` = flow from reach into the aquifer = **reach LOSES**
  - negative `q` = **reach GAINS** water from the aquifer
  - Verified on the canonical model: the mostly-gaining stream reports `q` in
    roughly `-864..150`, with the gaining reaches negative.
- myflopy keeps this RAW MF6 sign — it does **not** normalize (an earlier
  normalization was reverted on 2026-07-20 because it diverged from every MF6
  file and, when consulted after the flip, inverted the SFR map). Instead the
  ambiguity that a gaining stream and a gaining lake report opposite signs is
  resolved by **naming the reference frame in the column**:
  - SFR and the list BCs → **`q_gwf`** (aquifer frame; gaining is negative)
  - LAK → **`q_lake`** (feature frame; gaining is positive)
  - The accessor stays `results.q` for both; only the DataFrame column and its
    docstrings name the frame. See `ResultSpec.reference_frame`.
- Colours are oriented to the declared frame, via `_exchange_colorscale(frame)`:
  blue on whichever end is gaining (negative for `gwf`, positive for `feature`),
  red on losing. The combined `surface_water` map draws myflopy's own normalized
  `exchange_intensity` field (positive = gaining), so it uses the feature
  orientation.
- SFR exchange map values are more meaningful when normalized by reach length
  rather than plotted as raw volumetric `q`

### UZF Package Output

Relevant section:
- `mf6io.pdf` page 376, Table 54

For UZF package binary output:
- `GWF` means calculated flow from UZF cell (`ID1`) to GWF cell (`ID2`)
- auxiliary `FLOW-AREA` is written for this term

Implications for `myflopy`:
- `UZF-GWRCH`, `UZF-GWET`, and related result helpers should treat package
  output signs and identifiers according to the package-output table, not
  according to the main GWF budget assumptions

## Grid Interpretation Reminder

Relevant section:
- `mf6io.pdf` page 367, Table 48

For `DISV`:
- model budget array dimensions are `NCPL x 1 x NLAY`
- package-output `ID2` values that refer to GWF cells are still node/cell ids
  that need to be normalized carefully for zero-based internal use

## Guidance for `myflopy`

When implementing or debugging package explorers:

1. Treat the MF6 package-output tables in `mf6io.pdf` as the source of truth.
2. Do not assume package-output `ID2` means the same thing across packages.
3. Do not assume package-output rows include input-file identifiers like
   `iconn`; confirm from the MF6 table first.
4. Preserve `FLOW-AREA` from MF6 output whenever it is written.
5. Keep all public `myflopy` cell and period references zero-based, but
   document exactly where MF6 writes one-based ids.
