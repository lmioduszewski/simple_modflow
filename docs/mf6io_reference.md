# MF6 I/O Reference Notes

This note captures the `simple_modflow` assumptions that come directly from the
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

Implications for `simple_modflow`:
- LAK exchange sign convention:
  - positive `q` = lake losing to groundwater
  - negative `q` = groundwater gaining to lake
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

Implications for `simple_modflow`:
- SFR exchange sign convention:
  - positive `q` = stream losing to groundwater
  - negative `q` = groundwater gaining to stream
- raw SFR choropleths should therefore use:
  - blue for gaining reaches
  - red for losing reaches
- SFR exchange map values are more meaningful when normalized by reach length
  rather than plotted as raw volumetric `q`

### UZF Package Output

Relevant section:
- `mf6io.pdf` page 376, Table 54

For UZF package binary output:
- `GWF` means calculated flow from UZF cell (`ID1`) to GWF cell (`ID2`)
- auxiliary `FLOW-AREA` is written for this term

Implications for `simple_modflow`:
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

## Guidance for `simple_modflow`

When implementing or debugging package explorers:

1. Treat the MF6 package-output tables in `mf6io.pdf` as the source of truth.
2. Do not assume package-output `ID2` means the same thing across packages.
3. Do not assume package-output rows include input-file identifiers like
   `iconn`; confirm from the MF6 table first.
4. Preserve `FLOW-AREA` from MF6 output whenever it is written.
5. Keep all public `simple_modflow` cell and period references zero-based, but
   document exactly where MF6 writes one-based ids.
