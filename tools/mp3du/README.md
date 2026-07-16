# MP3DU executables (untracked)

Place the mod-PATH3DU binaries here (implementation plan 2.4) — they are no
longer shipped inside the `myflopy` package:

- `mp3du.exe`
- `writep3dgsf.exe`
- `writep3doutput.exe`

`ParticleTrackingInput` resolves them in this order:

1. explicit `*_path` constructor arguments
2. the `MYFLOPY_MP3DU_DIR` environment variable (a directory)
3. this directory (`<repo>/tools/mp3du/`)
4. deprecated fallback: inside the installed package (warns; will be removed)

The binaries are Windows executables from the S.S. Papadopulos & Associates
mod-PATH3DU distribution — download from https://www.sspa.com/software/mod-path3du.
