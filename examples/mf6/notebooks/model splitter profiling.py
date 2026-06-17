from pathlib import Path
import pandas as pd
import myflopy as mf

root = Path('../artifacts/canonical_prt_parallel')
config = mf.CanonicalModelConfig.validation()
model = mf.build_canonical_model(root / 'gwf', config=config)

mask = mf.canonical_partition_mask(model, 8)

split = model.parallel.split_model(
    workspace=Path("split_profile"),
    mask=mask,
    write=False,
)
