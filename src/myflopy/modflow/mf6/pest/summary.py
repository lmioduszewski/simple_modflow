"""Human-readable configuration summary for a :class:`PestProject`.

The calibration facade hides pyEMU boilerplate, but it must never hide *state*.
``cal.settings()`` returns one of these objects so the resolved configuration --
every parameter, observation, and forecast, plus the built control-file counts --
can be reviewed at a glance before or after building.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class PestSettings:
    """A printable snapshot of a calibration project's resolved configuration."""

    name: str
    model_name: str
    original_workspace: str
    template_workspace: str
    start_datetime: str
    parameters: list[dict] = field(default_factory=list)
    observations: list[dict] = field(default_factory=list)
    forecasts: list[dict] = field(default_factory=list)
    built: bool = False
    npar: int | None = None
    npar_groups: int | None = None
    nobs: int | None = None
    nnz_obs: int | None = None
    n_forecasts: int | None = None
    noptmax: int | None = None

    def parameter_frame(self) -> pd.DataFrame:
        """Return the declared parameters as a tidy table."""

        return pd.DataFrame(self.parameters)

    def observation_frame(self) -> pd.DataFrame:
        """Return the declared observations as a tidy table."""

        return pd.DataFrame(self.observations)

    def __str__(self) -> str:
        lines: list[str] = []
        lines.append(f"PEST calibration: {self.name}  (model: {self.model_name})")
        lines.append(f"  template : {self.template_workspace}")
        lines.append(f"  start    : {self.start_datetime}")

        lines.append(f"  parameters ({len(self.parameters)}):")
        if self.parameters:
            for par in self.parameters:
                phys = par.get("physical")
                phys_txt = f" physical={tuple(phys)}" if phys else ""
                lines.append(
                    f"    - {par['target']:<10} style={par['style']:<11} "
                    f"bounds={tuple(par['bounds'])}{phys_txt} "
                    f"transform={par['transform']}"
                )
        else:
            lines.append("    (none)")

        lines.append(f"  observations ({len(self.observations)}):")
        if self.observations:
            for obs in self.observations:
                lines.append(
                    f"    - {obs.get('prefix', '?'):<10} kind={obs.get('kind', '?'):<12} "
                    f"n={obs.get('n', '?')}"
                )
        else:
            lines.append("    (none)")

        if self.forecasts:
            lines.append(f"  forecasts ({len(self.forecasts)}):")
            for fore in self.forecasts:
                lines.append(f"    - {fore.get('prefix', '?'):<10} n={fore.get('n', '?')}")

        if self.built:
            lines.append("  built control file:")
            lines.append(
                f"    npar={self.npar} (groups={self.npar_groups})  "
                f"nobs={self.nobs} (nonzero weight={self.nnz_obs})  "
                f"forecasts={self.n_forecasts}  noptmax={self.noptmax}"
            )
        else:
            lines.append("  (not built yet -- call cal.build())")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.__str__()
