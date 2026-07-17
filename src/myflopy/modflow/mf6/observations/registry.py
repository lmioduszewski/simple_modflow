"""`TargetRegistry` — binds targets to a model (`model.targets`)."""

from __future__ import annotations

import pandas as pd

from myflopy.modflow.mf6.observations.drn import BoundDrnFlowTargets, DrnFlowTargets
from myflopy.modflow.mf6.observations.heads import BoundHeadTargets, HeadTargets
from myflopy.modflow.mf6.observations.lake import BoundLakeStageTargets, LakeStageTargets
from myflopy.modflow.mf6.observations.sfr import (
    BoundSfrFlowTargets,
    BoundSfrStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)


class TargetRegistry:
    """Lightweight model-bound registry for reusable calibration targets."""

    _INTERNAL_NAMES = {"model", "_targets"}

    def __init__(self, model):
        """Create an empty target registry bound to ``model``."""

        object.__setattr__(self, "model", model)
        object.__setattr__(self, "_targets", {})

    def _bind(self, target):
        """Wrap a bare target in its model-bound helper (pass through anything unrecognized)."""

        if isinstance(target, HeadTargets):
            return BoundHeadTargets(self.model, target)
        if isinstance(target, LakeStageTargets):
            return BoundLakeStageTargets(self.model, target)
        if isinstance(target, SfrStageTargets):
            return BoundSfrStageTargets(self.model, target)
        if isinstance(target, SfrFlowTargets):
            return BoundSfrFlowTargets(self.model, target)
        if isinstance(target, DrnFlowTargets):
            return BoundDrnFlowTargets(self.model, target)
        return target

    def _coerce_target(self, value):
        """Unwrap a bound helper back to its bare target, validating the type (raises otherwise)."""

        if isinstance(
            value,
            (
                BoundHeadTargets,
                BoundLakeStageTargets,
                BoundSfrStageTargets,
                BoundSfrFlowTargets,
                BoundDrnFlowTargets,
            ),
        ):
            return value.targets
        if isinstance(value, (HeadTargets, LakeStageTargets, SfrStageTargets, SfrFlowTargets, DrnFlowTargets)):
            return value
        raise TypeError(
            "Model targets currently support HeadTargets, LakeStageTargets, "
            "SfrStageTargets, SfrFlowTargets, and DrnFlowTargets."
        )

    def keys(self) -> list[str]:
        """The registered target-set names, sorted."""

        return sorted(self._targets)

    def summary(self) -> pd.DataFrame:
        """A ``name``/``type`` table of every registered target set."""

        rows = [
            {"name": name, "type": type(target).__name__}
            for name, target in sorted(self._targets.items())
        ]
        return pd.DataFrame(rows)

    def _named_target(self, name: str):
        """The model-bound helper for a registered target ``name`` (raises if absent)."""

        if name not in self._targets:
            raise AttributeError(f"{type(self).__name__!r} has no target set {name!r}")
        return self[name]

    @property
    def heads(self) -> BoundHeadTargets:
        """Return model-bound head targets with IDE-visible completion."""

        return self._named_target("heads")

    @heads.setter
    def heads(self, value: HeadTargets | BoundHeadTargets):
        """Register ``value`` as the ``heads`` target set."""

        self["heads"] = value

    @property
    def lake_stage(self) -> BoundLakeStageTargets:
        """Return model-bound lake-stage targets."""

        return self._named_target("lake_stage")

    @lake_stage.setter
    def lake_stage(self, value: LakeStageTargets | BoundLakeStageTargets):
        """Register ``value`` as the ``lake_stage`` target set."""

        self["lake_stage"] = value

    @property
    def sfr_stage(self) -> BoundSfrStageTargets:
        """Return model-bound SFR-stage targets."""

        return self._named_target("sfr_stage")

    @sfr_stage.setter
    def sfr_stage(self, value: SfrStageTargets | BoundSfrStageTargets):
        """Register ``value`` as the ``sfr_stage`` target set."""

        self["sfr_stage"] = value

    @property
    def sfr_flow(self) -> BoundSfrFlowTargets:
        """Return model-bound SFR-flow targets."""

        return self._named_target("sfr_flow")

    @sfr_flow.setter
    def sfr_flow(self, value: SfrFlowTargets | BoundSfrFlowTargets):
        """Register ``value`` as the ``sfr_flow`` target set."""

        self["sfr_flow"] = value

    @property
    def drn_flow(self) -> BoundDrnFlowTargets:
        """Return model-bound DRN seepage-flow targets."""

        return self._named_target("drn_flow")

    @drn_flow.setter
    def drn_flow(self, value: DrnFlowTargets | BoundDrnFlowTargets):
        """Register ``value`` as the ``drn_flow`` target set."""

        self["drn_flow"] = value

    def __getitem__(self, key: str):
        """Return the model-bound helper for the registered target ``key``."""

        return self._bind(self._targets[key])

    def __setitem__(self, key: str, value):
        """Register ``value`` (coerced to a bare target) under name ``key``."""

        self._targets[str(key)] = self._coerce_target(value)

    def __contains__(self, key: str) -> bool:
        """Whether a target set named ``key`` is registered."""

        return str(key) in self._targets

    def __delitem__(self, key: str):
        """Remove the registered target set named ``key``."""

        del self._targets[str(key)]

    def __getattr__(self, name: str):
        """Resolve ``registry.<name>`` to a registered target's bound helper (else ``AttributeError``)."""

        targets = object.__getattribute__(self, "__dict__").get("_targets", {})
        if name in targets:
            return self._bind(targets[name])
        raise AttributeError(f"{type(self).__name__!r} has no target set {name!r}")

    def __setattr__(self, name: str, value):
        """Assigning ``registry.<name> = targets`` registers a target set (internal fields pass through)."""

        if name in self._INTERNAL_NAMES:
            object.__setattr__(self, name, value)
            return
        targets = object.__getattribute__(self, "__dict__").get("_targets")
        if targets is None:
            object.__setattr__(self, name, value)
            return
        targets[name] = self._coerce_target(value)

    def __dir__(self):
        """Advertise the registered target-set names for tab-completion."""

        targets = object.__getattribute__(self, "__dict__").get("_targets", {})
        return sorted(set(super().__dir__()) | set(targets))
