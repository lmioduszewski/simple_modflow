"""Thin package-builder wrappers around common FloPy MF6 package classes.

These helpers keep model-assembly code concise and now also support optional
project-catalog artifact capture during package creation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import flopy
import numpy as np

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def _gwf_package_filename(model: "SimulationBase", suffix: str) -> str:
    """Return the standard MF6 package filename for ``model`` and ``suffix``."""

    return f"{model.name}.{suffix}"


def _build_gwf_package(
    constructor,
    model: "SimulationBase",
    *,
    pname: str,
    filename_suffix: str,
    **kwargs,
):
    """Create one FloPy GWF package using the shared model and naming conventions."""

    return constructor(
        model.gwf,
        pname=pname,
        filename=_gwf_package_filename(model, filename_suffix),
        **kwargs,
    )


class OutputControl:
    """Create the standard MF6 output-control package for a model."""

    def __init__(
            self,
            model: "SimulationBase",
            save_record=(("HEAD", "LAST"), ("BUDGET", "LAST")),
            print_record=None
    ):
        """Parameters
        ----------
        model
            Target model receiving the output-control package.
        save_record
            FloPy/MF6 save-record tuples, defaulting to final heads and budgets.
        print_record
            Optional MF6 print-record tuples.
        """
        head_file = f"{model.name}.hds"
        budget_file = f"{model.name}.cbc"
        package = _build_gwf_package(
            flopy.mf6.modflow.ModflowGwfoc,
            model,
            pname="oc",
            filename_suffix="oc",
            saverecord=save_record,
            head_filerecord=head_file,
            budget_filerecord=budget_file,
            printrecord=print_record,
        )
        _finalize_wrapper(self, attr_name="oc", package=package, model=model)


def _maybe_create_package_artifact(
    model: "SimulationBase",
    package_name: str,
    *,
    artifact_id: str | None = None,
    artifact_catalog=None,
    artifact_description: str | None = None,
    artifact_tags: list[str] | None = None,
    artifact_metadata: dict | None = None,
    artifact_overwrite: bool = False,
):
    """Optionally capture a reusable package artifact while building a package.

    If an ``artifact_id`` and explicit artifact store are provided, the package
    is captured after construction.
    """

    if artifact_id is None:
        return None

    catalog = artifact_catalog
    if catalog is None:
        raise ValueError("artifact_catalog is required when artifact_id is provided.")

    artifact = catalog.create_package_artifact(
        artifact_id,
        model=model,
        package_name=package_name,
        description=artifact_description,
        tags=artifact_tags,
        metadata=artifact_metadata,
        overwrite=artifact_overwrite,
    )
    return artifact


def _finalize_wrapper(
    instance,
    *,
    attr_name: str,
    package,
    model: "SimulationBase",
    package_name: str | None = None,
    artifact_id: str | None = None,
    artifact_catalog=None,
    artifact_description: str | None = None,
    artifact_tags: list[str] | None = None,
    artifact_metadata: dict | None = None,
    artifact_overwrite: bool = False,
):
    """Attach the built FloPy package and optionally capture a package artifact."""

    setattr(instance, attr_name, package)
    instance.package_artifact = None
    if package_name is None:
        return

    instance.package_artifact = _maybe_create_package_artifact(
        model,
        package_name,
        artifact_id=artifact_id,
        artifact_catalog=artifact_catalog,
        artifact_description=artifact_description,
        artifact_tags=artifact_tags,
        artifact_metadata=artifact_metadata,
        artifact_overwrite=artifact_overwrite,
    )


class InitialConditions:
    """Create the MF6 initial-conditions package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            vor: "Vor",
            botm_cells: list = None,
            initial_sat_thickness: float = 0.5,
            nlay=1,
            strt=None,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the package.
        vor
            Grid helper used to infer default cell counts.
        botm_cells
            Bottom elevations used when deriving default starting heads.
        initial_sat_thickness
            Added to ``botm_cells`` when ``strt`` is not supplied.
        nlay
            Number of layers represented in ``botm_cells``.
        strt
            Optional explicit starting heads.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        if botm_cells is None:
            botm_cells = [0 for _ in range(vor.ncpl * nlay)]
        if strt is None:
            strt = [cell_elev + initial_sat_thickness for cell_elev in botm_cells]
        package = _build_gwf_package(
            flopy.mf6.modflow.mfgwfic.ModflowGwfic,
            model,
            pname="ic",
            filename_suffix="ic",
            strt=strt,
        )
        _finalize_wrapper(
            self,
            attr_name="ic",
            package=package,
            model=model,
            package_name="ic",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class KFlow:
    """Create the MF6 NPF package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            k: list = None,
            k33_vert=None,
            icelltype=1,
            perched: bool = False,
            save_specific_discharge: bool = True,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,

    ):
        """Parameters
        ----------
        model
            Target model receiving the NPF package.
        k
            Horizontal hydraulic conductivity input.
        k33_vert
            Optional vertical conductivity input.
        icelltype
            Convertible/confined flag by layer or cell. The canonical four-layer
            model uses ``[1, 1, 0, 0]``.
        perched
            Whether perched conditions are enabled.
        save_specific_discharge
            Whether MF6 should save specific discharge outputs.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        package = _build_gwf_package(
            flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf,
            model,
            pname="npf",
            filename_suffix="npf",
            icelltype=icelltype,
            k=k,
            perched=perched,
            k33=k33_vert,
            save_flows=True,
            save_saturation=True,
            save_specific_discharge=save_specific_discharge,
        )
        _finalize_wrapper(
            self,
            attr_name="npf",
            package=package,
            model=model,
            package_name="npf",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class Storage:
    """Create the MF6 storage package for steady/transient simulations."""

    def __init__(
            self,
            model: "SimulationBase",
            specific_storage: float = 0.0001,
            specific_yield: float = 0.2,
            sto_steady: dict = None,
            sto_transient: dict = None,
            iconvert=1,

    ):
        """Parameters
        ----------
        model
            Target model receiving the storage package.
        specific_storage
            Specific storage value for active cells.
        specific_yield
            Specific yield value for convertible cells.
        sto_steady, sto_transient
            MF6 steady-state and transient stress-period flags.
        iconvert
            Convertible/confined storage flag by layer or cell.
        """
        if sto_steady is None and sto_transient is None:
            sto_steady = {0: True}
        if sto_transient is None:
            sto_transient = {1: True}
        package = _build_gwf_package(
            flopy.mf6.ModflowGwfsto,
            model,
            pname="sto",
            filename_suffix="sto",
            save_flows=True,
            iconvert=iconvert,
            ss=specific_storage,
            sy=specific_yield,
            steady_state=sto_steady,
            transient=sto_transient,
        )
        _finalize_wrapper(self, attr_name="sto", package=package, model=model)


class Recharge:
    """Create the MF6 recharge package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            vor: "Vor" = None,
            rch_dict: dict = None,
            auxiliary: list[str] = None,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the recharge package.
        vor
            Grid helper used to determine ``maxbound`` when needed.
        rch_dict
            MF6 recharge stress-period data dictionary.
        auxiliary
            Optional auxiliary variable names.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        vor = model.vor if vor is None else vor
        package = _build_gwf_package(
            flopy.mf6.ModflowGwfrch,
            model,
            pname="rch",
            filename_suffix="rch",
            print_input=False,
            print_flows=False,
            save_flows=True,
            maxbound=len(vor.iverts),
            stress_period_data=rch_dict,
            auxiliary=auxiliary,
        )
        _finalize_wrapper(
            self,
            attr_name="rch",
            package=package,
            model=model,
            package_name="rch",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class Wells:
    """Legacy WEL-package wrapper that attaches a well package to a model object.

    Part of the legacy OO package layer used with :class:`SimulationBase`:
    constructing it builds a :class:`flopy.mf6.ModflowGwfwel` package (with
    myflopy's shared naming/``save_flows`` conventions and the project artifact
    catalog), attaches it to ``model`` as ``model.wel``, and exposes the wrapped
    package. For new work prefer the package-first ``mf.wel(...)`` /
    :func:`~myflopy.advanced.wel_spec` path instead.

    Parameters
    ----------
    model
        The :class:`SimulationBase` to attach the WEL package to.
    stress_period_data
        FloPy WEL stress-period data ``{per: [(cellid, q, ...), ...]}``.
    auxiliary
        Optional auxiliary variable name(s).
    boundnames
        Enable named boundaries (default ``True``).
    artifact_id, artifact_catalog, artifact_description, artifact_tags, artifact_metadata, artifact_overwrite
        Optional project artifact-catalog bookkeeping for the built package.
    """

    def __init__(
            self,
            model: "SimulationBase",
            stress_period_data,
            auxiliary=None,
            boundnames: bool = True,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Build and wrap a GWF WEL package from ``stress_period_data`` (see the class docstring)."""

        package = _build_gwf_package(
            flopy.mf6.ModflowGwfwel,
            model,
            pname="wel",
            filename_suffix="wel",
            save_flows=True,
            stress_period_data=stress_period_data,
            auxiliary=auxiliary,
            boundnames=boundnames,
        )
        _finalize_wrapper(
            self,
            attr_name="wel",
            package=package,
            model=model,
            package_name="wel",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class Drains:
    """Create the MF6 drain package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            stress_period_data: list,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the drain package.
        stress_period_data
            MF6 drain stress-period data.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        package = _build_gwf_package(
            flopy.mf6.ModflowGwfdrn,
            model,
            pname="drn",
            filename_suffix="drn",
            save_flows=True,
            print_flows=False,
            print_input=False,
            stress_period_data=stress_period_data,
        )
        _finalize_wrapper(
            self,
            attr_name="drn",
            package=package,
            model=model,
            package_name="drn",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class GHB:
    """Create the MF6 general-head-boundary package, optionally as an artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            stress_period_data,
            auxiliary=None,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the GHB package.
        stress_period_data
            MF6 GHB stress-period data.
        auxiliary
            Optional auxiliary variable names.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        package = _build_gwf_package(
            flopy.mf6.ModflowGwfghb,
            model,
            pname="ghb",
            filename_suffix="ghb",
            print_input=False,
            print_flows=False,
            save_flows=True,
            stress_period_data=stress_period_data,
            auxiliary=auxiliary,
        )
        _finalize_wrapper(
            self,
            attr_name="ghb",
            package=package,
            model=model,
            package_name="ghb",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class CHD:
    """Create the MF6 constant-head package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            stress_period_data,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the CHD package.
        stress_period_data
            MF6 CHD stress-period data.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        package = _build_gwf_package(
            flopy.mf6.ModflowGwfchd,
            model,
            pname="chd",
            filename_suffix="chd",
            print_input=False,
            print_flows=False,
            save_flows=True,
            stress_period_data=stress_period_data,
        )
        _finalize_wrapper(
            self,
            attr_name="chd",
            package=package,
            model=model,
            package_name="chd",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )


class UZF:
    """Create the MF6 UZF package, optionally as a reusable artifact."""

    def __init__(
            self,
            model: "SimulationBase",
            packagedata=None,
            perioddata=None,
            print_input=False,
            print_flows=True,
            save_flows=True,
            mover=False,
            simulate_et=False,
            linear_gwet=False,
            square_gwet=False,
            simulate_gwseep=False,
            unsat_etwc=False,
            unsat_etae=False,
            nuzfcells=None,
            ntrailwaves=7,
            nwavesets=40,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model receiving the UZF package.
        packagedata, perioddata
            FloPy/MF6 UZF package inputs.
        print_input, print_flows, save_flows
            Standard MF6 UZF control flags.
        mover, simulate_et, linear_gwet, square_gwet, simulate_gwseep,
        unsat_etwc, unsat_etae
            MF6 UZF option flags.
        nuzfcells, ntrailwaves, nwavesets
            Core UZF sizing and wave-tracking parameters.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """

        if nuzfcells is None:
            nuzfcells = int(np.bincount(model.modelgrid.idomain[0])[1])

        self.uzf = flopy.mf6.ModflowGwfuzf(
            model=model.gwf,
            print_input=print_input,
            print_flows=print_flows,
            save_flows=save_flows,
            budget_filerecord=f'{model.name}_budget.uzf',
            budgetcsv_filerecord=f'{model.name}_uzf_budget.csv',
            package_convergence_filerecord=f'{model.name}_uzf_package_convergence.csv',
            mover=mover,
            simulate_et=simulate_et,
            linear_gwet=linear_gwet,
            square_gwet=square_gwet,
            simulate_gwseep=simulate_gwseep,
            unsat_etwc=unsat_etwc,
            unsat_etae=unsat_etae,
            nuzfcells=nuzfcells,
            ntrailwaves=ntrailwaves,
            nwavesets=nwavesets,
            packagedata=packagedata,
            perioddata=perioddata,
            filename=f'{model.name}.uzf',
            pname='uzf',
        )
        self.package_artifact = _maybe_create_package_artifact(
            model,
            "uzf",
            artifact_id=artifact_id,
            artifact_catalog=artifact_catalog,
            artifact_description=artifact_description,
            artifact_tags=artifact_tags,
            artifact_metadata=artifact_metadata,
            artifact_overwrite=artifact_overwrite,
        )
