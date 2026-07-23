from __future__ import annotations

from pathlib import Path
import io
import keyword
import token as token_mod
import tokenize

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import Paragraph
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "myflopy_api_pamphlet.pdf"


PALETTE = {
    "navy": colors.HexColor("#12324A"),
    "teal": colors.HexColor("#1E7A78"),
    "sand": colors.HexColor("#F5F1E8"),
    "slate": colors.HexColor("#475569"),
    "ink": colors.HexColor("#10212F"),
    "gold": colors.HexColor("#C79B3B"),
    "mist": colors.HexColor("#E8F0F2"),
    "rose": colors.HexColor("#E9D7D4"),
    "leaf": colors.HexColor("#DCEAD9"),
    "sky": colors.HexColor("#DDEBF7"),
    "code_bg": colors.HexColor("#1E1E1E"),
    "code_border": colors.HexColor("#2F3F4A"),
    "code_gutter": colors.HexColor("#252526"),
    "code_text": colors.HexColor("#D4D4D4"),
    "code_keyword": colors.HexColor("#C586C0"),
    "code_string": colors.HexColor("#CE9178"),
    "code_number": colors.HexColor("#B5CEA8"),
    "code_name": colors.HexColor("#9CDCFE"),
    "code_func": colors.HexColor("#DCDCAA"),
    "code_class": colors.HexColor("#4EC9B0"),
    "code_comment": colors.HexColor("#6A9955"),
    "code_operator": colors.HexColor("#D4D4D4"),
    "code_module": colors.HexColor("#4FC1FF"),
    "code_accent": colors.HexColor("#569CD6"),
}


PAGES = [
    {
        "title": "myflopy",
        "subtitle": "A one-page-at-a-time guide to the main workflows in the modern API.",
        "accent": "navy",
        "tagline": "Build, inspect, compare, calibrate, and reopen MF6 models from one consistent surface.",
        "cards": [
            ("Build", "Voronoi / DISV helpers, package builders, and a stable `SimulationBase` workflow."),
            ("Explore", "Heads, budgets, package explorers, cross-sections, and surface-water summaries."),
            ("Compare", "Model groups, run explorers, project catalogs, and consistent grouped plotting."),
            ("Calibrate", "Reusable targets, FloPy observations, PEST workspaces, and reopenable review."),
        ],
    },
    {
        "title": "1. Grid & Build",
        "subtitle": "Voronoi-first model setup with one durable model object at the center.",
        "accent": "teal",
        "narrative": "The build workflow starts with geometry, not boilerplate. `VoronoiGridPlus` carries the spatial intelligence, while `SimulationBase` becomes the long-lived object that packages, plotting, and review all attach to later.",
        "takeaway": "Treat the grid helper and the model object as the stable backbone of the whole workflow.",
        "bullets": [
            "Use `VoronoiGridPlus` for polygon geometry, cached adjacency, centroid sampling, and feature-to-cell mapping.",
            "Attach discretization and core packages with `DisvGrid`, `TemporalDiscretization`, `InitialConditions`, `KFlow`, `Storage`, and `OutputControl`.",
            "Keep GIS inputs central with builders like `DRNFromVector`, `GHBFromVector`, `CHDFromVector`, `RCHFromVector`, and `KFromVector`.",
            "Run once, then keep inspecting the same object through `model.packages`, `model.outputs`, `model.hds`, and `model.bud()`.",
        ],
        "code": [
            "# start with a grid helper you can keep using later",
            "vor = mf.VoronoiGridPlus(triangle)",
            "cells = vor.get_vor_cells_as_dict(locs, loc_name_field='name')",
            "",
            "# build one durable model object",
            "model = mf.SimulationBase(",
            "    name='demo',",
            "    vor=vor,",
            "    nper=2,",
            ")",
            "",
            "# attach discretization and core flow properties",
            "mf.DisvGrid(vor=vor, model=model, ...)",
            "mf.KFlow(model=model, k=k_array, k33_vert=k33)",
            "mf.Storage(model=model, sy=0.15, ss=1e-5)",
            "mf.OutputControl(model=model, budget_filerecord='demo.cbb')",
            "",
            "# run and keep using the same object",
            "success, buff = model.run_simulation()",
        ],
        "api_notes": [
            "`VoronoiGridPlus` is more than a mesh container; it is the shared geometry helper used by selection, plotting, and input builders.",
            "`SimulationBase` is the anchor object; the package-level builders mainly feed data into it rather than replacing it.",
        ],
    },
    {
        "title": "2. Heads & Results",
        "subtitle": "Use the improved heads/results viewers for maps, tables, and sections.",
        "accent": "gold",
        "narrative": "The results layer is meant to feel immediate: load heads once, inspect them as tables, maps, or sections, and move between spatial and tabular views without changing mental models.",
        "takeaway": "The heads workflow should make it easy to move from raw results to a map or cross-section in one or two calls.",
        "bullets": [
            "`model.hds` exposes the `HeadsPlus` helper, while `model.all_heads` provides the normalized heads table.",
            "Use `model.cor(...)` for choropleths and `model.xsect(...)` or `plot_model_cross_section(...)` for section views.",
            "Result viewers can now pivot between time selectors like `per` and explicit `kstpkper` pairs more consistently.",
            "The same model object supports both interactive review and file-backed reopened runs.",
        ],
        "code": [
            "# pull normalized heads from the model",
            "heads = model.hds",
            "table = model.all_heads.reset_index()",
            "",
            "# map one stress period",
            "fig = model.cor(per=0, layer=0, type='hds').plot()",
            "",
            "# build a section view through selected cells",
            "section = model.xsect(per=0, cells=[10, 20, 30])",
            "mpl = mf.plot_model_cross_section(model, line=section.line)",
            "",
            "# reopened runs use the same surfaces",
            "head_slice = model.all_heads.query('layer == 0')",
        ],
        "api_notes": [
            "The newer results surface is meant to shorten the path from file-backed heads to maps and sections.",
            "Because reopened runs use the same surfaces, `HeadsPlus`-style workflows carry through after calibration too.",
        ],
    },
    {
        "title": "3. Targets & Observations",
        "subtitle": "Define targets once, then reuse them for review, FloPy observations, and calibration.",
        "accent": "navy",
        "narrative": "Targets are now meant to be reusable domain objects rather than a format hidden inside one notebook or one calibration script. Define them once, then feed them everywhere else.",
        "takeaway": "The same target definition should support comparison, plotting, FloPy observation files, and PEST.",
        "bullets": [
            "`HeadTargets` supports GIS points, cell dictionaries, `Series`, record tables, and `from_cells(...)` / `from_records(...)` helpers.",
            "`LakeStageTargets` supports lake-id mappings, `Series`, records, and plain lists with `times=...`.",
            "Attach them to the model with `model.targets.heads = ...` and `model.targets.lake_stage = ...`.",
            "Generate MF6 observation objects directly with `to_flopy_obs(...)` and `attach_flopy_obs(...)`.",
        ],
        "code": [
            "# attach reusable head targets",
            "model.targets.heads = mf.HeadTargets(",
            "    locations={'OBS_01': 101},",
            "    values=series,",
            "    time_column='per',",
            ")",
            "",
            "# lake stages can come from simple sequences too",
            "model.targets.lake_stage = mf.LakeStageTargets(",
            "    locations={'lake_main': 0},",
            "    values=[766.1, 766.3],",
            "    times=[0, 1],",
            ")",
            "",
            "obs_pkg = model.targets.heads.attach_flopy_obs(",
            "    filename='head_targets.obs',",
            "    csv_name='head_targets.csv',",
            ")",
            "",
            "compare = model.targets.heads.compare()",
            "stats = model.targets.heads.stats()",
        ],
        "api_notes": [
            "The plain constructors now accept lightweight dict and series inputs, while the classmethods remain useful for more structured cases.",
            "The same target object can stay attached to the model and still generate lower-level FloPy observation packages when needed.",
        ],
    },
    {
        "title": "4. Budgets & Package Results",
        "subtitle": "Use the improved budget wrappers and package explorers for tabular and map views.",
        "accent": "teal",
        "narrative": "Budget inspection should be as direct as head inspection. The newer budget surfaces wrap both listing-style summaries and package-level results so you can stay in one API while moving from totals to cells.",
        "takeaway": "Use the top-level budget helper for general accounting, then drop into package explorers for cell-by-cell detail.",
        "bullets": [
            "`model.bud()` returns the improved `Budget` wrapper, while `model.bud('rch')` or similar scopes the view to one package.",
            "Package explorers expose the unified grammar -- `get()`, `summary()`, `map()`, `plot()`, `mosaic()`, `animate()` -- consistently across LAK, SFR, and other namespaces.",
            "Budget viewers and result explorers now share more predictable period selectors and normalized result tables.",
            "Use the package namespaces when you want connection-aware or term-aware output rather than one flat budget table.",
        ],
        "code": [
            "# work from the top-level budget helper",
            "gwf_budget = model.bud()",
            "rch_budget = model.bud('rch')",
            "",
            "# listing-style summaries stay easy to reach",
            "summary = gwf_budget.budget()",
            "plot = rch_budget.plot_budget_obs(obs_dict)",
            "",
            "# package explorers go deeper",
            "lak_q = model.packages.lak.results.q.get(per=0)",
            "model.packages.sfr.results.q.profile.plot()",
            "model.packages.lak.results.q.budget.plot(per=0)",
            "model.packages.lak.budget.gwf.summary()",
        ],
        "api_notes": [
            "The modern budget workflow tries to make totals, tables, and package detail feel connected instead of like separate tools.",
            "When you need richer package semantics, the package explorers are usually the better surface than the raw listing budget.",
        ],
    },
    {
        "title": "5. Surface Water",
        "subtitle": "Lakes, streams, and combined exchange now have clearer result and validation surfaces.",
        "accent": "gold",
        "narrative": "Surface-water diagnosis is no longer split across a handful of isolated helpers. Lakes, streams, and combined exchange can now be inspected through a more consistent package/results surface, with validation close by.",
        "takeaway": "The surface-water API should make it easy to move from one package to a combined system view without changing conventions.",
        "bullets": [
            "LAK and SFR package explorers now expose normalized result tables and budget-term helper namespaces.",
            "Combined `surface_water` exchange maps use a shared sign convention and shared scaling so lakes and streams can be read together.",
            "`model.validate_surface_water()` and related validators are meant to catch configuration issues before they become interpretation problems.",
            "Grouped surface-water subplots carry the same ideas into scenario comparison.",
        ],
        "code": [
            "# inspect lake and stream results directly",
            "stage = model.packages.lak.results.stage.get()",
            "model.packages.sfr.results.q.map(per=0)",
            "",
            "# use the combined exchange view when needed",
            "model.packages.surface_water.results.q.map(per=0)",
            "",
            "# validate configuration before deeper review",
            "report = model.validate_surface_water()",
            "",
            "# package budget helpers stay available too",
            "model.packages.sfr.budget.flow_ja_face.get(per=0)",
        ],
        "api_notes": [
            "The package-level helper namespaces matter because lake-stage, stream exchange, and routing terms do not all live at the same physical level.",
            "The combined `surface_water` view is often the fastest way to spot where the system is gaining or losing water overall.",
        ],
    },
    {
        "title": "6. Model Groups",
        "subtitle": "Compare multiple runs through one lazy, aligned API.",
        "accent": "navy",
        "narrative": "Scenario comparison should not require custom merge logic every time. `ModelGroup` is there to align models, preserve a reference run, and expose grouped heads, grouped packages, and grouped surface-water plots from one surface.",
        "takeaway": "Use groups when comparison is the primary task, not as an afterthought after single-model analysis.",
        "bullets": [
            "`ModelGroup` accepts model objects or run/workspace paths and can lazily open file-backed runs for comparison.",
            "Use `group.hds.get()` and `group.hds.compare()` for aligned head tables and reference-based differences.",
            "Grouped package accessors now support side-by-side maps and plots for inputs and results.",
            "Use `shared_grid=True` when models share the same discretization so the grid view can be reused efficiently.",
        ],
        "code": [
            "# build a comparison group",
            "group = mf.ModelGroup(models, reference='base', shared_grid=True)",
            "",
            "# compare heads against the reference model",
            "head_diff = group.hds.compare(per=0)",
            "",
            "# build grouped package panels (unified grammar)",
            "group.packages.surface_water.results.q.mosaic(by='model', per=0)",
            "group.packages.lak.results.stage.plot()",
            "",
            "# summarize the participating runs",
            "group.summary()",
        ],
        "api_notes": [
            "The grouped API is most useful when you want one reference model and one or more alternatives to stay aligned automatically.",
            "Grouped surfaces are meant to reuse the same ideas as the single-model surfaces, not invent a parallel dialect.",
        ],
    },
    {
        "title": "7. Project & Run Management",
        "subtitle": "The same codebase also supports archive-style and project-managed workflows.",
        "accent": "teal",
        "narrative": "Not every useful model starts inside one notebook. Some are long-lived archives, some are organized around related scenarios, and some need richer artifact lineage. The project layer is there for that broader lifecycle.",
        "takeaway": "Use the lightweight exploration path for existing runs, and the catalog path when you want managed scenario workflows.",
        "bullets": [
            "Use `Project`, `Run`, and reusable specs when you want managed related runs and reproducible scenario lineages.",
            "Use `Project.discover_native_runs(...)` or `load_mf6_run(...)` when you want to inspect existing archives without rebuilding models.",
            "The long-term aim is a smooth path from model build to calibration to results review without leaving the same conceptual API.",
            "Notebook and script examples in `examples/mf6/` are the fastest way to see those patterns together.",
        ],
        "code": [
            "# lazy-open an existing MF6 archive",
            "loaded = mf.load_mf6_run(workspace)",
            "run = mf.Project(project_root).discover_native_runs(root, load=True)[0]",
            "",
            "# organize related runs when lineage matters",
            "project = mf.Project(project_root)",
            "records = explorer.records",
            "",
            "# reopen a completed calibration run",
            "run = mf.open_pest_run(pest_workspace)",
        ],
        "api_notes": [
            "Use the lazy loader when you want to inspect existing run folders without rebuilding anything first.",
            "Use the project/catalog path when scenario lineage and artifacts matter more than ad hoc experimentation.",
        ],
    },
    {
        "title": "8. Plotting & Review",
        "subtitle": "Targets, calibration plots, and model-bound review surfaces stay close together.",
        "accent": "gold",
        "narrative": "Plotting is easier to discover when it lives next to the thing being plotted. The newer model-bound target surfaces are meant to shorten the path from a target definition to maps, cross-plots, and time-series review.",
        "takeaway": "Use `model.targets...` for the common review calls, and let the lower-level plotting tools stay secondary.",
        "bullets": [
            "Use `model.targets.heads.plot.locations()` for target maps and `plot.obs_vs_sim()` / `plot.timeseries(...)` for review.",
            "The same targets can generate `CalibrationPlot` views without re-ingesting spreadsheets or shapefiles.",
            "This target layer is intentionally independent of PEST, even though the calibration workflow also consumes it.",
            "Lake-stage targets and head targets can both stay attached to the model while still generating FloPy observation packages.",
        ],
        "code": [
            "# compare target heads against one model",
            "compare = model.targets.heads.compare()",
            "stats = model.targets.heads.stats()",
            "",
            "# stay close to the target object while plotting",
            "model.targets.heads.plot.locations()",
            "model.targets.heads.plot.obs_vs_sim()",
            "",
            "# transient review stays one call away",
            "model.targets.heads.plot.timeseries('OBS_01')",
            "",
            "plot = mf.CalibrationPlot.from_targets(model.targets.heads)",
            "mf.CalibrationPlot.from_compare(compare)",
        ],
        "api_notes": [
            "The long-term role of `CalibrationPlot` is to consume the target system rather than act as a separate ingestion layer.",
            "If a review workflow still needs custom joins or one-off plotting glue, that is usually an API smell worth fixing.",
        ],
    },
    {
        "title": "9. Calibration / PEST",
        "subtitle": "The reusable pyEMU layer is designed to sit on top of the general model workflow.",
        "accent": "navy",
        "narrative": "Calibration should feel like an extension of the model workflow, not a separate universe. Parameter specs and observation specs should read like structured inputs, while saved metadata should make the resulting runs easy to reopen later.",
        "takeaway": "Keep target definitions general and let the PEST layer adapt them instead of redefining them.",
        "bullets": [
            "Create a `PestProject(model=..., workspace=...)` and add parameter specs such as `KPilotPointParameter` and `DrainConductanceParameter`.",
            "Use `HeadTargetObservationSpec(targets=model.targets.heads.targets)` so the calibration layer consumes the shared target object.",
            "Built workspaces save normalized metadata and target snapshots so completed runs can be reopened automatically later.",
            "Result review is meant to flow back into the normal model/review API rather than staying trapped in the PEST workspace.",
        ],
        "code": [
            "# build a reusable calibration workspace",
            "pest = mf.PestProject(",
            "    model=model,",
            "    name='calibration_run',",
            "    workspace=model.workspace / 'pest',",
            "    start_datetime='2024-01-01',",
            ")",
            "",
            "pest.add_parameter(mf.KPilotPointParameter(...))",
            "pest.add_observation(",
            "    mf.HeadTargetObservationSpec(",
            "        targets=model.targets.heads.targets,",
            "    )",
            ")",
            "",
            "pst = pest.build_pst('run.pst')",
            "pst.control_data.noptmax = 3",
        ],
        "api_notes": [
            "The calibration layer adapts model-side targets rather than redefining them, which keeps the mental model cleaner.",
            "Saved run metadata is what later allows completed PEST workspaces to reopen smoothly through the review API.",
        ],
    },
    {
        "title": "10. Reopen & Reuse",
        "subtitle": "Completed calibration workspaces should come back as normal reviewable model objects.",
        "accent": "teal",
        "narrative": "A calibrated run is most useful when it comes back as a normal model object. The post-run review path should feel like continuing the same analysis, not starting over from a folder full of artifacts.",
        "takeaway": "Reopened calibrated models should be first-class `SimulationBase`-style objects, not just piles of files on disk.",
        "bullets": [
            "Open a completed workspace with `mf.open_pest_run(...)` and load the baseline model, calibrated model, and saved targets.",
            "`load_calibrated_model()` reapplies final parameters and materializes fresh MF6 outputs before comparison.",
            "Use `review()`, `compare_head_targets()`, `plot_k()`, `plot_k_ratio()`, and `plot_well_timeseries(...)` for post-run diagnosis.",
            "This same reopen surface is what lets the broader project and notebook workflows stay smooth after long runs complete.",
        ],
        "code": [
            "# reopen a completed calibration workspace",
            "run = mf.open_pest_run(workspace)",
            "",
            "# load baseline and calibrated model states",
            "baseline = run.load_baseline_model()",
            "calibrated = run.load_calibrated_model()",
            "targets = run.load_head_targets()",
            "",
            "# collect common review products",
            "review = run.review()",
            "run.plot_k_ratio()",
            "run.compare_head_target_stats()",
        ],
        "api_notes": [
            "The reopened calibrated model is meant to behave like a normal `SimulationBase`-style object, not just a file bundle.",
            "The cleanest workflow is the one that keeps the same concepts alive from first build to final review.",
        ],
    },
]


def build_styles():
    sample = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "PamTitle",
            parent=sample["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=24,
            leading=28,
            textColor=PALETTE["ink"],
            spaceAfter=8,
        ),
        "subtitle": ParagraphStyle(
            "PamSubtitle",
            parent=sample["BodyText"],
            fontName="Helvetica",
            fontSize=10.5,
            leading=14,
            textColor=PALETTE["slate"],
        ),
        "body": ParagraphStyle(
            "PamBody",
            parent=sample["BodyText"],
            fontName="Helvetica",
            fontSize=10,
            leading=13,
            textColor=PALETTE["ink"],
        ),
        "body_compact": ParagraphStyle(
            "PamBodyCompact",
            parent=sample["BodyText"],
            fontName="Helvetica",
            fontSize=9.1,
            leading=11.2,
            textColor=PALETTE["ink"],
        ),
        "kicker": ParagraphStyle(
            "PamKicker",
            parent=sample["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=8.5,
            leading=11,
            textColor=PALETTE["teal"],
            spaceAfter=4,
        ),
        "narrative": ParagraphStyle(
            "PamNarrative",
            parent=sample["BodyText"],
            fontName="Helvetica",
            fontSize=9.7,
            leading=12.2,
            textColor=PALETTE["ink"],
        ),
        "takeaway": ParagraphStyle(
            "PamTakeaway",
            parent=sample["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=8.8,
            leading=10.8,
            textColor=PALETTE["slate"],
        ),
        "code": ParagraphStyle(
            "PamCode",
            parent=sample["Code"],
            fontName="Courier",
            fontSize=8.5,
            leading=11,
            textColor=PALETTE["ink"],
        ),
        "cover_tag": ParagraphStyle(
            "PamCoverTag",
            parent=sample["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=18,
            textColor=colors.white,
        ),
    }


def draw_paragraph(c: canvas.Canvas, text: str, style: ParagraphStyle, x: float, y_top: float, width: float):
    para = Paragraph(text, style)
    w, h = para.wrap(width, 1000)
    para.drawOn(c, x, y_top - h)
    return y_top - h, h


def draw_bullet_list(c: canvas.Canvas, items: list[str], styles, x: float, y_top: float, width: float):
    y = y_top
    bullet_w = 12
    for item in items:
        c.setFillColor(PALETTE["teal"])
        c.circle(x + 4, y - 7, 2.3, fill=1, stroke=0)
        y, _ = draw_paragraph(c, item, styles["body_compact"], x + bullet_w, y, width - bullet_w)
        y -= 5
    return y


def tokenize_code_line(line: str):
    """Tokenize one Python-ish code line into styled segments."""

    def classify_name(name: str, previous_sig: str | None, next_sig: str | None):
        if keyword.iskeyword(name):
            return PALETTE["code_keyword"], "Courier-Bold"
        if name == "mf":
            return PALETTE["code_module"], "Courier-Bold"
        if name.isupper():
            return PALETTE["code_number"], "Courier-Bold"
        if name and name[0].isupper():
            return PALETTE["code_class"], "Courier-Bold"
        if next_sig == "(":
            return PALETTE["code_func"], "Courier"
        return PALETTE["code_name"], "Courier"

    try:
        raw_tokens = [
            tok for tok in tokenize.generate_tokens(io.StringIO(line).readline)
            if tok.type not in {token_mod.ENDMARKER, token_mod.NEWLINE, token_mod.NL, token_mod.INDENT, token_mod.DEDENT}
        ]
    except tokenize.TokenError:
        raw_tokens = []
    if not raw_tokens:
        return [(line, PALETTE["code_text"], "Courier")]

    segments: list[tuple[str, colors.Color, str]] = []
    last_col = 0
    significant = [tok for tok in raw_tokens if tok.type not in {tokenize.INDENT, tokenize.DEDENT}]
    for idx, tok in enumerate(raw_tokens):
        start_col = tok.start[1]
        end_col = tok.end[1]
        if start_col > last_col:
            gap = line[last_col:start_col]
            if gap:
                segments.append((gap, PALETTE["code_text"], "Courier"))
        text = line[start_col:end_col]
        if tok.type == token_mod.NAME:
            sig_index = significant.index(tok) if tok in significant else -1
            prev_sig = significant[sig_index - 1].string if sig_index > 0 else None
            next_sig = significant[sig_index + 1].string if 0 <= sig_index < len(significant) - 1 else None
            color, font = classify_name(text, prev_sig, next_sig)
        elif tok.type == token_mod.STRING:
            color, font = PALETTE["code_string"], "Courier"
        elif tok.type == token_mod.NUMBER:
            color, font = PALETTE["code_number"], "Courier"
        elif tok.type == token_mod.COMMENT:
            color, font = PALETTE["code_comment"], "Courier-Oblique"
        else:
            color, font = PALETTE["code_operator"], "Courier"
        segments.append((text, color, font))
        last_col = end_col
    if last_col < len(line):
        tail = line[last_col:]
        if tail:
            segments.append((tail, PALETTE["code_text"], "Courier"))
    return segments


def segment_width(segment: tuple[str, colors.Color, str], font_size: float) -> float:
    text, _, font = segment
    return stringWidth(text, font, font_size)


def split_segment(segment: tuple[str, colors.Color, str], max_width: float, font_size: float):
    """Split one long segment to fit the requested width."""

    text, color, font = segment
    if segment_width(segment, font_size) <= max_width:
        return [segment]
    pieces: list[tuple[str, colors.Color, str]] = []
    current = ""
    for char in text:
        candidate = current + char
        if stringWidth(candidate, font, font_size) <= max_width or not current:
            current = candidate
            continue
        pieces.append((current, color, font))
        current = char
    if current:
        pieces.append((current, color, font))
    return pieces


def wrap_code_segments(lines: list[str], *, font_size: float, width: float):
    """Wrap tokenized code segments to fit inside the code card."""

    wrapped: list[list[tuple[str, colors.Color, str]]] = []
    for line in lines:
        segments = tokenize_code_line(line)
        current: list[tuple[str, colors.Color, str]] = []
        current_w = 0.0
        for segment in segments:
            seg_w = segment_width(segment, font_size)
            if current_w + seg_w <= width:
                current.append(segment)
                current_w += seg_w
                continue
            if seg_w > width:
                for piece in split_segment(segment, width if not current else width - current_w, font_size):
                    piece_w = segment_width(piece, font_size)
                    if current and current_w + piece_w > width:
                        wrapped.append(current)
                        current = [piece]
                        current_w = piece_w
                    else:
                        current.append(piece)
                        current_w += piece_w
            else:
                if current:
                    wrapped.append(current)
                current = [(segment[0].lstrip(), segment[1], segment[2])]
                current_w = segment_width(current[0], font_size)
        wrapped.append(current if current else [("", PALETTE["code_text"], "Courier")])
    return wrapped


def draw_code_block(c: canvas.Canvas, lines: list[str], x: float, y_top: float, width: float):
    font_size = 8.15
    line_h = 10.8
    padding = 10
    gutter_w = 24
    code_width = width - padding * 2 - gutter_w - 4
    wrapped_lines = wrap_code_segments(lines, font_size=font_size, width=code_width)
    height = padding * 2 + line_h * len(wrapped_lines) + 18
    c.setFillColor(PALETTE["code_bg"])
    c.roundRect(x, y_top - height, width, height, 10, fill=1, stroke=0)
    c.setStrokeColor(PALETTE["code_border"])
    c.roundRect(x, y_top - height, width, height, 10, fill=0, stroke=1)
    c.setFillColor(PALETTE["code_gutter"])
    c.roundRect(x, y_top - height, gutter_w + 12, height, 10, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#FF5F56"))
    c.circle(x + 14, y_top - 12, 3, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#FFBD2E"))
    c.circle(x + 24, y_top - 12, 3, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#27C93F"))
    c.circle(x + 34, y_top - 12, 3, fill=1, stroke=0)
    c.setFillColor(PALETTE["code_accent"])
    c.setFont("Helvetica-Bold", 7)
    c.drawRightString(x + width - 10, y_top - 15, "python")
    current_y = y_top - padding - 22
    line_no = 1
    for line_segments in wrapped_lines:
        c.setFillColor(colors.HexColor("#6E7681"))
        c.setFont("Courier", 7.5)
        c.drawRightString(x + gutter_w, current_y, str(line_no))
        current_x = x + gutter_w + 10
        for text, color, font in line_segments:
            if not text:
                continue
            c.setFillColor(color)
            c.setFont(font, font_size)
            c.drawString(current_x, current_y, text)
            current_x += stringWidth(text, font, font_size)
        current_y -= line_h
        line_no += 1
    return y_top - height


def measure_code_block_height(lines: list[str], width: float):
    """Return the rendered height for one code block."""

    font_size = 8.15
    line_h = 10.8
    padding = 10
    gutter_w = 24
    code_width = width - padding * 2 - gutter_w - 4
    wrapped_lines = wrap_code_segments(lines, font_size=font_size, width=code_width)
    return padding * 2 + line_h * len(wrapped_lines) + 18


def draw_header(c: canvas.Canvas, page, page_no: int, styles):
    width, height = letter
    accent = PALETTE[page["accent"]]
    c.setFillColor(accent)
    c.rect(0, height - 72, width, 72, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 22)
    c.drawString(48, height - 42, page["title"])
    c.setFont("Helvetica", 10)
    c.drawString(48, height - 58, page["subtitle"])
    c.setFillColor(PALETTE["sand"])
    c.rect(0, 0, width, 18, fill=1, stroke=0)
    c.setFillColor(PALETTE["slate"])
    c.setFont("Helvetica", 8)
    c.drawRightString(width - 36, 6, f"Page {page_no}")
    c.drawString(36, 6, "myflopy workflow pamphlet")


def draw_cover(c: canvas.Canvas, page, styles):
    width, height = letter
    c.setFillColor(PALETTE["navy"])
    c.rect(0, 0, width, height, fill=1, stroke=0)
    c.setFillColor(PALETTE["sand"])
    c.roundRect(40, 90, width - 80, height - 180, 20, fill=1, stroke=0)
    c.setFillColor(PALETTE["navy"])
    c.setFont("Helvetica-Bold", 30)
    c.drawString(68, height - 150, page["title"])
    c.setFont("Helvetica", 13)
    c.setFillColor(PALETTE["slate"])
    c.drawString(68, height - 175, page["subtitle"])

    c.setFillColor(PALETTE["teal"])
    c.roundRect(68, height - 265, width - 136, 64, 16, fill=1, stroke=0)
    y_top, _ = draw_paragraph(c, page["tagline"], styles["cover_tag"], 88, height - 220, width - 176)

    start_y = height - 325
    card_w = (width - 176) / 2
    card_h = 92
    positions = [(68, start_y), (68 + card_w + 20, start_y), (68, start_y - card_h - 16), (68 + card_w + 20, start_y - card_h - 16)]
    fills = [PALETTE["mist"], PALETTE["leaf"], PALETTE["sky"], PALETTE["rose"]]
    for (title, body), (x, y), fill in zip(page["cards"], positions, fills):
        c.setFillColor(fill)
        c.roundRect(x, y - card_h, card_w, card_h, 14, fill=1, stroke=0)
        c.setFillColor(PALETTE["ink"])
        c.setFont("Helvetica-Bold", 13)
        c.drawString(x + 14, y - 22, title)
        draw_paragraph(c, body, styles["body"], x + 14, y - 34, card_w - 28)

    c.setFillColor(PALETTE["gold"])
    c.setFont("Helvetica-Bold", 10)
    c.drawString(68, 64, "Suggested companion files")
    c.setFillColor(PALETTE["slate"])
    c.setFont("Helvetica", 9)
    c.drawString(68, 49, "docs/preferred_api.md")
    c.drawString(68, 37, "examples/mf6/notebooks/targets_api_quickstart.ipynb")


def draw_topic_page(c: canvas.Canvas, page, styles, page_no: int):
    width, height = letter
    draw_header(c, page, page_no, styles)
    content_left = 44
    content_right = width - 44
    top = height - 96
    content_width = content_right - content_left

    c.setFillColor(colors.white)
    c.rect(0, 18, width, height - 90, fill=1, stroke=0)
    c.setFillColor(PALETTE["mist"])
    c.circle(width - 80, height - 120, 46, fill=1, stroke=0)
    c.setFillColor(PALETTE["rose"] if page["accent"] == "gold" else PALETTE["leaf"])
    c.circle(width - 120, height - 155, 20, fill=1, stroke=0)

    # Measure the upper text stack first so the cards can be sized to fit
    # without overlapping the code section below.
    measure_y = top - 14
    para = Paragraph("Overview & flow", styles["kicker"])
    _w, h = para.wrap(content_width - 32, 1000)
    measure_y -= h
    para = Paragraph(page["narrative"], styles["narrative"])
    _w, h = para.wrap(content_width - 32, 1000)
    measure_y -= h + 2
    take_y = measure_y - 26
    para = Paragraph("Why it matters", styles["kicker"])
    _w, h = para.wrap(content_width - 56, 1000)
    take_y -= h
    para = Paragraph(page["takeaway"], styles["takeaway"])
    _w, h = para.wrap(content_width - 56, 1000)
    _take_bottom = take_y - h - 2
    text_end = measure_y - 90
    for item in page["bullets"]:
        para = Paragraph(item, styles["body_compact"])
        _w, h = para.wrap(content_width - 44, 1000)
        text_end -= h + 5
    notes = page.get("api_notes", [])
    if notes:
        text_end -= 2
        para = Paragraph("Reading notes", styles["kicker"])
        _w, h = para.wrap(content_width - 32, 1000)
        text_end -= h
        for item in notes:
            para = Paragraph(item, styles["body_compact"])
            _w, h = para.wrap(content_width - 44, 1000)
            text_end -= h + 5

    code_x = content_left + 16
    code_width = content_width - 32
    code_height = measure_code_block_height(page["code"], code_width)
    code_bottom_margin = 56
    code_top = max(code_bottom_margin + code_height + 14, text_end - 28)
    top_card_bottom = max(code_top + 8, text_end - 10)
    top_card_height = top - top_card_bottom - 12
    lower_card_top = code_top + 10
    lower_card_height = max(140, lower_card_top - 54)

    c.setFillColor(PALETTE["sand"])
    c.roundRect(content_left, top_card_bottom, content_width, top_card_height, 14, fill=1, stroke=0)
    c.setFillColor(colors.HexColor("#EEF3F7"))
    c.roundRect(content_left, 54, content_width, lower_card_height, 14, fill=1, stroke=0)

    # Redraw the text and code now that the dynamic card frames are in place.
    y = top - 14
    y, _ = draw_paragraph(c, "Overview & flow", styles["kicker"], content_left + 16, y, content_width - 32)
    y, _ = draw_paragraph(c, page["narrative"], styles["narrative"], content_left + 16, y - 2, content_width - 32)

    c.setFillColor(PALETTE["leaf"] if page["accent"] != "gold" else PALETTE["rose"])
    c.roundRect(content_left + 16, y - 64, content_width - 32, 52, 10, fill=1, stroke=0)
    y_take = y - 26
    y_take, _ = draw_paragraph(c, "Why it matters", styles["kicker"], content_left + 28, y_take, content_width - 56)
    draw_paragraph(c, page["takeaway"], styles["takeaway"], content_left + 28, y_take - 2, content_width - 56)

    text_x = content_left + 16
    text_y = y - 90
    c.setFillColor(PALETTE["ink"])
    c.setFont("Helvetica-Bold", 12)
    c.drawString(text_x, text_y, "Suggested flow")
    text_y = draw_bullet_list(c, page["bullets"], styles, text_x, text_y - 18, content_width - 32)
    if notes:
        text_y -= 2
        text_y, _ = draw_paragraph(c, "Reading notes", styles["kicker"], text_x, text_y, content_width - 32)
        draw_bullet_list(c, notes, styles, text_x, text_y - 4, content_width - 32)

    c.setFillColor(PALETTE["ink"])
    c.setFont("Helvetica-Bold", 12)
    c.drawString(code_x, code_top, "Key API surface")
    draw_code_block(c, page["code"], code_x, code_top - 14, code_width)


def build_pdf():
    styles = build_styles()
    c = canvas.Canvas(str(OUTPUT), pagesize=letter)
    draw_cover(c, PAGES[0], styles)
    c.showPage()
    for page_no, page in enumerate(PAGES[1:], start=2):
        draw_topic_page(c, page, styles, page_no)
        c.showPage()
    c.save()


if __name__ == "__main__":
    build_pdf()
    print(OUTPUT)
