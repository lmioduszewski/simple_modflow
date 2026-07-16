import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import seaborn as sns


def get_mplfig(
        title=None,
        xlabel=None,
        ylabel=None,
        major_color="#C0C0C0",
        minor_color="#E0E0E0",
        figsize=(11, 8.5),
        font_sanserif="Calibri",
        font_size=9,
        sns_theme="whitegrid",
):
    """
    Creates a Matplotlib figure and axis with customized appearance and styling.

    This function sets up a Matplotlib figure and axis with specified titles, labels,
    gridline appearances, and font properties. It configures the style using seaborn
    and updates several Matplotlib parameters to ensure consistent visual appearance.
    Axis gridlines for both major and minor ticks are styled separately, and several
    layout options are pre-configured for easy visualization.

    :param title: Title for the chart to be displayed above the plot
    :type title: str, optional
    :param xlabel: Label text for the x-axis
    :type xlabel: str, optional
    :param ylabel: Label text for the y-axis
    :type ylabel: str, optional
    :param major_color: Color for major gridlines and ticks in hex format
    :type major_color: str, optional
    :param minor_color: Color for minor gridlines and ticks in hex format
    :type minor_color: str, optional
    :param figsize: Tuple defining the figure size (width, height) in inches
    :type figsize: tuple, optional
    :param font_sanserif: Font family to be used for sans-serif font
    :type font_sanserif: str, optional
    :param font_size: Font size for the text elements in the figure
    :type font_size: int, optional
    :param sns_theme: seaborn theme to be used for styling the figure
    :type sns_theme: str, optional
    :return: A tuple containing the Matplotlib figure and axis objects
    :rtype: tuple
    """

    sns.set_context("paper")
    sns.set_theme(style=sns_theme)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font_sanserif],
        "font.size": font_size,
        "figure.figsize": figsize,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        'axes.linewidth': 1,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        'grid.linewidth': 0.5,
    })

    fig, ax = plt.subplots(constrained_layout=True)

    # --- Titles and labels ---
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Force tick marks on
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(axis="both", which="both", bottom=True, left=True, top=False, right=False)

    ax.tick_params(axis='both', which='major', labelsize=8, length=4, width=0.6, direction='out', color=major_color)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='minor', length=2, width=0.4, direction='out', color=minor_color, )
    ax.grid(True, which='major', axis='both', linewidth=0.6, color=major_color)
    ax.grid(True, which='minor', axis='both', linewidth=0.4, linestyle=':', alpha=1, color=minor_color)

    return fig, ax


def set_axis_scale(
    obj,
    x_units_per_inch,
    y_units_per_inch,
    ax=None,
    x_anchor="center",
    y_anchor="center",
    x_value=None,
    y_value=None,
    set_major_ticks=True,
    show_major_grid=True,
    tick_inches=1.0,
):
    """
    Set axis limits so the current axes box matches a desired engineering scale,
    while leaving the figure size and axes position unchanged.

    You can anchor each axis independently at its min, max, or center.

    example usage:

        set_axis_scale(ax, 2000, 100, y_anchor='max', y_value=750).savefig('lake_sawyer_xs.png', dpi=300)

    Parameters
    ----------
    obj : matplotlib.axes.Axes or matplotlib.figure.Figure
        Target axes, or figure containing the target axes.
    x_units_per_inch : float
        Data units per inch in x.
    y_units_per_inch : float
        Data units per inch in y.
    ax : matplotlib.axes.Axes, optional
        If obj is a Figure, which axes to use. Defaults to first axes.
    x_anchor : {"min", "max", "center"}, default "center"
        Which x reference point to keep fixed.
    y_anchor : {"min", "max", "center"}, default "center"
        Which y reference point to keep fixed.
    x_value : float, optional
        Value to hold fixed for x_anchor.
        If omitted, the current visible min/max/center is used.
    y_value : float, optional
        Value to hold fixed for y_anchor.
        If omitted, the current visible min/max/center is used.
    set_major_ticks : bool, default True
        Whether to set major tick spacing from the requested scale.
    show_major_grid : bool, default True
        Whether to turn on major gridlines.
    tick_inches : float, default 1.0
        Physical inches between major ticks/gridlines.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The same figure, after updating limits/ticks/grid.
    """
    if isinstance(obj, Axes):
        ax = obj
        fig = ax.figure
    elif isinstance(obj, Figure):
        fig = obj
        if ax is None:
            if not fig.axes:
                raise ValueError("Figure has no axes.")
            ax = fig.axes[0]
    else:
        raise TypeError("obj must be a matplotlib Axes or Figure.")

    if x_units_per_inch <= 0 or y_units_per_inch <= 0:
        raise ValueError("x_units_per_inch and y_units_per_inch must be positive.")
    if tick_inches <= 0:
        raise ValueError("tick_inches must be positive.")
    if x_anchor not in {"min", "max", "center"}:
        raise ValueError("x_anchor must be 'min', 'max', or 'center'.")
    if y_anchor not in {"min", "max", "center"}:
        raise ValueError("y_anchor must be 'min', 'max', or 'center'.")

    # Ensure final layout is applied before measuring axes size
    fig.canvas.draw()

    # Axes size in inches
    pos = ax.get_position()
    fig_w, fig_h = fig.get_size_inches()
    ax_w_in = pos.width * fig_w
    ax_h_in = pos.height * fig_h

    # Span implied by requested scale
    new_x_span = ax_w_in * x_units_per_inch
    new_y_span = ax_h_in * y_units_per_inch

    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()

    x_reversed = cur_xlim[1] < cur_xlim[0]
    y_reversed = cur_ylim[1] < cur_ylim[0]

    x0, x1 = sorted(cur_xlim)
    y0, y1 = sorted(cur_ylim)

    cur_x_center = 0.5 * (x0 + x1)
    cur_y_center = 0.5 * (y0 + y1)

    # Choose fixed reference values
    if x_value is None:
        if x_anchor == "min":
            x_value = x0
        elif x_anchor == "max":
            x_value = x1
        else:
            x_value = cur_x_center

    if y_value is None:
        if y_anchor == "min":
            y_value = y0
        elif y_anchor == "max":
            y_value = y1
        else:
            y_value = cur_y_center

    # Compute new limits from fixed reference
    if x_anchor == "min":
        new_xlim = (x_value, x_value + new_x_span)
    elif x_anchor == "max":
        new_xlim = (x_value - new_x_span, x_value)
    else:  # center
        new_xlim = (x_value - new_x_span / 2, x_value + new_x_span / 2)

    if y_anchor == "min":
        new_ylim = (y_value, y_value + new_y_span)
    elif y_anchor == "max":
        new_ylim = (y_value - new_y_span, y_value)
    else:  # center
        new_ylim = (y_value - new_y_span / 2, y_value + new_y_span / 2)

    # Restore inverted axes if needed
    if x_reversed:
        new_xlim = new_xlim[::-1]
    if y_reversed:
        new_ylim = new_ylim[::-1]

    ax.set_xlim(new_xlim)
    ax.set_ylim(new_ylim)

    if set_major_ticks:
        x_major = x_units_per_inch * tick_inches
        y_major = y_units_per_inch * tick_inches
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
        ax.yaxis.set_major_locator(MultipleLocator(y_major))

    if show_major_grid:
        ax.grid(True, which="major")

    return fig
