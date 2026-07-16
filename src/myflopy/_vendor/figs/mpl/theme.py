from dataclasses import dataclass
from contextlib import contextmanager
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

@dataclass(frozen=True)
class Theme:
    name: str
    rc: dict

    @contextmanager
    def context(self):
        with plt.rc_context(self.rc):
            yield

def apply_date_axis(ax, major="year", minor="month", fmt="%Y-%m"):
    if major == "year":
        ax.xaxis.set_major_locator(mdates.YearLocator())
    elif major == "month":
        ax.xaxis.set_major_locator(mdates.MonthLocator())
    elif major == "day":
        ax.xaxis.set_major_locator(mdates.DayLocator())

    if minor == "month":
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
    elif minor == "day":
        ax.xaxis.set_minor_locator(mdates.DayLocator())

    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))

REPORT = Theme("report", {
    "font.family": "sans-serif",
    "font.sans-serif": ["Calibri"],
    "font.size": 9,
    "figure.figsize": (6.5, 4),
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
