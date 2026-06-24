"""Journal-style Matplotlib defaults used by abacuskit plotting commands."""

from __future__ import annotations

from pathlib import Path

FIGSIZE_MM = {
    "single": (85, 65),
    "single_square": (85, 85),
    "single_tall": (85, 110),
    "double": (178, 110),
    "double_low": (178, 75),
    "double_tall": (178, 135),
    "triple_panel": (178, 65),
}

JOURNAL_DPI = 300
JOURNAL_SAVE_DPI = 600
JOURNAL_PAD_INCHES = 0.03
JOURNAL_AXES_LINEWIDTH = 1.2
JOURNAL_TICK_WIDTH = 1.0
EFERMI_RELATIVE_LABEL = r"$\mathbf{E}-\mathbf{E}_{\mathbf{F}}$ (eV)"


def mm_to_inch(mm):
    """Convert millimetres to inches; accepts a scalar or an iterable."""
    if isinstance(mm, (str, bytes)):
        raise TypeError("mm_to_inch expects a number or iterable of numbers")
    try:
        return tuple(float(value) / 25.4 for value in mm)
    except TypeError:
        return float(mm) / 25.4


def set_journal_style() -> None:
    """Apply compact publication-oriented Matplotlib rcParams."""
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.titlesize": 9,
            "axes.linewidth": JOURNAL_AXES_LINEWIDTH,
            "lines.linewidth": 1.2,
            "lines.markersize": 4,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": JOURNAL_TICK_WIDTH,
            "ytick.major.width": JOURNAL_TICK_WIDTH,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "xtick.minor.size": 2,
            "ytick.minor.size": 2,
            "legend.frameon": False,
            "figure.dpi": JOURNAL_DPI,
            "savefig.dpi": JOURNAL_SAVE_DPI,
            "savefig.bbox": "tight",
            "savefig.pad_inches": JOURNAL_PAD_INCHES,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def apply_export_axis_style(fig) -> None:
    """Apply final axis styling to every axes before writing a figure."""
    for ax in fig.axes:
        ax.grid(False, which="both", axis="both")
        ax.tick_params(width=JOURNAL_TICK_WIDTH)
        for spine in ax.spines.values():
            spine.set_linewidth(JOURNAL_AXES_LINEWIDTH)
        ax.xaxis.label.set_fontweight("bold")
        ax.yaxis.label.set_fontweight("bold")
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight("bold")


def get_figsize(kind: str) -> tuple[float, float]:
    """Return a journal figure size in inches for a named plot kind."""
    try:
        return mm_to_inch(FIGSIZE_MM[kind])
    except KeyError as exc:
        known = ", ".join(sorted(FIGSIZE_MM))
        raise ValueError(f"unknown figure size kind {kind!r}; choose from {known}") from exc


def save_journal_figure(fig, out_png, export_pdf: bool = True, dpi: int | None = None, bbox_tight: bool = True):
    """Save the existing PNG output and optionally a same-stem PDF."""
    import matplotlib as mpl

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs = {}
    if bbox_tight:
        save_kwargs.update({"bbox_inches": "tight", "pad_inches": JOURNAL_PAD_INCHES})
    else:
        save_kwargs.update({"bbox_inches": None, "pad_inches": 0.0})
    if dpi is not None:
        save_kwargs["dpi"] = dpi
    rc_bbox = "tight" if bbox_tight else None
    apply_export_axis_style(fig)
    with mpl.rc_context({"savefig.bbox": rc_bbox}):
        fig.savefig(out_path, **save_kwargs)
    written = [out_path]
    if export_pdf and out_path.suffix.lower() != ".pdf":
        pdf_path = out_path.with_suffix(".pdf")
        with mpl.rc_context({"savefig.bbox": rc_bbox}):
            if bbox_tight:
                fig.savefig(pdf_path, bbox_inches="tight", pad_inches=JOURNAL_PAD_INCHES)
            else:
                fig.savefig(pdf_path, bbox_inches=None, pad_inches=0.0)
        written.append(pdf_path)
    return written


def style_colorbar(colorbar, label: str | None = None, direction: str = "in") -> None:
    """Use consistent colorbar typography without changing plotted data."""
    colorbar.ax.tick_params(labelsize=8, width=0.8, length=3, direction=direction)
    if label is not None:
        colorbar.set_label(label, fontsize=9)
    elif colorbar.ax.get_ylabel():
        colorbar.ax.yaxis.label.set_size(9)


def _format_map_tick(value: float, _position: int) -> str:
    abs_value = abs(value)
    if abs_value >= 10:
        text = f"{value:.0f}"
    else:
        text = f"{value:.1f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return "0" if text in {"-0", "-0.0", ""} else text


def style_map_axes(ax, matched_ticks: bool = False, tick_count: int = 5) -> None:
    """Make ticks readable on filled color maps."""
    if matched_ticks:
        from matplotlib.ticker import FuncFormatter, LinearLocator

        ax.xaxis.set_major_locator(LinearLocator(tick_count))
        ax.yaxis.set_major_locator(LinearLocator(tick_count))
        ax.xaxis.set_major_formatter(FuncFormatter(_format_map_tick))
        ax.yaxis.set_major_formatter(FuncFormatter(_format_map_tick))
    ax.tick_params(direction="out", length=3.5, width=0.8)
    ax.grid(False)


def add_square_map_axes(
    fig,
    left: float = 0.16,
    bottom: float = 0.17,
    height: float = 0.66,
    pad: float = 0.035,
    colorbar_width: float = 0.035,
):
    """Add a physically square map axes and a same-height colorbar axes."""
    fig_width, fig_height = fig.get_size_inches()
    width = height * fig_height / fig_width
    ax = fig.add_axes([left, bottom, width, height])
    cax = fig.add_axes([left + width + pad, bottom, colorbar_width, height])
    ax.set_box_aspect(1)
    return ax, cax


def add_panel_label(ax, label: str, x: float = 0.02, y: float = 0.98) -> None:
    """Place a bold panel label inside the upper-left corner of an axes."""
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        fontweight="bold",
    )


def style_grid(ax, enabled: bool = True) -> None:
    ax.grid(False, which="both", axis="both")


def style_legend(ax, **kwargs):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    legend_kwargs = {"frameon": False, "fontsize": 8, "loc": "best"}
    legend_kwargs.update(kwargs)
    return ax.legend(**legend_kwargs)
