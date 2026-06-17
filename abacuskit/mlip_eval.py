"""DeepMD/MLIP evaluation plots for ABACUS reference data."""

from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

EV_PER_A3_TO_GPA = 160.21766208
DEFAULT_COMPONENTS = ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")


@dataclass
class MLIPEvalConfig:
    root: Path
    out: Path
    prefix: str = "valid"
    quantity: str = "all"
    natoms: int | None = None
    data_dir: Path | None = None
    title: str | None = None
    dpi: int = 300
    fmt: str = "png"
    outlier_sigma: float = 4.0
    top_outliers: int = 20


@dataclass
class QuantityData:
    key: str
    title: str
    ref_label: str
    pred_label: str
    axis_unit: str
    metric_unit: str
    ref: np.ndarray
    pred: np.ndarray
    components: list[str] | None = None

    @property
    def residual(self) -> np.ndarray:
        return self.pred - self.ref

    @property
    def count(self) -> int:
        return int(self.ref.size)


def _loadtxt(path: Path) -> np.ndarray:
    try:
        data = np.loadtxt(path)
    except OSError as exc:
        raise FileNotFoundError(path) from exc
    if data.size == 0:
        raise ValueError(f"{path} is empty")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _candidate_data_dirs(root: Path, detail_file: Path | None, data_dir: Path | None) -> list[Path]:
    candidates: list[Path] = []
    if data_dir is not None:
        candidates.append(data_dir)
    if detail_file is not None:
        candidates.extend([detail_file.parent / "data", detail_file.parent.parent / "data"])
    candidates.extend([root / "data", root.parent / "data"])

    result: list[Path] = []
    seen: set[Path] = set()
    for item in candidates:
        expanded = item.expanduser()
        try:
            key = expanded.resolve()
        except OSError:
            key = expanded.absolute()
        if key in seen:
            continue
        seen.add(key)
        result.append(expanded)
    return result


def _infer_natoms(root: Path, detail_file: Path | None, data_dir: Path | None) -> int | None:
    for base in _candidate_data_dirs(root, detail_file, data_dir):
        for raw in [base / "type.raw", base / "valid" / "type.raw", base / "train" / "type.raw"]:
            if raw.is_file():
                values = [line.strip() for line in raw.read_text().splitlines() if line.strip()]
                if values:
                    return len(values)
        for coord in [
            base / "set.000" / "coord.npy",
            base / "valid" / "set.000" / "coord.npy",
            base / "train" / "set.000" / "coord.npy",
        ]:
            if coord.is_file():
                arr = np.load(coord, mmap_mode="r")
                if arr.ndim == 3:
                    return int(arr.shape[1])
                if arr.ndim == 2 and arr.shape[1] % 3 == 0:
                    return int(arr.shape[1] // 3)
    return None


def _find_detail_file(root: Path, prefix: str, suffix: str) -> Path | None:
    direct = root / f"detail.{prefix}.{suffix}.out"
    if direct.is_file():
        return direct
    matches = sorted(root.glob(f"**/detail.{prefix}.{suffix}.out"))
    return matches[0] if matches else None


def _load_energy(root: Path, prefix: str, natoms: int | None, data_dir: Path | None) -> QuantityData | None:
    per_atom_file = _find_detail_file(root, prefix, "e_peratom")
    if per_atom_file is not None:
        data = _loadtxt(per_atom_file)
        if data.shape[1] >= 2:
            ref = data[:, 0]
            pred = data[:, 1]
            return QuantityData(
                key="energy",
                title="Energy",
                ref_label="DFT energy",
                pred_label="DP energy",
                axis_unit="eV/atom",
                metric_unit="meV/atom",
                ref=ref,
                pred=pred,
            )

    energy_file = _find_detail_file(root, prefix, "e")
    if energy_file is None:
        return None
    atoms = natoms or _infer_natoms(root, energy_file, data_dir)
    if not atoms:
        raise ValueError("cannot infer atom count for energy per-atom conversion; pass --natoms")
    data = _loadtxt(energy_file)
    if data.shape[1] < 2:
        raise ValueError(f"{energy_file} must contain reference and predicted energy columns")
    return QuantityData(
        key="energy",
        title="Energy",
        ref_label="DFT energy",
        pred_label="DP energy",
        axis_unit="eV/atom",
        metric_unit="meV/atom",
        ref=data[:, 0] / atoms,
        pred=data[:, 1] / atoms,
    )


def _load_force(root: Path, prefix: str) -> QuantityData | None:
    force_file = _find_detail_file(root, prefix, "f")
    if force_file is None:
        return None
    data = _loadtxt(force_file)
    if data.shape[1] < 6:
        raise ValueError(f"{force_file} must contain data_fx,data_fy,data_fz,pred_fx,pred_fy,pred_fz")
    ref = data[:, 0:3].reshape(-1)
    pred = data[:, 3:6].reshape(-1)
    components = [axis for _ in range(data.shape[0]) for axis in ("x", "y", "z")]
    return QuantityData(
        key="force",
        title="Force",
        ref_label="DFT force",
        pred_label="DP force",
        axis_unit="eV/Å",
        metric_unit="meV/Å",
        ref=ref,
        pred=pred,
        components=components,
    )


def _load_boxes(root: Path, detail_file: Path | None, data_dir: Path | None, frames: int) -> np.ndarray | None:
    for base in _candidate_data_dirs(root, detail_file, data_dir):
        for box in [
            base / "set.000" / "box.npy",
            base / "valid" / "set.000" / "box.npy",
            base / "train" / "set.000" / "box.npy",
        ]:
            if not box.is_file():
                continue
            arr = np.load(box)
            if arr.ndim == 3 and arr.shape[1:] == (3, 3):
                boxes = arr
            elif arr.ndim == 2 and arr.shape[1] == 9:
                boxes = arr.reshape(-1, 3, 3)
            else:
                continue
            if boxes.shape[0] >= frames:
                return boxes[:frames]
    return None


def _load_stress(root: Path, prefix: str, data_dir: Path | None) -> QuantityData | None:
    virial_file = _find_detail_file(root, prefix, "v")
    if virial_file is None:
        return None
    data = _loadtxt(virial_file)
    if data.shape[1] < 18:
        raise ValueError(f"{virial_file} must contain 9 reference and 9 predicted virial columns")
    ref_virial = data[:, 0:9]
    pred_virial = data[:, 9:18]
    if np.nanmax(np.abs(ref_virial)) <= 1.0e-12:
        return None

    boxes = _load_boxes(root, virial_file, data_dir, data.shape[0])
    if boxes is None:
        return QuantityData(
            key="virial",
            title="Virial",
            ref_label="DFT virial",
            pred_label="DP virial",
            axis_unit="eV",
            metric_unit="eV",
            ref=ref_virial.reshape(-1),
            pred=pred_virial.reshape(-1),
            components=[comp for _ in range(data.shape[0]) for comp in DEFAULT_COMPONENTS],
        )

    volumes = np.abs(np.linalg.det(boxes)).reshape(-1, 1)
    if np.any(volumes <= 0):
        raise ValueError("box volumes must be positive for virial-to-stress conversion")
    ref_stress = ref_virial / volumes * EV_PER_A3_TO_GPA
    pred_stress = pred_virial / volumes * EV_PER_A3_TO_GPA
    return QuantityData(
        key="stress",
        title="Stress",
        ref_label="DFT stress",
        pred_label="DP stress",
        axis_unit="GPa",
        metric_unit="GPa",
        ref=ref_stress.reshape(-1),
        pred=pred_stress.reshape(-1),
        components=[comp for _ in range(data.shape[0]) for comp in DEFAULT_COMPONENTS],
    )


def _metric_scale(quantity: QuantityData) -> float:
    return 1000.0 if quantity.metric_unit.startswith("meV") else 1.0


def _rmse(residual: np.ndarray) -> float:
    return float(np.sqrt(np.mean(residual**2)))


def _mae(residual: np.ndarray) -> float:
    return float(np.mean(np.abs(residual)))


def _r2(ref: np.ndarray, pred: np.ndarray) -> float:
    ss_res = float(np.sum((pred - ref) ** 2))
    ss_tot = float(np.sum((ref - np.mean(ref)) ** 2))
    if ss_tot <= 0.0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def _metrics(quantity: QuantityData) -> dict[str, float | int | str]:
    residual = quantity.residual
    scale = _metric_scale(quantity)
    return {
        "quantity": quantity.key,
        "points": quantity.count,
        "unit": quantity.metric_unit,
        "r2": _r2(quantity.ref, quantity.pred),
        "mae": _mae(residual) * scale,
        "rmse": _rmse(residual) * scale,
        "residual_mean": float(np.mean(residual) * scale),
        "residual_std": float(np.std(residual) * scale),
    }


def _axis_limits(ref: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    finite = np.concatenate([ref[np.isfinite(ref)], pred[np.isfinite(pred)]])
    if finite.size == 0:
        return -1.0, 1.0
    lo = float(np.min(finite))
    hi = float(np.max(finite))
    if math.isclose(lo, hi):
        pad = max(abs(lo) * 0.05, 1.0e-6)
    else:
        pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


def _hist_range(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0
    lo = float(np.min(finite))
    hi = float(np.max(finite))
    if math.isclose(lo, hi):
        pad = max(abs(lo) * 0.05, 1.0e-6)
    else:
        pad = 0.02 * (hi - lo)
    return lo - pad, hi + pad


def _panel(fig, spec, quantity: QuantityData, color: str, letter: str | None = None) -> None:
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    sub = spec.subgridspec(2, 2, height_ratios=(0.23, 1.0), width_ratios=(1.0, 0.23), hspace=0.03, wspace=0.03)
    ax_top = fig.add_subplot(sub[0, 0])
    ax = fig.add_subplot(sub[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(sub[1, 1], sharey=ax)

    ref = quantity.ref
    pred = quantity.pred
    residual = quantity.residual
    lo, hi = _axis_limits(ref, pred)

    ax.scatter(ref, pred, s=12 if quantity.key == "energy" else 7, color=color, alpha=0.45, linewidths=0)
    ax.plot([lo, hi], [lo, hi], color="0.5", lw=1.5, ls="--")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{quantity.ref_label} ({quantity.axis_unit})")
    ax.set_ylabel(f"{quantity.pred_label} ({quantity.axis_unit})")
    ax.tick_params(direction="in")
    ax.grid(True, color="0.88", lw=0.5)

    metrics = _metrics(quantity)
    ax.text(
        0.05,
        0.95,
        "\n".join(
            [
                f"$R^2$ = {metrics['r2']:.4f}" if np.isfinite(metrics["r2"]) else "$R^2$ = n/a",
                f"MAE = {metrics['mae']:.2f} {quantity.metric_unit}",
                f"RMSE = {metrics['rmse']:.2f} {quantity.metric_unit}",
                f"N = {metrics['points']}",
            ]
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
    )
    if letter:
        ax.text(-0.15, 1.18, f"({letter})", transform=ax.transAxes, fontsize=13, fontweight="bold")

    bins = min(60, max(12, int(np.sqrt(ref.size))))
    ax_top.hist(ref, bins=bins, range=_hist_range(ref), color=color, alpha=0.45)
    ax_top.axis("off")
    ax_right.hist(pred, bins=bins, range=_hist_range(pred), orientation="horizontal", color=color, alpha=0.45)
    ax_right.axis("off")

    inset = inset_axes(ax, width="32%", height="28%", loc="lower right", borderpad=1.2)
    res_scale = _metric_scale(quantity)
    residual_plot = residual * res_scale
    inset.hist(
        residual_plot,
        bins=min(50, max(10, int(np.sqrt(residual_plot.size)))),
        range=_hist_range(residual_plot),
        color=color,
        alpha=0.55,
    )
    inset.axvline(0.0, color="black", lw=0.9, ls="--")
    inset.set_title("Residual", fontsize=8, pad=1)
    inset.tick_params(labelsize=7, direction="in")
    for spine in inset.spines.values():
        spine.set_linewidth(0.8)


def _write_figure(path: Path, quantities: list[QuantityData], title: str | None, dpi: int) -> Path:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"energy": "#1f77b4", "force": "#2ca02c", "stress": "#ff7f0e", "virial": "#ff7f0e"}
    fig = plt.figure(figsize=(5.2 * len(quantities), 5.2), dpi=dpi)
    grid = fig.add_gridspec(1, len(quantities), wspace=0.34)
    for idx, quantity in enumerate(quantities):
        _panel(fig, grid[0, idx], quantity, colors.get(quantity.key, f"C{idx}"), chr(ord("a") + idx))
    if title:
        fig.suptitle(title, fontsize=13)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_single_figure(path: Path, quantity: QuantityData, title: str | None, dpi: int) -> Path:
    single_title = title or f"{quantity.title} parity"
    return _write_figure(path, [quantity], single_title, dpi)


def _outlier_rows(quantity: QuantityData, sigma: float, top_n: int) -> list[dict[str, object]]:
    residual = quantity.residual
    scaled = residual * _metric_scale(quantity)
    std = float(np.std(scaled))
    abs_res = np.abs(scaled)
    if std > 0.0 and sigma > 0.0:
        mask = abs_res >= sigma * std
    else:
        mask = np.zeros(abs_res.shape, dtype=bool)
    order = np.argsort(abs_res)[::-1]
    selected: list[int] = []
    for idx in order:
        if len(selected) < top_n or mask[idx]:
            selected.append(int(idx))
        if len(selected) >= top_n and (not mask.any() or all(not mask[i] for i in order[len(selected) :])):
            break

    rows: list[dict[str, object]] = []
    for idx in selected:
        z = float(abs_res[idx] / std) if std > 0 else float("nan")
        rows.append(
            {
                "quantity": quantity.key,
                "index": idx,
                "component": quantity.components[idx] if quantity.components and idx < len(quantity.components) else "",
                "ref": float(quantity.ref[idx]),
                "pred": float(quantity.pred[idx]),
                "residual": float(scaled[idx]),
                "abs_residual": float(abs_res[idx]),
                "residual_unit": quantity.metric_unit,
                "z_score": z,
                "outlier": bool(mask[idx]),
            }
        )
    return rows


def _write_summary(out: Path, quantities: list[QuantityData]) -> tuple[Path, Path]:
    csv_path = out / "summary.csv"
    json_path = out / "summary.json"
    rows = [_metrics(q) for q in quantities]
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["quantity", "points", "unit", "r2", "mae", "rmse", "residual_mean", "residual_std"])
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2) + "\n")
    return csv_path, json_path


def _write_outliers(out: Path, quantities: list[QuantityData], sigma: float, top_n: int) -> Path:
    path = out / "outliers.csv"
    rows: list[dict[str, object]] = []
    for quantity in quantities:
        rows.extend(_outlier_rows(quantity, sigma, top_n))
    fields = ["quantity", "index", "component", "ref", "pred", "residual", "abs_residual", "residual_unit", "z_score", "outlier"]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _filter_quantities(quantities: Iterable[QuantityData], choice: str) -> list[QuantityData]:
    all_quantities = list(quantities)
    if choice == "all":
        return all_quantities
    aliases = {"virial": {"virial", "stress"}, "stress": {"virial", "stress"}}
    wanted = aliases.get(choice, {choice})
    return [quantity for quantity in all_quantities if quantity.key in wanted]


def load_mlip_eval_data(config: MLIPEvalConfig) -> list[QuantityData]:
    root = config.root.expanduser()
    quantities: list[QuantityData] = []
    energy = _load_energy(root, config.prefix, config.natoms, config.data_dir)
    force = _load_force(root, config.prefix)
    stress = _load_stress(root, config.prefix, config.data_dir)
    for item in (energy, force, stress):
        if item is not None:
            quantities.append(item)
    selected = _filter_quantities(quantities, config.quantity)
    if not selected:
        raise ValueError(f"no {config.quantity!r} evaluation data found under {root}")
    return selected


def run_mlip_eval(config: MLIPEvalConfig) -> dict[str, object]:
    out = config.out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    quantities = load_mlip_eval_data(config)

    figure_paths: list[Path] = []
    if config.quantity == "all":
        figure_paths.append(_write_figure(out / f"{config.prefix}_mlip_eval_overview.{config.fmt}", quantities, config.title, config.dpi))

    for quantity in quantities:
        subdir = out / quantity.key
        subdir.mkdir(parents=True, exist_ok=True)
        figure_paths.append(
            _write_single_figure(
                subdir / f"{config.prefix}_{quantity.key}_parity.{config.fmt}",
                quantity,
                f"{quantity.title} parity",
                config.dpi,
            )
        )

    summary_csv, summary_json = _write_summary(out, quantities)
    outliers_csv = _write_outliers(out, quantities, config.outlier_sigma, config.top_outliers)
    return {
        "figures": [str(path) for path in figure_paths],
        "summary_csv": str(summary_csv),
        "summary_json": str(summary_json),
        "outliers_csv": str(outliers_csv),
        "quantities": [q.key for q in quantities],
    }
