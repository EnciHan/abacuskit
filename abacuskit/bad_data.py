"""Bad-data screening helpers for ABACUS and DeepMD datasets."""

from __future__ import annotations

import csv
import json
import math
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


@dataclass
class FrameRecord:
    system: Path
    set_name: str
    local_index: int
    global_index: int
    natoms: int
    energy: float | None
    energy_per_atom: float | None
    max_force: float | None


def unique_output_path(path: Path) -> Path:
    path = path.expanduser()
    if not path.exists():
        return path
    for idx in range(1, 1000):
        candidate = path.parent / f"{path.name}_{idx:03d}"
        if not candidate.exists():
            return candidate
    raise FileExistsError(f"cannot find unused output path for {path}")


def robust_z(values: Iterable[float | None]) -> np.ndarray:
    arr = np.array([np.nan if value is None else float(value) for value in values], dtype=float)
    finite = np.isfinite(arr)
    z = np.zeros(arr.shape, dtype=float)
    if finite.sum() < 3:
        return z
    vals = arr[finite]
    median = float(np.median(vals))
    mad = float(np.median(np.abs(vals - median)))
    if mad > 0:
        z[finite] = 0.67448975 * (vals - median) / mad
        return z
    std = float(np.std(vals))
    if std > 0:
        z[finite] = (vals - float(np.mean(vals))) / std
    return z


def _as_frame_array(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 0:
        return arr.reshape(1)
    return arr


def _frame_count(set_dir: Path) -> int:
    for name in ("energy.npy", "coord.npy", "force.npy", "box.npy", "virial.npy"):
        path = set_dir / name
        if path.is_file():
            return int(np.load(path, mmap_mode="r").shape[0])
    return 0


def _read_natoms(system: Path) -> int | None:
    type_raw = system / "type.raw"
    if type_raw.is_file():
        values = [line.strip() for line in type_raw.read_text().splitlines() if line.strip()]
        if values:
            return len(values)
    first_coord = next(iter(sorted(system.glob("set.*/coord.npy"))), None)
    if first_coord is not None:
        coord = np.load(first_coord, mmap_mode="r")
        if coord.ndim >= 2 and coord.shape[1] % 3 == 0:
            return int(coord.shape[1] // 3)
    return None


def is_deepmd_system(path: Path) -> bool:
    return path.is_dir() and (path / "type.raw").is_file() and any(path.glob("set.*/coord.npy"))


def find_deepmd_systems(paths: Iterable[Path]) -> list[Path]:
    systems: list[Path] = []
    seen: set[Path] = set()
    for raw_path in paths:
        path = raw_path.expanduser()
        candidates: list[Path] = []
        if is_deepmd_system(path):
            candidates.append(path)
        elif path.is_dir():
            for child in sorted(path.iterdir()):
                if is_deepmd_system(child):
                    candidates.append(child)
                elif child.is_dir():
                    candidates.extend(sorted(p for p in child.iterdir() if is_deepmd_system(p)))
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved not in seen:
                seen.add(resolved)
                systems.append(candidate)
    return systems


def load_deepmd_records(system: Path, global_start: int = 0) -> list[FrameRecord]:
    natoms = _read_natoms(system)
    if not natoms:
        raise ValueError(f"cannot infer natoms for DeepMD system {system}")
    records: list[FrameRecord] = []
    global_index = global_start
    for set_dir in sorted(system.glob("set.*")):
        if not set_dir.is_dir():
            continue
        nframes = _frame_count(set_dir)
        if nframes <= 0:
            continue
        energy_path = set_dir / "energy.npy"
        force_path = set_dir / "force.npy"
        energies = np.load(energy_path) if energy_path.is_file() else None
        forces = np.load(force_path) if force_path.is_file() else None
        for local_idx in range(nframes):
            energy = None
            energy_pa = None
            if energies is not None:
                energy = float(np.asarray(energies[local_idx]).reshape(-1)[0])
                energy_pa = energy / natoms
            max_force = None
            if forces is not None:
                force_frame = np.asarray(forces[local_idx], dtype=float)
                max_force = float(np.nanmax(np.abs(force_frame)))
            records.append(
                FrameRecord(
                    system=system,
                    set_name=set_dir.name,
                    local_index=local_idx,
                    global_index=global_index,
                    natoms=natoms,
                    energy=energy,
                    energy_per_atom=energy_pa,
                    max_force=max_force,
                )
            )
            global_index += 1
    return records


def _reasons_for_frame(
    record: FrameRecord,
    energy_z: float,
    force_z: float,
    max_force: float | None,
    max_abs_energy_per_atom: float | None,
    energy_sigma: float,
    force_sigma: float,
) -> list[str]:
    reasons: list[str] = []
    values = [record.energy_per_atom, record.max_force]
    if any(value is not None and not math.isfinite(value) for value in values):
        reasons.append("nonfinite_label")
    if max_force is not None and record.max_force is not None and record.max_force > max_force:
        reasons.append("max_force")
    if (
        max_abs_energy_per_atom is not None
        and record.energy_per_atom is not None
        and abs(record.energy_per_atom) > max_abs_energy_per_atom
    ):
        reasons.append("max_abs_energy_per_atom")
    if energy_sigma > 0 and record.energy_per_atom is not None and abs(energy_z) > energy_sigma:
        reasons.append("energy_outlier")
    if force_sigma > 0 and record.max_force is not None and abs(force_z) > force_sigma:
        reasons.append("force_outlier")
    return reasons


def filter_records_by_labels(
    records: list[FrameRecord],
    max_force: float | None,
    max_abs_energy_per_atom: float | None,
    energy_sigma: float,
    force_sigma: float,
) -> dict[int, list[str]]:
    energy_z = robust_z(record.energy_per_atom for record in records)
    force_z = robust_z(record.max_force for record in records)
    rejected: dict[int, list[str]] = {}
    for idx, record in enumerate(records):
        reasons = _reasons_for_frame(
            record,
            float(energy_z[idx]),
            float(force_z[idx]),
            max_force,
            max_abs_energy_per_atom,
            energy_sigma,
            force_sigma,
        )
        if reasons:
            rejected[record.global_index] = reasons
    return rejected


def copy_deepmd_system_filtered(system: Path, out: Path, rejected_global: set[int], global_start: int = 0) -> tuple[int, int]:
    if out.exists():
        raise FileExistsError(f"output system already exists: {out}")
    out.mkdir(parents=True)
    for item in system.iterdir():
        if item.is_file():
            shutil.copy2(item, out / item.name)

    kept = 0
    total = 0
    global_index = global_start
    for set_dir in sorted(system.glob("set.*")):
        if not set_dir.is_dir():
            continue
        nframes = _frame_count(set_dir)
        if nframes <= 0:
            continue
        keep_indices = [idx for idx in range(nframes) if global_index + idx not in rejected_global]
        total += nframes
        global_index += nframes
        if not keep_indices:
            continue
        keep_arr = np.array(keep_indices, dtype=int)
        out_set = out / set_dir.name
        out_set.mkdir()
        for npy in sorted(set_dir.glob("*.npy")):
            arr = _as_frame_array(np.load(npy))
            if arr.shape[0] == nframes:
                np.save(out_set / npy.name, arr[keep_arr])
            else:
                np.save(out_set / npy.name, arr)
        kept += len(keep_indices)
    return kept, total


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def _detail_file(root: Path, prefix: str, suffix: str) -> Path | None:
    candidate = root / f"detail.{prefix}.{suffix}.out"
    if candidate.is_file():
        return candidate
    matches = sorted(root.glob(f"**/detail.{prefix}.{suffix}.out"))
    return matches[0] if matches else None


def _load_detail_pair(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def find_detail_prefix(root: Path, requested: str) -> str | None:
    if requested != "auto":
        return requested
    for prefix in ("train", "valid"):
        if _detail_file(root, prefix, "e") or _detail_file(root, prefix, "e_peratom") or _detail_file(root, prefix, "f"):
            return prefix
    return None


def _detail_prefixes(root: Path, requested: str) -> list[str]:
    if requested != "auto":
        return [requested]
    prefixes = []
    for prefix in ("train", "valid"):
        if _detail_file(root, prefix, "e") or _detail_file(root, prefix, "e_peratom") or _detail_file(root, prefix, "f"):
            prefixes.append(prefix)
    return prefixes


def detail_rejections(
    records: list[FrameRecord],
    detail_root: Path,
    prefix: str,
    max_energy_error: float | None,
    max_force_error: float | None,
    energy_error_sigma: float,
    force_error_sigma: float,
) -> dict[int, list[str]]:
    rejected: dict[int, list[str]] = {}
    nframes = len(records)

    efile = _detail_file(detail_root, prefix, "e_peratom") or _detail_file(detail_root, prefix, "e")
    if efile is not None:
        edata = _load_detail_pair(efile)
        if edata.shape[1] >= 2 and edata.shape[0] == nframes:
            residual = edata[:, 1] - edata[:, 0]
            if not efile.name.endswith("e_peratom.out"):
                natoms = np.array([record.natoms for record in records], dtype=float)
                residual = residual / natoms
            residual_mev = residual * 1000.0
            z = robust_z(np.abs(residual_mev))
            for idx, value in enumerate(np.abs(residual_mev)):
                reasons = []
                if max_energy_error is not None and value > max_energy_error:
                    reasons.append("energy_error")
                if energy_error_sigma > 0 and abs(float(z[idx])) > energy_error_sigma:
                    reasons.append("energy_error_outlier")
                if reasons:
                    rejected.setdefault(records[idx].global_index, []).extend(reasons)

    ffile = _detail_file(detail_root, prefix, "f")
    if ffile is not None:
        fdata = _load_detail_pair(ffile)
        if fdata.shape[1] >= 6:
            expected_atoms = sum(record.natoms for record in records)
            if fdata.shape[0] == expected_atoms:
                atom_residual = np.linalg.norm(fdata[:, 3:6] - fdata[:, 0:3], axis=1)
                frame_errors: list[float] = []
                cursor = 0
                for record in records:
                    chunk = atom_residual[cursor : cursor + record.natoms]
                    cursor += record.natoms
                    frame_errors.append(float(np.max(chunk)) if chunk.size else 0.0)
                z = robust_z(frame_errors)
                for idx, value in enumerate(frame_errors):
                    reasons = []
                    if max_force_error is not None and value > max_force_error:
                        reasons.append("force_error")
                    if force_error_sigma > 0 and abs(float(z[idx])) > force_error_sigma:
                        reasons.append("force_error_outlier")
                    if reasons:
                        rejected.setdefault(records[idx].global_index, []).extend(reasons)
    return rejected


def _record_rows(records: list[FrameRecord], rejected: dict[int, list[str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        reasons = rejected.get(record.global_index, [])
        rows.append(
            {
                "system": str(record.system),
                "set": record.set_name,
                "local_index": record.local_index,
                "global_index": record.global_index,
                "natoms": record.natoms,
                "energy": record.energy,
                "energy_per_atom": record.energy_per_atom,
                "max_force": record.max_force,
                "keep": not reasons,
                "reasons": ";".join(sorted(set(reasons))),
            }
        )
    return rows


def screen_trained_deepmd(
    systems: list[Path],
    out: Path,
    detail_root: Path | None = None,
    detail_prefix: str = "auto",
    max_force: float | None = 50.0,
    max_abs_energy_per_atom: float | None = None,
    energy_sigma: float = 6.0,
    force_sigma: float = 6.0,
    max_energy_error: float | None = None,
    max_force_error: float | None = None,
    energy_error_sigma: float = 6.0,
    force_error_sigma: float = 6.0,
) -> dict[str, Any]:
    out = unique_output_path(out)
    out.mkdir(parents=True)

    records: list[FrameRecord] = []
    system_starts: dict[Path, int] = {}
    cursor = 0
    for system in systems:
        system_starts[system] = cursor
        loaded = load_deepmd_records(system, cursor)
        records.extend(loaded)
        cursor += len(loaded)
    if not records:
        raise ValueError("no DeepMD frames found")

    rejected = filter_records_by_labels(records, max_force, max_abs_energy_per_atom, energy_sigma, force_sigma)
    used_prefixes: list[str] = []
    if detail_root is not None:
        for used_prefix in _detail_prefixes(detail_root, detail_prefix):
            detail_found = False
            global_detail_rejected = detail_rejections(
                records,
                detail_root,
                used_prefix,
                max_energy_error,
                max_force_error,
                energy_error_sigma,
                force_error_sigma,
            )
            if global_detail_rejected:
                detail_found = True
            for idx, reasons in global_detail_rejected.items():
                rejected.setdefault(idx, []).extend(reasons)
            for system in systems:
                label = f"{system.parent.name}/{system.name}".lower()
                if used_prefix not in label:
                    continue
                subset = [record for record in records if record.system == system]
                if not subset or len(subset) == len(records):
                    continue
                subset_detail_rejected = detail_rejections(
                    subset,
                    detail_root,
                    used_prefix,
                    max_energy_error,
                    max_force_error,
                    energy_error_sigma,
                    force_error_sigma,
                )
                if subset_detail_rejected:
                    detail_found = True
                for idx, reasons in subset_detail_rejected.items():
                    rejected.setdefault(idx, []).extend(reasons)
            if detail_found:
                used_prefixes.append(used_prefix)

    clean_root = out / "cleaned"
    clean_root.mkdir()
    system_rows: list[dict[str, Any]] = []
    for system in systems:
        target = clean_root / system.name
        kept, total = copy_deepmd_system_filtered(system, target, set(rejected), system_starts[system])
        system_rows.append({"system": str(system), "cleaned": str(target), "frames": total, "kept": kept, "removed": total - kept})

    rows = _record_rows(records, rejected)
    fields = [
        "system",
        "set",
        "local_index",
        "global_index",
        "natoms",
        "energy",
        "energy_per_atom",
        "max_force",
        "keep",
        "reasons",
    ]
    _write_csv(out / "frames.csv", rows, fields)
    _write_csv(out / "rejected_frames.csv", [row for row in rows if not row["keep"]], fields)
    summary = {
        "mode": "trained_deepmd",
        "out": str(out),
        "detail_root": str(detail_root) if detail_root else None,
        "detail_prefix": ",".join(used_prefixes) if used_prefixes else None,
        "systems": system_rows,
        "frames": len(records),
        "kept": len(records) - len(rejected),
        "removed": len(rejected),
    }
    _write_json(out / "summary.json", summary)
    return summary


def _load_dpdata_system(job: Path, fmt: str):
    import dpdata

    return dpdata.LabeledSystem(str(job), fmt=fmt)


def _dpdata_frame_records(system: Any, job: Path, global_start: int) -> list[FrameRecord]:
    natoms = int(system.get_natoms())
    nframes = len(system)
    energies = system.data.get("energies")
    forces = system.data.get("forces")
    records: list[FrameRecord] = []
    for idx in range(nframes):
        energy = float(energies[idx]) if energies is not None else None
        energy_pa = energy / natoms if energy is not None else None
        max_force = None
        if forces is not None:
            max_force = float(np.nanmax(np.abs(forces[idx])))
        records.append(
            FrameRecord(
                system=job,
                set_name="dpdata",
                local_index=idx,
                global_index=global_start + idx,
                natoms=natoms,
                energy=energy,
                energy_per_atom=energy_pa,
                max_force=max_force,
            )
        )
    return records


def screen_raw_abacus(
    jobs: list[dict[str, Any]],
    out: Path,
    max_force: float | None = 50.0,
    max_abs_energy_per_atom: float | None = None,
    energy_sigma: float = 6.0,
    force_sigma: float = 6.0,
    set_size: int = 5000,
    keep_unconverged: bool = False,
) -> dict[str, Any]:
    out = unique_output_path(out)
    out.mkdir(parents=True)

    loaded: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    cursor = 0
    for item in jobs:
        job = Path(item["job"])
        status = item.get("status") or {}
        fmt = item["fmt"]
        if not keep_unconverged and not status.get("converged"):
            job_rows.append(
                {
                    "job": str(job),
                    "fmt": fmt,
                    "frames": 0,
                    "kept": 0,
                    "removed": 0,
                    "status": "rejected",
                    "reason": "not_converged",
                    "message": status.get("message", ""),
                }
            )
            continue
        try:
            system = _load_dpdata_system(job, fmt)
            records = _dpdata_frame_records(system, job, cursor)
        except Exception as exc:
            job_rows.append(
                {
                    "job": str(job),
                    "fmt": fmt,
                    "frames": 0,
                    "kept": 0,
                    "removed": 0,
                    "status": "rejected",
                    "reason": "load_failed",
                    "message": str(exc),
                }
            )
            continue
        if not records:
            job_rows.append(
                {
                    "job": str(job),
                    "fmt": fmt,
                    "frames": 0,
                    "kept": 0,
                    "removed": 0,
                    "status": "rejected",
                    "reason": "no_frames",
                    "message": "",
                }
            )
            continue
        loaded.append({"job": job, "fmt": fmt, "system": system, "records": records, "start": cursor})
        cursor += len(records)

    records = [record for item in loaded for record in item["records"]]
    rejected = filter_records_by_labels(records, max_force, max_abs_energy_per_atom, energy_sigma, force_sigma)

    deepmd_root = out / "deepmd"
    deepmd_root.mkdir()
    for item in loaded:
        system = item["system"]
        records_item = item["records"]
        keep_indices = [idx for idx, record in enumerate(records_item) if record.global_index not in rejected]
        target = deepmd_root / item["job"].name
        if keep_indices:
            sub_system = system.sub_system(keep_indices)
            sub_system.to("deepmd/npy", str(target), set_size=set_size)
        kept = len(keep_indices)
        total = len(records_item)
        job_rows.append(
            {
                "job": str(item["job"]),
                "fmt": item["fmt"],
                "frames": total,
                "kept": kept,
                "removed": total - kept,
                "status": "accepted" if kept else "rejected",
                "reason": "frame_filter" if kept < total else "",
                "message": "",
            }
        )

    frame_rows = _record_rows(records, rejected)
    fields = [
        "system",
        "set",
        "local_index",
        "global_index",
        "natoms",
        "energy",
        "energy_per_atom",
        "max_force",
        "keep",
        "reasons",
    ]
    _write_csv(out / "jobs.csv", job_rows, ["job", "fmt", "frames", "kept", "removed", "status", "reason", "message"])
    _write_csv(out / "frames.csv", frame_rows, fields)
    _write_csv(out / "rejected_frames.csv", [row for row in frame_rows if not row["keep"]], fields)
    summary = {
        "mode": "raw_abacus",
        "out": str(out),
        "jobs": len(jobs),
        "loaded_jobs": len(loaded),
        "frames": len(records),
        "kept": len(records) - len(rejected),
        "removed": len(rejected),
        "deepmd": str(deepmd_root),
    }
    _write_json(out / "summary.json", summary)
    return summary
