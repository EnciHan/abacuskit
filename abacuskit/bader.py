"""Bader charge post-processing helpers for ABACUS charge-density cubes."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Iterable

import numpy as np
from ase.data import chemical_symbols


@dataclass
class CubeAtom:
    index: int
    atomic_number: int
    symbol: str
    valence_electrons: float
    x: float
    y: float
    z: float


@dataclass
class AcfAtom:
    index: int
    x: float
    y: float
    z: float
    bader_electrons: float
    min_dist: float | None
    atomic_volume: float | None


@dataclass
class PreparedCube:
    charge_cube: Path
    source_cubes: list[Path]
    generated_total_cube: bool


def natural_key(text: str) -> tuple:
    import re

    parts = re.split(r"(\d+)", text)
    return tuple(int(part) if part.isdigit() else part.lower() for part in parts)


def _first_file(root: Path, names: Iterable[str], patterns: Iterable[str]) -> Path | None:
    for name in names:
        candidate = root / name
        if candidate.is_file():
            return candidate
    for pattern in patterns:
        found = sorted([path for path in root.glob(pattern) if path.is_file()], key=lambda path: natural_key(path.name))
        if found:
            return found[0]
    return None


def find_charge_cubes(root: Path) -> dict[str, Path]:
    """Find likely ABACUS charge-density cube files under an OUT.* directory."""
    root = Path(root)
    cubes: dict[str, Path] = {}
    total = _first_file(root, ["TOTAL_CHG.cube", "CHG.cube"], ["*TOTAL*CHG*.cube", "*CHG.cube"])
    spin1 = _first_file(root, ["SPIN1_CHG.cube"], ["*SPIN1*CHG*.cube"])
    spin2 = _first_file(root, ["SPIN2_CHG.cube"], ["*SPIN2*CHG*.cube"])
    if total:
        cubes["total"] = total
    if spin1:
        cubes["spin1"] = spin1
    if spin2:
        cubes["spin2"] = spin2
    return cubes


def read_cube_grid(path: Path) -> tuple[list[str], list[CubeAtom], np.ndarray]:
    lines = Path(path).read_text(errors="ignore").splitlines()
    if len(lines) < 6:
        raise ValueError(f"cube file too short: {path}")
    try:
        natoms = abs(int(float(lines[2].split()[0])))
        nx = int(float(lines[3].split()[0]))
        ny = int(float(lines[4].split()[0]))
        nz = int(float(lines[5].split()[0]))
    except (IndexError, ValueError) as exc:
        raise ValueError(f"cannot parse cube header: {path}") from exc

    start = 6 + natoms
    if len(lines) < start:
        raise ValueError(f"cube file misses atom lines: {path}")

    atoms: list[CubeAtom] = []
    for offset, line in enumerate(lines[6:start], start=1):
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(f"cannot parse cube atom line {offset}: {path}")
        atomic_number = int(float(fields[0]))
        symbol = chemical_symbols[atomic_number] if 0 <= atomic_number < len(chemical_symbols) else str(atomic_number)
        atoms.append(
            CubeAtom(
                index=offset,
                atomic_number=atomic_number,
                symbol=symbol,
                valence_electrons=float(fields[1]),
                x=float(fields[2]),
                y=float(fields[3]),
                z=float(fields[4]),
            )
        )

    values: list[float] = []
    for line in lines[start:]:
        for part in line.split():
            try:
                values.append(float(part))
            except ValueError:
                pass
    need = nx * ny * nz
    if len(values) < need:
        raise ValueError(f"cube file has {len(values)} values, expected {need}: {path}")
    return lines[:start], atoms, np.array(values[:need], dtype=float).reshape((nx, ny, nz))


def write_cube_grid(path: Path, header: list[str], values: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    flat = values.reshape(-1)
    lines = list(header)
    for start in range(0, len(flat), 6):
        chunk = flat[start : start + 6]
        lines.append(" ".join(f"{value: .10e}" for value in chunk))
    path.write_text("\n".join(lines) + "\n")


def _headers_match(path_a: Path, header_a: list[str], path_b: Path, header_b: list[str], shape_a, shape_b) -> None:
    if shape_a != shape_b:
        raise ValueError(f"cube grid shapes differ: {path_a} {shape_a} vs {path_b} {shape_b}")
    if len(header_a) != len(header_b):
        raise ValueError(f"cube headers have different atom counts: {path_a} vs {path_b}")
    if header_a[2:] != header_b[2:]:
        raise ValueError(f"cube geometry/atom headers differ: {path_a} vs {path_b}")


def prepare_charge_cube(path: Path, work_dir: Path, cube: Path | None = None, total_cube: Path | None = None) -> PreparedCube:
    """Return a total charge-density cube suitable for Bader analysis."""
    if cube is not None:
        cube = Path(cube).expanduser()
        if not cube.is_file():
            raise FileNotFoundError(f"cannot find charge-density cube: {cube}")
        return PreparedCube(charge_cube=cube, source_cubes=[cube], generated_total_cube=False)

    root = Path(path).expanduser()
    if root.is_file():
        return PreparedCube(charge_cube=root, source_cubes=[root], generated_total_cube=False)
    if not root.is_dir():
        raise FileNotFoundError(f"cannot find ABACUS output directory or cube file: {root}")

    cubes = find_charge_cubes(root)
    spin1 = cubes.get("spin1")
    spin2 = cubes.get("spin2")
    if spin1 and spin2:
        header1, _, values1 = read_cube_grid(spin1)
        header2, _, values2 = read_cube_grid(spin2)
        _headers_match(spin1, header1, spin2, header2, values1.shape, values2.shape)
        out = Path(total_cube).expanduser() if total_cube else Path(work_dir).expanduser() / "TOTAL_CHG.cube"
        header = list(header1)
        if len(header) >= 2:
            header[1] = f"Total charge density generated by abacuskit from {spin1.name} + {spin2.name}"
        write_cube_grid(out, header, values1 + values2)
        return PreparedCube(charge_cube=out, source_cubes=[spin1, spin2], generated_total_cube=True)

    total = cubes.get("total") or spin1
    if total:
        return PreparedCube(charge_cube=total, source_cubes=[total], generated_total_cube=False)
    raise FileNotFoundError(f"cannot find CHG/SPIN*_CHG cube files under {root}")


def resolve_bader_executable(bader: str | Path | None = None) -> str:
    requested = str(bader) if bader else os.environ.get("ABACUSKIT_BADER", "bader")
    if os.sep in requested or (os.altsep and os.altsep in requested):
        path = Path(requested).expanduser()
        if path.is_file():
            return str(path)
        raise FileNotFoundError(f"cannot find bader executable: {path}")
    found = shutil.which(requested)
    if found:
        return found
    raise FileNotFoundError(
        "cannot find bader executable; install Henkelman bader, put it in PATH, "
        "set ABACUSKIT_BADER, or pass --bader /path/to/bader"
    )


def run_bader_program(
    charge_cube: Path,
    work_dir: Path,
    bader: str | Path | None = None,
    reference_cube: Path | None = None,
    extra_args: list[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    executable = resolve_bader_executable(bader)
    work_dir = Path(work_dir).expanduser()
    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [executable, str(Path(charge_cube).resolve())]
    if reference_cube is not None:
        ref = Path(reference_cube).expanduser()
        if not ref.is_file():
            raise FileNotFoundError(f"cannot find Bader reference cube: {ref}")
        cmd.extend(["-ref", str(ref.resolve())])
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, cwd=work_dir, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        details = "\n".join(part for part in [result.stdout.strip(), result.stderr.strip()] if part)
        raise RuntimeError(f"bader failed with exit code {result.returncode}: {' '.join(cmd)}\n{details}")
    return result


def parse_acf(path: Path) -> tuple[list[AcfAtom], dict[str, float]]:
    atoms: list[AcfAtom] = []
    summary: dict[str, float] = {}
    for raw in Path(path).read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or set(line) <= {"-"}:
            continue
        fields = line.split()
        if fields[0].isdigit() and len(fields) >= 5:
            min_dist = float(fields[5]) if len(fields) > 5 else None
            atomic_volume = float(fields[6]) if len(fields) > 6 else None
            atoms.append(
                AcfAtom(
                    index=int(fields[0]),
                    x=float(fields[1]),
                    y=float(fields[2]),
                    z=float(fields[3]),
                    bader_electrons=float(fields[4]),
                    min_dist=min_dist,
                    atomic_volume=atomic_volume,
                )
            )
            continue
        if ":" in line:
            key, value = line.split(":", 1)
            try:
                summary[key.strip().lower().replace(" ", "_")] = float(value.split()[0])
            except (IndexError, ValueError):
                pass
    if not atoms:
        raise ValueError(f"cannot parse atom rows from ACF.dat: {path}")
    return atoms, summary


def combine_bader_rows(cube_atoms: list[CubeAtom], acf_atoms: list[AcfAtom]) -> list[dict[str, object]]:
    if len(cube_atoms) != len(acf_atoms):
        raise ValueError(f"cube atom count ({len(cube_atoms)}) differs from ACF atom count ({len(acf_atoms)})")
    rows: list[dict[str, object]] = []
    for cube_atom, acf_atom in zip(cube_atoms, acf_atoms):
        if cube_atom.index != acf_atom.index:
            raise ValueError(f"atom index mismatch between cube and ACF: {cube_atom.index} vs {acf_atom.index}")
        charge = cube_atom.valence_electrons - acf_atom.bader_electrons
        rows.append(
            {
                "atom_index": cube_atom.index,
                "symbol": cube_atom.symbol,
                "atomic_number": cube_atom.atomic_number,
                "valence_electrons": cube_atom.valence_electrons,
                "bader_electrons": acf_atom.bader_electrons,
                "charge": charge,
                "x": cube_atom.x,
                "y": cube_atom.y,
                "z": cube_atom.z,
                "min_dist": acf_atom.min_dist,
                "atomic_volume": acf_atom.atomic_volume,
            }
        )
    return rows


def run_bader_analysis(
    path: Path,
    work_dir: Path,
    cube: Path | None = None,
    total_cube: Path | None = None,
    bader: str | Path | None = None,
    reference_cube: Path | None = None,
    extra_args: list[str] | None = None,
) -> dict[str, object]:
    prepared = prepare_charge_cube(path=path, cube=cube, work_dir=work_dir, total_cube=total_cube)
    completed = run_bader_program(
        charge_cube=prepared.charge_cube,
        work_dir=work_dir,
        bader=bader,
        reference_cube=reference_cube,
        extra_args=extra_args,
    )
    acf_path = Path(work_dir).expanduser() / "ACF.dat"
    if not acf_path.is_file():
        raise FileNotFoundError(f"bader finished but ACF.dat was not written under {work_dir}")
    _, cube_atoms, _ = read_cube_grid(prepared.charge_cube)
    acf_atoms, summary = parse_acf(acf_path)
    rows = combine_bader_rows(cube_atoms, acf_atoms)
    return {
        "rows": rows,
        "summary": summary,
        "files": {
            "charge_cube": str(prepared.charge_cube),
            "source_cubes": [str(path) for path in prepared.source_cubes],
            "acf": str(acf_path),
            "work_dir": str(Path(work_dir).expanduser()),
        },
        "generated_total_cube": prepared.generated_total_cube,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def write_bader_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "atom_index",
        "symbol",
        "atomic_number",
        "valence_electrons",
        "bader_electrons",
        "charge",
        "x",
        "y",
        "z",
        "min_dist",
        "atomic_volume",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_bader_json(path: Path, metadata: dict[str, object]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2) + "\n")
