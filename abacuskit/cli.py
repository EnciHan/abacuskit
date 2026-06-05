#!/usr/bin/env python3
"""abacuskit: integrated ABACUS + DeepMD workflow toolkit.

The tool intentionally keeps every generated file plain and editable:
STRU/INPUT for ABACUS, shell scripts for running jobs, and JSON for DeepMD.

Author: Han Enci, Zhong Lisheng, Yu Yutong, Xu Mengting, Chen Jingyuan, Xi'an University of Technology.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

try:
    from . import __affiliation__, __author__, __version__
except ImportError:
    __version__ = "v1.1.1"
    __author__ = "Han Enci, Zhong Lisheng, Yu Yutong, Xu Mengting, Chen Jingyuan"
    __affiliation__ = "Xi'an University of Technology"

BOHR_PER_ANGSTROM = 1.88972612546


def env_path(name: str, default: str | Path) -> Path:
    return Path(os.environ.get(name, str(default))).expanduser()


DEFAULT_PSEUDO_DIR = env_path("ABACUSKIT_PSEUDO_DIR", "pseudopotentials")
DEFAULT_ORBITAL_DIRS = {
    "efficiency": env_path("ABACUSKIT_ORBITAL_EFFICIENCY_DIR", "orbitals/efficiency"),
    "precision": env_path("ABACUSKIT_ORBITAL_PRECISION_DIR", "orbitals/precision"),
}
DEFAULT_ABACUS_ROOT = env_path("ABACUSKIT_ABACUS_ROOT", "~/apps/abacus")
DEFAULT_ABACUS_ENV = env_path(
    "ABACUSKIT_ABACUS_ENV", DEFAULT_ABACUS_ROOT / "default" / "toolchain" / "abacus_env.sh"
)
DEFAULT_ABACUS_LABEL = os.environ.get("ABACUSKIT_ABACUS_LABEL", DEFAULT_ABACUS_ENV.parent.parent.name)
DEFAULT_DEEPMD_PYTHON = env_path("ABACUSKIT_DEEPMD_PYTHON", sys.executable)
DEFAULT_DP = env_path("ABACUSKIT_DP", "dp")
ORBITAL_LABEL_TO_L = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
L_TO_ORBITAL_LABEL = {v: k for k, v in ORBITAL_LABEL_TO_L.items()}
DFTU_KEYS = (
    "dft_plus_u",
    "orbital_corr",
    "hubbard_u",
    "yukawa_potential",
    "yukawa_lambda",
    "uramping",
    "omc",
    "onsite_radius",
)
DFTU_MIXING_KEYS = ("mixing_restart", "mixing_dmr", "uramping")
TERMINAL_LOGO = [
    r"  ___   ____    ___    ____  _   _  ____  _  __  ___  _____ ",
    r" / _ \ | __ )  / _ \  / ___|| | | |/ ___|| |/ / |_ _||_   _|",
    r"| /_\ ||  _ \ | /_\ || |    | | | |\___ \| ' /   | |   | |  ",
    r"|  _  || |_) ||  _  || |___ | |_| | ___) | . \   | |   | |  ",
    r"|_| |_||____/ |_| |_| \____| \___/ |____/|_|\_\ |___|  |_|  ",
]


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def terminal_color(text: str, code: str) -> str:
    if not sys.stdout.isatty() or os.environ.get("NO_COLOR"):
        return text
    return f"\033[{code}m{text}\033[0m"


def print_terminal_logo() -> None:
    print()
    for line in TERMINAL_LOGO:
        print(terminal_color(line, "96;1"))


def natural_key(text: str) -> tuple:
    parts = re.split(r"(\d+)", text)
    return tuple(int(p) if p.isdigit() else p.lower() for p in parts)


def unique_symbols(atoms) -> list[str]:
    symbols = []
    for sym in atoms.get_chemical_symbols():
        if sym not in symbols:
            symbols.append(sym)
    return symbols


def element_from_filename(path: Path) -> str | None:
    match = re.match(r"([A-Z][a-z]?)", path.name)
    if not match:
        return None
    return match.group(1)


def score_pseudo(path: Path, element: str) -> tuple[int, int, str]:
    name = path.name
    score = 0
    if name == f"{element}.upf":
        score += 100
    if name == f"{element}.UPF":
        score += 95
    if "ONCV" in name:
        score += 40
    if "PBE" in name.upper():
        score += 20
    if "FR" in name:
        score -= 5
    if "core" in name.lower():
        score -= 15
    return (-score, len(name), name)


def score_orbital(path: Path, element: str) -> tuple[int, int, str]:
    name = path.name
    score = 0
    if name.startswith(f"{element}_"):
        score += 100
    if "100Ry" in name:
        score += 15
    if "150Ry" in name:
        score += 10
    radius = re.search(r"_(\d+)au_", name)
    if radius:
        score += int(radius.group(1))
    # Prefer compact bases in the efficiency library, richer bases in precision.
    basis_tokens = re.findall(r"(\d+)[spdfg]", name)
    if basis_tokens:
        score -= sum(int(x) for x in basis_tokens)
    return (-score, len(name), name)


def index_library(directory: Path, suffixes: tuple[str, ...]) -> dict[str, list[Path]]:
    if not directory.is_dir():
        die(f"library directory not found: {directory}")
    by_element: dict[str, list[Path]] = defaultdict(list)
    for path in directory.iterdir():
        if path.is_file() and path.name.lower().endswith(suffixes):
            elem = element_from_filename(path)
            if elem:
                by_element[elem].append(path)
    return dict(by_element)


def parse_element_paths(items: list[str] | None, suffix: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for item in items or []:
        if "=" not in item:
            die(f"bad element file override {item!r}, expected Element=/path/to/file.{suffix}")
        sym, value = item.split("=", 1)
        path = Path(value).expanduser()
        if not path.is_file():
            die(f"override file not found for {sym}: {path}")
        if not path.name.lower().endswith(f".{suffix}"):
            die(f"override file for {sym} is not a .{suffix} file: {path}")
        result[sym.strip()] = path
    return result


def parse_element_quality(items: list[str] | None) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            die("bad --element-orbital-quality value, expected Element=efficiency or Element=precision")
        sym, quality = item.split("=", 1)
        quality = quality.strip().lower()
        if quality not in DEFAULT_ORBITAL_DIRS:
            die(f"unknown orbital quality {quality!r}; choose one of {sorted(DEFAULT_ORBITAL_DIRS)}")
        result[sym.strip()] = quality
    return result


def orbital_dir_for_element(
    element: str,
    global_orbital_dir: Path,
    element_quality: dict[str, str],
) -> Path:
    if element in element_quality:
        return DEFAULT_ORBITAL_DIRS[element_quality[element]]
    return global_orbital_dir


def choose_library_files(
    symbols: Iterable[str],
    pseudo_dir: Path,
    orbital_dir: Path,
    basis_type: str,
    pseudo_overrides: dict[str, Path] | None = None,
    orbital_overrides: dict[str, Path] | None = None,
    element_orbital_quality: dict[str, str] | None = None,
) -> tuple[dict[str, Path], dict[str, Path]]:
    pseudo_index = index_library(pseudo_dir, (".upf",))
    orbital_indexes: dict[Path, dict[str, list[Path]]] = {}
    pseudo_overrides = pseudo_overrides or {}
    orbital_overrides = orbital_overrides or {}
    element_orbital_quality = element_orbital_quality or {}

    pseudos: dict[str, Path] = {}
    orbitals: dict[str, Path] = {}
    missing = []
    for sym in symbols:
        if sym in pseudo_overrides:
            pseudos[sym] = pseudo_overrides[sym]
        else:
            pp_candidates = [p for p in pseudo_index.get(sym, []) if p.name.startswith(sym)]
            if not pp_candidates:
                missing.append(f"pseudo for {sym}")
            else:
                pseudos[sym] = sorted(pp_candidates, key=lambda p: score_pseudo(p, sym))[0]

        if basis_type == "lcao":
            if sym in orbital_overrides:
                orbitals[sym] = orbital_overrides[sym]
            else:
                elem_orbital_dir = orbital_dir_for_element(sym, orbital_dir, element_orbital_quality)
                if elem_orbital_dir not in orbital_indexes:
                    orbital_indexes[elem_orbital_dir] = index_library(elem_orbital_dir, (".orb",))
                orbital_index = orbital_indexes[elem_orbital_dir]
                orb_candidates = [p for p in orbital_index.get(sym, []) if p.name.startswith(f"{sym}_")]
                if not orb_candidates:
                    missing.append(f"orbital for {sym} in {elem_orbital_dir}")
                else:
                    orbitals[sym] = sorted(orb_candidates, key=lambda p: score_orbital(p, sym))[0]

    if missing:
        die("missing library files: " + ", ".join(missing))
    return pseudos, orbitals


def write_stru(
    atoms,
    output: Path,
    pseudo_dir: Path,
    orbital_dir: Path,
    basis_type: str = "lcao",
    magnetism: dict[str, float] | None = None,
    move_flags: list[tuple[int, int, int]] | None = None,
    pseudo_overrides: dict[str, Path] | None = None,
    orbital_overrides: dict[str, Path] | None = None,
    element_orbital_quality: dict[str, str] | None = None,
) -> dict:
    symbols = unique_symbols(atoms)
    pseudos, orbitals = choose_library_files(
        symbols,
        pseudo_dir,
        orbital_dir,
        basis_type,
        pseudo_overrides=pseudo_overrides,
        orbital_overrides=orbital_overrides,
        element_orbital_quality=element_orbital_quality,
    )
    magnetism = magnetism or {}
    if move_flags is None:
        move_flags = [(1, 1, 1)] * len(atoms)
    if len(move_flags) != len(atoms):
        die(f"move flag count {len(move_flags)} does not match atom count {len(atoms)}")

    output.parent.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[int]] = {sym: [] for sym in symbols}
    for idx, sym in enumerate(atoms.get_chemical_symbols()):
        grouped[sym].append(idx)

    lines: list[str] = ["ATOMIC_SPECIES"]
    for sym in symbols:
        mass = atomic_masses[atomic_numbers[sym]]
        lines.append(f"{sym} {mass:.6f} {pseudos[sym].name}")

    if basis_type == "lcao":
        lines += ["", "NUMERICAL_ORBITAL"]
        lines += [orbitals[sym].name for sym in symbols]

    lines += [
        "",
        "LATTICE_CONSTANT",
        f"{BOHR_PER_ANGSTROM:.11f}",
        "",
        "LATTICE_VECTORS",
    ]
    for vec in atoms.cell.array:
        lines.append(f"{vec[0]:18.10f} {vec[1]:18.10f} {vec[2]:18.10f}")

    lines += ["", "ATOMIC_POSITIONS", "Cartesian"]
    positions = atoms.get_positions()
    for sym in symbols:
        lines += ["", sym, f"{magnetism.get(sym, 0.0):.6f}", str(len(grouped[sym]))]
        for idx in grouped[sym]:
            x, y, z = positions[idx]
            mx, my, mz = move_flags[idx]
            lines.append(
                f"{x:18.10f} {y:18.10f} {z:18.10f} {mx} {my} {mz}"
            )

    output.write_text("\n".join(lines) + "\n")
    return {
        "symbols": symbols,
        "pseudo_dir": str(pseudo_dir),
        "orbital_dir": str(orbital_dir),
        "element_orbital_quality": element_orbital_quality or {},
        "pseudos": {k: str(v) for k, v in pseudos.items()},
        "orbitals": {k: str(v) for k, v in orbitals.items()},
        "magnetism": magnetism,
    }


def parse_magnetism(items: list[str] | None) -> dict[str, float]:
    result: dict[str, float] = {}
    for item in items or []:
        if "=" not in item:
            die(f"bad --mag value {item!r}, expected Element=value")
        sym, val = item.split("=", 1)
        result[sym] = float(val)
    return result


def direction_mask(dirs: str) -> tuple[int, int, int]:
    axes = dirs.lower().replace(",", "").replace(" ", "")
    if not axes or any(ch not in "xyz" for ch in axes):
        die(f"bad fixed directions {dirs!r}; use any combination of x, y, z, for example z or xy")
    mask = [1, 1, 1]
    for ch in axes:
        mask["xyz".index(ch)] = 0
    return tuple(mask)


def parse_index_selector(selector: str, natoms: int) -> list[int]:
    result: set[int] = set()
    for part in selector.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            left, right = part.split("-", 1)
            start = int(left)
            stop = int(right)
            if start > stop:
                start, stop = stop, start
            result.update(range(start - 1, stop))
        else:
            result.add(int(part) - 1)
    bad = [i + 1 for i in result if i < 0 or i >= natoms]
    if bad:
        die(f"atom index out of range: {bad}")
    return sorted(result)


def apply_mask(old: tuple[int, int, int], mask: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(min(a, b) for a, b in zip(old, mask))


def build_move_flags(
    atoms,
    fix_element: list[str] | None = None,
    fix_index: list[str] | None = None,
    fix_below: list[str] | None = None,
    fix_above: list[str] | None = None,
) -> list[tuple[int, int, int]]:
    flags = [(1, 1, 1)] * len(atoms)
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    for item in fix_element or []:
        if "=" not in item:
            die("bad --fix-element value, expected Element=dirs, for example Ni=z")
        sym, dirs = item.split("=", 1)
        mask = direction_mask(dirs)
        for idx, atom_sym in enumerate(symbols):
            if atom_sym == sym.strip():
                flags[idx] = apply_mask(flags[idx], mask)

    for item in fix_index or []:
        if "=" not in item:
            die("bad --fix-index value, expected 1,2,5-9=dirs")
        selector, dirs = item.split("=", 1)
        mask = direction_mask(dirs)
        for idx in parse_index_selector(selector, len(atoms)):
            flags[idx] = apply_mask(flags[idx], mask)

    axis_to_col = {"x": 0, "y": 1, "z": 2}
    for values, relation in ((fix_below or [], "<"), (fix_above or [], ">")):
        for item in values:
            if "=" not in item or ":" not in item:
                flag = "--fix-below" if relation == "<" else "--fix-above"
                die(f"bad {flag} value, expected axis=value:dirs, for example z=3.0:xy")
            left, dirs = item.split(":", 1)
            axis, value = left.split("=", 1)
            axis = axis.strip().lower()
            if axis not in axis_to_col:
                die(f"bad fixed coordinate axis {axis!r}; use x, y, or z")
            cutoff = float(value)
            mask = direction_mask(dirs)
            col = axis_to_col[axis]
            for idx, pos in enumerate(positions):
                matched = pos[col] < cutoff if relation == "<" else pos[col] > cutoff
                if matched:
                    flags[idx] = apply_mask(flags[idx], mask)
    return flags


def read_atoms(path: Path, supercell: tuple[int, int, int] = (1, 1, 1)):
    atoms = read(path)
    if any(n != 1 for n in supercell):
        atoms = atoms.repeat(supercell)
    return atoms


def make_input(
    suffix: str,
    pseudo_dir: Path,
    orbital_dir: Path,
    calculation: str,
    basis_type: str,
    device: str,
    ks_solver: str,
    kspacing: float,
    ecutwfc: float,
    nspin: int,
    cal_stress: bool,
    extra: dict[str, str],
) -> str:
    params: list[tuple[str, object]] = [
        ("suffix", suffix),
        ("calculation", calculation),
        ("stru_file", "STRU"),
        ("pseudo_dir", pseudo_dir),
    ]
    if basis_type == "lcao":
        params.append(("orbital_dir", orbital_dir))
    params += [
        ("basis_type", basis_type),
        ("ks_solver", ks_solver),
        ("device", device),
        ("symmetry", 0),
        ("gamma_only", 0),
        ("kspacing", kspacing),
        ("dft_functional", "PBE"),
        ("ecutwfc", ecutwfc),
        ("nspin", nspin),
        ("scf_thr", "1e-6"),
        ("scf_nmax", 300),
        ("smearing_method", "gauss"),
        ("smearing_sigma", 0.015),
        ("mixing_type", "broyden"),
        ("mixing_beta", 0.10),
        ("mixing_ndim", 20),
        ("cal_force", 1),
        ("cal_stress", 1 if cal_stress else 0),
        ("out_wfc_lcao", 0),
    ]
    if "out_chg" not in extra:
        params.append(("out_chg", 0))
    if calculation == "md":
        params += [
            ("md_nstep", extra.pop("md_nstep", "1000")),
            ("md_dt", extra.pop("md_dt", "1.0")),
            ("md_type", extra.pop("md_type", "nvt")),
            ("md_tfirst", extra.pop("md_tfirst", "300")),
            ("md_tlast", extra.pop("md_tlast", "300")),
            ("md_dumpfreq", extra.pop("md_dumpfreq", "1")),
            ("dump_force", extra.pop("dump_force", "1")),
            ("dump_virial", extra.pop("dump_virial", "1" if cal_stress else "0")),
        ]
    if nspin == 2:
        append_spin2_mixing_defaults(params, extra)
    params.extend((k, v) for k, v in extra.items())

    width = max(len(k) for k, _ in params) + 2
    body = ["INPUT_PARAMETERS"]
    body.extend(f"{k:<{width}}{v}" for k, v in params)
    return "\n".join(body) + "\n"


def input_template_params(args) -> list[tuple[str, object]]:
    orbital_dir = args.orbital_dir or DEFAULT_ORBITAL_DIRS[args.orbital_quality]
    extra = parse_key_values(args.set)
    params: list[tuple[str, object]] = [
        ("suffix", args.suffix),
        ("calculation", args.kind),
        ("stru_file", "STRU"),
        ("pseudo_dir", args.pseudo_dir),
    ]
    if args.basis_type == "lcao":
        params.append(("orbital_dir", orbital_dir))
    params += [
        ("basis_type", args.basis_type),
        ("ks_solver", args.ks_solver),
        ("device", args.device),
        ("symmetry", 0),
        ("gamma_only", 0),
        ("kspacing", args.kspacing),
        ("dft_functional", "PBE"),
        ("ecutwfc", args.ecutwfc),
        ("nspin", args.nspin),
        ("scf_thr", "1e-6"),
        ("scf_nmax", 300),
        ("smearing_method", "gauss"),
        ("smearing_sigma", 0.015),
        ("mixing_type", "broyden"),
        ("mixing_beta", 0.10),
        ("mixing_ndim", 20),
        ("cal_force", 1),
        ("cal_stress", 1 if args.kind == "relax" or args.cal_stress else 0),
        ("out_wfc_lcao", 0),
    ]
    if "out_chg" not in extra:
        params.append(("out_chg", 0))
    if args.kind == "relax":
        params += [
            ("relax_method", "cg"),
            ("relax_nmax", args.relax_nmax),
            ("force_thr_ev", args.force_thr_ev),
            ("stress_thr", args.stress_thr),
        ]
    if args.dos:
        if "out_dos" not in extra:
            params.append(("out_dos", 2 if args.basis_type == "lcao" else 1))
        if "dos_sigma" not in extra:
            params.append(("dos_sigma", 0.07))
        if "dos_edelta_ev" not in extra:
            params.append(("dos_edelta_ev", 0.01))
    if args.nspin == 2:
        append_spin2_mixing_defaults(params, extra)
    params.extend((k, v) for k, v in extra.items())
    return params


def format_input_params(params: list[tuple[str, object]], with_comments: bool = True) -> str:
    width = max(len(k) for k, _ in params) + 2
    lines = ["INPUT_PARAMETERS"]
    if with_comments:
        lines += [
            "# Template generated from ABACUS documented INPUT keywords.",
            "# Check pseudo_dir/orbital_dir, k-point density, ecutwfc, and spin settings before production.",
        ]
    lines.extend(f"{k:<{width}}{v}" for k, v in params)
    return "\n".join(lines) + "\n"


def cmd_input_template(args) -> None:
    text = format_input_params(input_template_params(args), with_comments=not args.no_comments)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"wrote {args.kind} INPUT template to {args.out}")


def write_kpt(path: Path, mesh: tuple[int, int, int], shift: tuple[int, int, int], model: str) -> None:
    model_text = "Gamma" if model.lower() == "gamma" else "MP"
    lines = [
        "K_POINTS",
        "0",
        model_text,
        f"{mesh[0]} {mesh[1]} {mesh[2]} {shift[0]} {shift[1]} {shift[2]}",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def cmd_kpt(args) -> None:
    write_kpt(args.out, tuple(args.mesh), tuple(args.shift), args.model)
    print(f"wrote {args.model} KPT to {args.out}")


def default_cpu_bind(mpi_np: int) -> str:
    if mpi_np <= 1:
        return "0-11"
    if mpi_np == 2:
        return "1-24"
    if mpi_np == 4:
        return "1-51"
    return "0-51"


def build_abacus_run_cmd(mpi_np: int, use_numactl: bool, cpu_bind: str | None, mem_bind: str | None) -> str:
    command = f"mpirun -np {mpi_np} --bind-to none abacus"
    if not use_numactl:
        return command
    cpu_bind = cpu_bind or default_cpu_bind(mpi_np)
    mem_bind = "0" if mem_bind is None else mem_bind
    return f"numactl --physcpubind={cpu_bind} --membind={mem_bind} {command}"


def write_run_script(
    path: Path,
    abacus_env: Path,
    mpi_np: int,
    gpu_ids: str | None,
    omp_threads: int = 12,
    use_numactl: bool = True,
    cpu_bind: str | None = None,
    mem_bind: str | None = "0",
) -> None:
    gpu_ids = gpu_ids or "0"
    run_cmd = build_abacus_run_cmd(mpi_np, use_numactl, cpu_bind, mem_bind)
    text = f"""#!/usr/bin/env bash
set -euo pipefail

source "{abacus_env}"
export CUDA_DEVICE_ORDER="${{CUDA_DEVICE_ORDER:-PCI_BUS_ID}}"
export CUDA_VISIBLE_DEVICES="${{CUDA_VISIBLE_DEVICES:-{gpu_ids}}}"
export OMP_NUM_THREADS="${{OMP_NUM_THREADS:-{omp_threads}}}"
export OPENBLAS_NUM_THREADS="${{OPENBLAS_NUM_THREADS:-1}}"
export MKL_NUM_THREADS="${{MKL_NUM_THREADS:-1}}"
export BLIS_NUM_THREADS="${{BLIS_NUM_THREADS:-1}}"
export NUMEXPR_NUM_THREADS="${{NUMEXPR_NUM_THREADS:-1}}"
unset OMP_PROC_BIND OMP_PLACES KMP_AFFINITY
ulimit -s unlimited

RUN_CMD="${{RUN_CMD:-{run_cmd}}}"
echo "$RUN_CMD" > run_cmd.txt
{{
  echo "START_TIME=$(date '+%F %T')"
  echo "ABACUS=$(command -v abacus)"
  echo "ABACUS_ENV={abacus_env}"
  echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
  echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
  echo "RUN_CMD=$RUN_CMD"
  /usr/bin/time -p bash -lc "$RUN_CMD"
  echo "END_TIME=$(date '+%F %T')"
}} > abacus.log 2>&1
touch DONE
"""
    path.write_text(text)
    path.chmod(0o755)


def parse_key_values(items: list[str] | None) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in items or []:
        if "=" not in item:
            die(f"bad key=value item: {item}")
        k, v = item.split("=", 1)
        result[k.strip()] = v.strip()
    return result


def params_have_key(params: list[tuple[str, object]], key: str) -> bool:
    return any(name == key for name, _ in params)


def parse_float_value(value: object, fallback: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return fallback


def append_spin2_mixing_defaults(
    params: list[tuple[str, object]],
    extra: dict[str, str],
) -> None:
    mixing_beta = parse_float_value(extra.get("mixing_beta", 0.10), 0.10)
    defaults: list[tuple[str, object]] = [
        ("mixing_beta_mag", min(4.0 * mixing_beta, 1.6)),
        ("mixing_gg0", 1.0),
        ("mixing_gg0_mag", 0.0),
        ("mixing_gg0_min", 0.1),
        ("mixing_restart", 0),
        ("mixing_dmr", "false"),
    ]
    for key, value in defaults:
        if key not in extra and not params_have_key(params, key):
            params.append((key, value))


def structure_options(args) -> dict:
    return {
        "magnetism": parse_magnetism(args.mag),
        "move_flags": None,
        "pseudo_overrides": parse_element_paths(args.pseudo_file, "upf"),
        "orbital_overrides": parse_element_paths(args.orbital_file, "orb"),
        "element_orbital_quality": parse_element_quality(args.element_orbital_quality),
    }


def cmd_cif2stru(args) -> None:
    orbital_dir = args.orbital_dir or DEFAULT_ORBITAL_DIRS[args.orbital_quality]
    atoms = read_atoms(args.cif, args.supercell)
    options = structure_options(args)
    options["move_flags"] = build_move_flags(
        atoms,
        fix_element=args.fix_element,
        fix_index=args.fix_index,
        fix_below=args.fix_below,
        fix_above=args.fix_above,
    )
    meta = write_stru(
        atoms=atoms,
        output=args.output,
        pseudo_dir=args.pseudo_dir,
        orbital_dir=orbital_dir,
        basis_type=args.basis_type,
        **options,
    )
    args.output.with_suffix(".abacus_lib.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"wrote {args.output}")


def iter_cifs(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(sorted(path.glob("*.cif"), key=lambda p: natural_key(p.name)))
        elif path.is_file():
            files.append(path)
        else:
            die(f"input path not found: {path}")
    if not files:
        die("no CIF files found")
    return files


def cmd_prepare_abacus(args) -> None:
    orbital_dir = args.orbital_dir or DEFAULT_ORBITAL_DIRS[args.orbital_quality]
    cifs = iter_cifs(args.cifs)
    args.out.mkdir(parents=True, exist_ok=True)
    extra = parse_key_values(args.set)
    for i, cif in enumerate(cifs):
        job_dir = args.out / f"{i:06d}"
        job_dir.mkdir(parents=True, exist_ok=True)
        atoms = read_atoms(cif, args.supercell)
        options = structure_options(args)
        options["move_flags"] = build_move_flags(
            atoms,
            fix_element=args.fix_element,
            fix_index=args.fix_index,
            fix_below=args.fix_below,
            fix_above=args.fix_above,
        )
        meta = write_stru(
            atoms=atoms,
            output=job_dir / "STRU",
            pseudo_dir=args.pseudo_dir,
            orbital_dir=orbital_dir,
            basis_type=args.basis_type,
            **options,
        )
        input_text = make_input(
            suffix=args.suffix,
            pseudo_dir=args.pseudo_dir,
            orbital_dir=orbital_dir,
            calculation=args.calculation,
            basis_type=args.basis_type,
            device=args.device,
            ks_solver=args.ks_solver,
            kspacing=args.kspacing,
            ecutwfc=args.ecutwfc,
            nspin=args.nspin,
            cal_stress=args.cal_stress,
            extra=dict(extra),
        )
        (job_dir / "INPUT").write_text(input_text)
        write_run_script(
            job_dir / "run_abacus.sh",
            args.abacus_env,
            args.mpi_np,
            args.gpu_ids,
            omp_threads=args.omp_threads,
            use_numactl=not args.no_numactl,
            cpu_bind=args.cpu_bind,
            mem_bind=args.mem_bind,
        )
        shutil.copy2(cif, job_dir / "source.cif")
        (job_dir / "metadata.json").write_text(
            json.dumps({"source": str(cif), "library": meta}, indent=2) + "\n"
        )
    write_job_array(args.out, len(cifs))
    print(f"prepared {len(cifs)} ABACUS jobs under {args.out}")


def write_job_array(out: Path, count: int) -> None:
    text = f"""#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
for job in "$root"/[0-9][0-9][0-9][0-9][0-9][0-9]; do
  [ -d "$job" ] || continue
  echo "running $(basename "$job")"
  (cd "$job" && bash run_abacus.sh)
done
"""
    path = out / "run_all_abacus.sh"
    path.write_text(text)
    path.chmod(0o755)
    (out / "jobs.txt").write_text("\n".join(f"{i:06d}" for i in range(count)) + "\n")


def abacus_env_to_root(env_path: Path) -> Path:
    return env_path.parent.parent


def discover_abacus_versions(root: Path = DEFAULT_ABACUS_ROOT) -> list[dict[str, object]]:
    versions = []
    for path in sorted(root.glob("abacus-*"), key=lambda p: natural_key(p.name)):
        if not path.is_dir():
            continue
        env = path / "toolchain" / "abacus_env.sh"
        binary = path / "bin" / "abacus"
        basic = path / "bin" / "abacus_basic_gpu"
        if not env.is_file():
            continue
        versions.append(
            {
                "label": path.name,
                "root": str(path),
                "env": str(env),
                "binary": str(binary if binary.exists() else basic if basic.exists() else ""),
                "is_default": env.resolve() == DEFAULT_ABACUS_ENV.resolve(),
                "mtime": int(path.stat().st_mtime),
            }
        )
    return versions


def cmd_abacus_versions(args) -> None:
    versions = discover_abacus_versions(args.root)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(versions, indent=2) + "\n")
    if not versions:
        die(f"no ABACUS versions found under {args.root}")
    print(f"{'default':<7}  {'label':<58}  env")
    print(f"{'-' * 7}  {'-' * 58}  {'-' * 20}")
    for item in versions:
        marker = "*" if item["is_default"] else ""
        print(f"{marker:<7}  {item['label']:<58}  {item['env']}")


def write_absolute_job_array(path: Path, jobs: list[Path]) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "JOBS=(",
    ]
    lines += [f'  "{job.resolve()}"' for job in jobs]
    lines += [
        ")",
        "",
        'for job in "${JOBS[@]}"; do',
        '  [ -d "$job" ] || { echo "missing job directory: $job"; continue; }',
        '  echo "running $job"',
        '  (cd "$job" && bash run_abacus.sh)',
        "done",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    path.chmod(0o755)


def cmd_launch_script(args) -> None:
    jobs = iter_job_dirs(args.jobs)
    if not jobs:
        die("no ABACUS job directories found")
    for job in jobs:
        write_run_script(
            job / args.script_name,
            args.abacus_env,
            args.mpi_np,
            args.gpu_ids,
            omp_threads=args.omp_threads,
            use_numactl=not args.no_numactl,
            cpu_bind=args.cpu_bind,
            mem_bind=args.mem_bind,
        )
    if args.array_script:
        write_absolute_job_array(args.array_script, jobs)
    print(
        f"wrote {args.script_name} for {len(jobs)} jobs using "
        f"{abacus_env_to_root(args.abacus_env).name}"
    )
    if args.array_script:
        print(f"wrote launcher {args.array_script}")


def random_strained(atoms, strain: float, rng: random.Random):
    new = atoms.copy()
    eps = np.array([[rng.uniform(-strain, strain) for _ in range(3)] for _ in range(3)])
    eps = 0.5 * (eps + eps.T)
    cell = np.dot(np.eye(3) + eps, atoms.cell.array)
    new.set_cell(cell, scale_atoms=True)
    return new


def cmd_make_candidates(args) -> None:
    rng = random.Random(args.seed)
    atoms0 = read_atoms(args.cif, args.supercell)
    args.out.mkdir(parents=True, exist_ok=True)
    for i in range(args.count):
        atoms = random_strained(atoms0, args.strain, rng) if args.strain > 0 else atoms0.copy()
        if args.rattle > 0:
            disp = np.array(
                [[rng.gauss(0.0, args.rattle) for _ in range(3)] for _ in range(len(atoms))]
            )
            atoms.set_positions(atoms.get_positions() + disp)
        write(args.out / f"candidate_{i:06d}.cif", atoms)
    print(f"wrote {args.count} candidate CIFs to {args.out}")


def find_abacus_input(job: Path) -> Path | None:
    job = job.resolve() if job.exists() else job
    candidates = []
    if job.is_file() and job.name == "INPUT":
        candidates.append(job)
    elif job.is_dir():
        candidates.append(job / "INPUT")
        if job.name.startswith("OUT."):
            candidates.append(job.parent / "INPUT")
    for input_file in candidates:
        if input_file.is_file():
            return input_file
    return None


def read_input_param(job: Path, key: str) -> str | None:
    input_file = find_abacus_input(job)
    if not input_file:
        return None
    key = key.lower()
    for raw in input_file.read_text(errors="ignore").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        words = line.split()
        if len(words) >= 2 and words[0].lower() == key:
            return words[1]
    return None


def read_input_calculation(job: Path) -> str:
    return (read_input_param(job, "calculation") or "unknown").lower()


def find_abacus_outdir(job: Path) -> Path | None:
    job = job.resolve() if job.exists() else job
    if job.name.startswith("OUT.") and job.is_dir():
        return job
    input_file = find_abacus_input(job)
    if input_file:
        suffix = read_input_param(job, "suffix") or "ABACUS"
        candidate = job / f"OUT.{suffix}"
        if candidate.is_dir():
            return candidate
    outdirs = sorted([p for p in job.glob("OUT.*") if p.is_dir()], key=lambda p: natural_key(p.name))
    return outdirs[0] if outdirs else None


def find_running_log(job: Path) -> Path | None:
    job = job.resolve() if job.exists() else job
    outdir = find_abacus_outdir(job)
    if outdir:
        logs = sorted(outdir.glob("running_*.log"))
        if logs:
            return max(logs, key=lambda p: (p.stat().st_mtime, natural_key(p.name)))
    direct = sorted(job.glob("running_*.log"))
    return max(direct, key=lambda p: (p.stat().st_mtime, natural_key(p.name))) if direct else None


def last_pattern_position(text: str, patterns: list[str]) -> int:
    last = -1
    for pattern in patterns:
        matches = list(re.finditer(pattern, text, flags=re.IGNORECASE))
        if matches:
            last = max(last, matches[-1].start())
    return last


def parse_last_iter_energy(text: str) -> float | None:
    etot_index = None
    last_energy = None
    for line in text.splitlines():
        words = line.split()
        if not words:
            continue
        if "ETOT/eV" in words:
            etot_index = words.index("ETOT/eV")
            continue
        if etot_index is None or len(words) <= etot_index:
            continue
        try:
            last_energy = float(words[etot_index])
        except ValueError:
            continue
    return last_energy


def parse_final_energy(text: str) -> float | None:
    energy_patterns = [
        r"!FINAL_ETOT_IS\s+([-+0-9.eE]+)",
        r"final\s+etot\s+is\s+([-+0-9.eE]+)",
        r"final\s+energy\s+is\s+([-+0-9.eE]+)",
        r"TOTAL\s+ENERGY\s*=\s*([-+0-9.eE]+)",
        r"final_etot\s*[:=]\s*([-+0-9.eE]+)",
        r"\bETOT\s*[:=]\s*([-+0-9.eE]+)\s*eV",
    ]
    for pattern in energy_patterns:
        matches = re.findall(pattern, text, flags=re.IGNORECASE)
        if matches:
            return float(matches[-1])
    return None


def parse_abacus_status(job: Path) -> dict:
    job = job.resolve() if job.exists() else job
    outdir = find_abacus_outdir(job)
    log = find_running_log(job)
    row = {
        "job": str(job),
        "calculation": read_input_calculation(job),
        "outdir": str(outdir or ""),
        "log": str(log or ""),
        "finished": False,
        "converged": False,
        "failed": False,
        "energy_ev": None,
        "message": "",
    }
    if not log or not log.is_file():
        row["message"] = "missing running_*.log"
        return row
    text = log.read_text(errors="ignore")
    lower = text.lower()
    row["finished"] = "finish time" in lower or "total  time" in lower
    negative = [
        r"scf\s+is\s+not\s+converged",
        r"\bnot\s+converged\b",
        r"convergence\s+has\s+not\s+been\s+achieved",
        r"convergence\s+has\s+not\s+achieved",
        r"convergence\s+has\s+not\s+been\s+reached",
    ]
    hard_failed = any(pat in lower for pat in ["warning_quit", "segmentation fault"]) or bool(
        re.search(r"(^|\n)\s*(error|fatal)\b", text, flags=re.IGNORECASE)
    )
    positive = [
        r"charge\s+density\s+convergence\s+is\s+achieved",
        r"\bconvergence\s+is\s+achieved",
        r"relaxation\s+is\s+converged",
        r"cell\s+relaxation\s+is\s+converged",
        r"cell-relax\s+is\s+converged",
        r"lattice\s+relaxation\s+is\s+converged",
        r"geometry\s+optimization\s+is\s+converged",
        r"force\s+convergence\s+is\s+achieved",
        r"stress\s+convergence\s+is\s+achieved",
    ]
    last_negative = last_pattern_position(text, negative)
    last_positive = last_pattern_position(text, positive)
    row["converged"] = last_positive >= 0 and not hard_failed and last_positive > last_negative
    row["failed"] = hard_failed or (last_negative >= 0 and last_negative > last_positive)

    row["energy_ev"] = parse_final_energy(text)
    if row["energy_ev"] is None and row["converged"]:
        row["energy_ev"] = parse_last_iter_energy(text)
    if row["failed"]:
        row["message"] = "failed or not converged"
    elif row["finished"] and row["converged"]:
        row["message"] = "finished and converged"
    elif row["finished"]:
        row["message"] = "finished but convergence was not detected"
    else:
        row["message"] = "not finished"
    return row


def iter_job_dirs(paths: list[Path]) -> list[Path]:
    jobs: list[Path] = []
    for path in paths:
        path = path.expanduser()
        if find_abacus_input(path) or find_abacus_outdir(path) or find_running_log(path):
            jobs.append(path)
        elif path.is_dir():
            children = [
                p
                for p in path.iterdir()
                if p.is_dir() and (find_abacus_input(p) or find_abacus_outdir(p) or find_running_log(p))
            ]
            jobs.extend(sorted(children, key=lambda p: natural_key(p.name)))
        else:
            die(f"path not found: {path}")
    return jobs


def cmd_check_abacus(args) -> None:
    jobs = iter_job_dirs(args.jobs)
    if not jobs:
        die("no ABACUS jobs found")
    rows = [parse_abacus_status(job) for job in jobs]
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2) + "\n")
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    width = max(len(Path(r["job"]).name or str(r["job"])) for r in rows)
    out_width = max(6, max(len(Path(r["outdir"]).name) if r["outdir"] else 0 for r in rows))
    print(
        f"{'job':<{width}}  {'type':<8}  {'outdir':<{out_width}}  "
        "finished  converged  failed  energy_ev       message"
    )
    for row in rows:
        energy = "" if row["energy_ev"] is None else f"{row['energy_ev']:.8f}"
        job_name = Path(row["job"]).name or str(row["job"])
        outdir = Path(row["outdir"]).name if row["outdir"] else ""
        print(
            f"{job_name:<{width}}  {row['calculation']:<8}  {outdir:<{out_width}}  "
            f"{str(row['finished']):<8}  "
            f"{str(row['converged']):<9}  {str(row['failed']):<6}  {energy:<14}  {row['message']}"
        )


def read_abacus_input(path: Path) -> dict[str, str]:
    input_file = path if path.is_file() else path / "INPUT"
    if not input_file.is_file():
        outdir = find_abacus_outdir(path)
        if outdir and (outdir / "INPUT").is_file():
            input_file = outdir / "INPUT"
        else:
            return {}
    params: dict[str, str] = {}
    started = False
    for raw in input_file.read_text(errors="ignore").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line == "INPUT_PARAMETERS":
            started = True
            continue
        if not started:
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            params[parts[0].lower()] = parts[1].strip()
    return params


def read_kpt_mesh(path: Path) -> str | None:
    kpt = path if path.is_file() else path / "KPT"
    if not kpt.is_file():
        return None
    lines = [line.strip() for line in kpt.read_text(errors="ignore").splitlines() if line.strip()]
    if len(lines) >= 4 and lines[1].split()[0] == "0":
        return " ".join(lines[3].split()[:3])
    return None


def parse_last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.IGNORECASE)
    if not matches:
        return None
    value = matches[-1]
    if isinstance(value, tuple):
        value = next((x for x in value if x), value[-1])
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_time_seconds(text: str, label: str) -> float | None:
    match = re.search(rf"{re.escape(label)}\s*:\s*([0-9.]+)", text, flags=re.IGNORECASE)
    if match:
        return float(match.group(1))
    return None


def parse_abacus_metrics(job: Path) -> dict:
    row = parse_abacus_status(job)
    text = ""
    log = Path(row["log"]) if row["log"] else None
    if log and log.is_file():
        text = log.read_text(errors="ignore")
    input_params = read_abacus_input(job)
    outdir = Path(row["outdir"]) if row["outdir"] else None
    metrics = dict(row)
    metrics["normal_end"] = bool(row["finished"])
    metrics["converge"] = bool(row["converged"])
    metrics["INPUT"] = input_params
    for key in ("calculation", "basis_type", "ks_solver", "device", "ecutwfc", "kspacing", "nspin"):
        metrics[key] = input_params.get(key)
    metrics["kpt"] = read_kpt_mesh(job)
    if outdir:
        metrics["kpt"] = metrics["kpt"] or read_kpt_mesh(outdir)

    version = "unknown"
    commit = "unknown"
    for idx, line in enumerate(text.splitlines()):
        if "WELCOME TO ABACUS" in line:
            version = line.split()[-1]
        if "Commit:" in line:
            commit = line.split(":", 1)[1].strip()
        if version != "unknown" and commit != "unknown":
            break
    metrics["version"] = version if commit == "unknown" else f"{version}({commit})"
    metrics["ncore"] = int(parse_last_float(r"DSIZE\s*=\s*([-+0-9.eE]+)", text) or 0) or None
    metrics["nbands"] = int(parse_last_float(r"NBANDS\s*=\s*([-+0-9.eE]+)", text) or 0) or None
    metrics["nkstot"] = int(parse_last_float(r"nkstot\s*=\s*([-+0-9.eE]+)", text) or 0) or None
    metrics["ibzk"] = int(parse_last_float(r"nkstot_ibz\s*=\s*([-+0-9.eE]+)", text) or 0) or None
    metrics["natom"] = None
    natom_total = 0
    for match in re.findall(r"number of atom for this type\s*=\s*([0-9]+)", text, flags=re.IGNORECASE):
        natom_total += int(match)
    if natom_total:
        metrics["natom"] = natom_total

    metrics["energy"] = row["energy_ev"]
    metrics["energy_ks"] = parse_last_float(r"E_KohnSham\s+([-+0-9.eE]+)", text)
    metrics["volume"] = parse_last_float(r"Volume\s*\(A\^3\)\s*=\s*([-+0-9.eE]+)", text)
    metrics["efermi"] = parse_last_float(r"E_Fermi(?:_up|_dw)?\s+([-+0-9.eE]+)", text)
    if metrics["natom"] and metrics["energy"] is not None:
        metrics["energy_per_atom"] = metrics["energy"] / metrics["natom"]
    else:
        metrics["energy_per_atom"] = None

    metrics["scf_steps"] = len(re.findall(r"^\s*ITER\s+", text, flags=re.MULTILINE))
    if metrics["scf_steps"] == 0:
        metrics["scf_steps"] = len(re.findall(r"charge density convergence", text, flags=re.IGNORECASE)) or None
    metrics["total_time"] = parse_time_seconds(text, "Total  Time")
    metrics["scf_time"] = parse_time_seconds(text, "SCF  Time")
    if metrics["scf_time"] is None:
        metrics["scf_time"] = parse_time_seconds(text, "SCF Time")

    metrics["total_mag"] = parse_last_float(r"total magnetism \(Bohr mag/cell\)\s*=\s*([-+0-9.eE]+)", text)
    metrics["absolute_mag"] = parse_last_float(r"absolute magnetism\s*=\s*([-+0-9.eE]+)", text)
    largest_force = re.findall(
        r"(?:Largest gradient in force|Largest gradient is)\s*(?:is)?\s*([-+0-9.eE]+)",
        text,
        flags=re.IGNORECASE,
    )
    metrics["largest_gradient"] = float(largest_force[-1]) if largest_force else None
    largest_stress = re.findall(r"Largest gradient in stress is\s*([-+0-9.eE]+)", text, flags=re.IGNORECASE)
    metrics["largest_gradient_stress"] = float(largest_stress[-1]) if largest_stress else None
    stress_blocks = re.findall(
        r"TOTAL-STRESS \(KBAR\).*?\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)",
        text,
        flags=re.IGNORECASE | re.S,
    )
    if stress_blocks:
        stress = [float(x) for x in stress_blocks[-1]]
        metrics["pressure"] = (stress[0] + stress[4] + stress[8]) / 3.0
    else:
        metrics["pressure"] = None
    return metrics


def cmd_collect_metrics(args) -> None:
    jobs = iter_job_dirs(args.jobs)
    if not jobs:
        die("no ABACUS jobs found")
    rows = [parse_abacus_metrics(job) for job in jobs]
    keys = args.metrics or [
        "job",
        "normal_end",
        "converge",
        "energy",
        "energy_per_atom",
        "natom",
        "ecutwfc",
        "kspacing",
        "kpt",
        "nspin",
        "total_mag",
        "efermi",
        "total_time",
        "scf_steps",
        "message",
    ]
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2) + "\n")
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    print_metrics_table(rows, keys)


def format_metric_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.8g}"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def print_metrics_table(rows: list[dict], keys: list[str]) -> None:
    widths = {}
    for key in keys:
        widths[key] = max(len(key), *(len(format_metric_value(row.get(key))) for row in rows))
        widths[key] = min(widths[key], 28)
    header = "  ".join(f"{key:<{widths[key]}}" for key in keys)
    print(header)
    print("  ".join("-" * widths[key] for key in keys))
    for row in rows:
        values = []
        for key in keys:
            value = format_metric_value(row.get(key))
            if len(value) > widths[key]:
                value = value[: widths[key] - 1] + "~"
            values.append(f"{value:<{widths[key]}}")
        print("  ".join(values))


def cmd_report_metrics(args) -> None:
    rows = json.loads(args.metrics.read_text())
    if isinstance(rows, dict):
        rows = list(rows.values())
    if not isinstance(rows, list):
        die("metrics JSON must contain a list of result rows")
    keys = args.keys or [
        "job",
        "normal_end",
        "converge",
        "energy",
        "energy_per_atom",
        "natom",
        "ecutwfc",
        "kspacing",
        "kpt",
        "nspin",
        "total_mag",
        "efermi",
        "total_time",
        "scf_steps",
        "message",
    ]
    total = len(rows)
    converged = sum(1 for row in rows if row.get("converge"))
    failed = sum(1 for row in rows if row.get("failed"))
    head = "".join(f"<th>{html.escape(key)}</th>" for key in keys)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(format_metric_value(row.get(key)))}</td>" for key in keys)
        cls = "ok" if row.get("converge") else "bad" if row.get("failed") else ""
        body.append(f"<tr class='{cls}'>{cells}</tr>")
    doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>abacuskit metrics report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #1f2937; }}
    h1 {{ margin: 0 0 8px; }}
    .summary {{ display: flex; gap: 12px; margin: 16px 0 24px; }}
    .pill {{ border: 1px solid #d1d5db; border-radius: 8px; padding: 10px 14px; background: #f9fafb; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
    th, td {{ border: 1px solid #d1d5db; padding: 6px 8px; text-align: left; vertical-align: top; }}
    th {{ background: #eef2ff; position: sticky; top: 0; }}
    tr.ok td:first-child {{ border-left: 4px solid #16a34a; }}
    tr.bad td:first-child {{ border-left: 4px solid #dc2626; }}
  </style>
</head>
<body>
  <h1>abacuskit metrics report</h1>
  <div class="summary">
    <div class="pill">Total jobs: <strong>{total}</strong></div>
    <div class="pill">Converged: <strong>{converged}</strong></div>
    <div class="pill">Failed: <strong>{failed}</strong></div>
  </div>
  <table>
    <thead><tr>{head}</tr></thead>
    <tbody>
      {''.join(body)}
    </tbody>
  </table>
</body>
</html>
"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(doc)
    print(f"wrote metrics report to {args.out}")


def update_input_key(path: Path, key: str, value: str) -> None:
    input_file = path / "INPUT"
    if not input_file.is_file():
        die(f"INPUT not found in {path}")
    lines = input_file.read_text(errors="ignore").splitlines()
    output = []
    changed = False
    for line in lines:
        stripped = line.split("#", 1)[0].strip()
        parts = stripped.split(None, 1)
        if len(parts) >= 1 and parts[0].lower() == key.lower():
            output.append(f"{key:<18}{value}")
            changed = True
        else:
            output.append(line)
    if not changed:
        output.append(f"{key:<18}{value}")
    input_file.write_text("\n".join(output) + "\n")


def sanitize_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.+-]+", "_", value.strip())
    return label.strip("_") or "value"


def copy_job_template(src: Path, dst: Path, force: bool) -> None:
    if dst.exists():
        if not force:
            die(f"target exists: {dst}; pass --force to overwrite")
        shutil.rmtree(dst)
    ignore = shutil.ignore_patterns("OUT.*", "DONE", "abacus.log", "run_cmd.txt", "*.tmp")
    shutil.copytree(src, dst, ignore=ignore)
    run_script = dst / "run_abacus.sh"
    if not run_script.is_file():
        write_run_script(run_script, DEFAULT_ABACUS_ENV, 1, None)


def write_directory_array(out: Path, dirs: list[Path]) -> None:
    text = """#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
while IFS= read -r job; do
  [ -n "$job" ] || continue
  echo "running $job"
  (cd "$root/$job" && bash run_abacus.sh)
done < "$root/jobs.txt"
"""
    path = out / "run_all_abacus.sh"
    path.write_text(text)
    path.chmod(0o755)
    (out / "jobs.txt").write_text("\n".join(p.name for p in dirs) + "\n")


def cmd_conv_test(args) -> None:
    jobs = iter_job_dirs(args.jobs)
    if not jobs:
        die("no ABACUS template jobs found")
    args.out.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    for job in jobs:
        base = job.name
        for value in args.values:
            target = args.out / f"{base}_{args.key}_{sanitize_label(value)}"
            copy_job_template(job, target, args.force)
            if args.key.lower() == "kpt":
                parts = [int(x) for x in value.replace(",", " ").split()]
                if len(parts) not in {3, 6}:
                    die("kpt values must be like '3 3 1' or '3 3 1 0 0 0'")
                mesh = tuple(parts[:3])
                shift = tuple(parts[3:] if len(parts) == 6 else [0, 0, 0])
                write_kpt(target / "KPT", mesh, shift, args.kpt_model)
                update_input_key(target, "gamma_only", 0)
            else:
                update_input_key(target, args.key, value)
            created.append(target)
    write_directory_array(args.out, created)
    manifest = {
        "key": args.key,
        "values": args.values,
        "source_jobs": [str(job) for job in jobs],
        "created_jobs": [str(job) for job in created],
    }
    (args.out / "conv_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"prepared {len(created)} convergence-test jobs under {args.out}")


def collect_one_with_dpdata(job: Path, fmt: str, out: Path, set_size: int) -> int:
    import dpdata

    system = dpdata.LabeledSystem(str(job), fmt=fmt)
    if len(system) == 0:
        return 0
    system.to("deepmd/npy", str(out), set_size=set_size)
    return len(system)


def cmd_collect_deepmd(args) -> None:
    python = args.python
    already_respawned = os.environ.get("ABACUS_DEEPMD_FLOW_RESPAWNED") == "1"
    if python and Path(python).resolve() != Path(sys.executable).resolve() and not already_respawned:
        cmd = [str(python), str(Path(__file__).resolve()), "collect-deepmd"]
        passthrough = [str(x) for x in args.jobs]
        cmd += passthrough + ["--out", str(args.out), "--set-size", str(args.set_size)]
        if args.fmt:
            cmd += ["--fmt", args.fmt]
        if args.split_ratio != 0.0:
            cmd += ["--split-ratio", str(args.split_ratio)]
        env = dict(os.environ)
        env["ABACUS_DEEPMD_FLOW_RESPAWNED"] = "1"
        raise SystemExit(subprocess.call(cmd, env=env))

    jobs: list[Path] = []
    for path in args.jobs:
        if (path / "INPUT").is_file():
            jobs.append(path)
        else:
            jobs.extend(sorted([p for p in path.iterdir() if (p / "INPUT").is_file()]))
    if not jobs:
        die("no ABACUS job directories found")

    train_dir = args.out / "train"
    valid_dir = args.out / "valid"
    train_dir.mkdir(parents=True, exist_ok=True)
    if args.split_ratio > 0:
        valid_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)
    rng.shuffle(jobs)
    n_valid = math.floor(len(jobs) * args.split_ratio)
    valid_jobs = set(jobs[:n_valid])
    rows = []
    for job in jobs:
        calc = read_input_calculation(job)
        fmt = args.fmt or f"abacus/{'md' if calc == 'md' else 'relax' if calc == 'relax' else 'scf'}"
        target_root = valid_dir if job in valid_jobs else train_dir
        target = target_root / job.name
        try:
            frames = collect_one_with_dpdata(job, fmt, target, args.set_size)
            rows.append({"job": str(job), "target": str(target), "fmt": fmt, "frames": frames})
        except Exception as exc:
            rows.append({"job": str(job), "target": str(target), "fmt": fmt, "error": str(exc)})

    (args.out / "collect_report.json").write_text(json.dumps(rows, indent=2) + "\n")
    ok = sum(1 for r in rows if r.get("frames", 0) > 0)
    print(f"converted {ok}/{len(rows)} jobs to DeepMD npy under {args.out}")


def resolve_out_path(path: Path) -> Path:
    if path.is_file():
        return path
    outdir = find_abacus_outdir(path)
    return outdir or path


def find_first_file(root: Path, names: list[str], patterns: list[str]) -> Path | None:
    for name in names:
        candidate = root / name
        if candidate.is_file():
            return candidate
    for pattern in patterns:
        found = sorted([p for p in root.glob(pattern) if p.is_file()], key=lambda p: natural_key(p.name))
        if found:
            return found[0]
    return None


def read_numeric_table(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("<"):
            continue
        parts = line.split()
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    if not rows:
        die(f"no numeric data found in {path}")
    max_cols = max(len(r) for r in rows)
    rows = [r for r in rows if len(r) == max_cols]
    return np.array(rows, dtype=float)


def parse_pdos_file(path: Path) -> tuple[np.ndarray, dict[tuple[str, str], np.ndarray]]:
    text = path.read_text(errors="ignore")
    energy_match = re.search(r"<energy_values[^>]*>(.*?)</energy_values>", text, re.S)
    if not energy_match:
        die(f"cannot find energy_values in PDOS file {path}")
    energies = np.array([float(x) for x in energy_match.group(1).split()], dtype=float)
    groups: dict[tuple[str, str], np.ndarray] = {}
    orbital_re = re.compile(r"<orbital\s+(.*?)>\s*<data>(.*?)</data>", re.S)
    for header, data_text in orbital_re.findall(text):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', header))
        species = attrs.get("species", "").strip()
        try:
            l_value = int(attrs.get("l", "0"))
        except ValueError:
            l_value = 0
        label = L_TO_ORBITAL_LABEL.get(l_value, f"l{l_value}")
        values = []
        for line in data_text.splitlines():
            nums = [float(x) for x in line.split()] if line.split() else []
            if nums:
                values.append(sum(nums))
        if len(values) != len(energies):
            continue
        key = (species, label)
        groups[key] = groups.get(key, np.zeros_like(energies)) + np.array(values, dtype=float)
    if not groups:
        die(f"no PDOS orbital data parsed from {path}")
    return energies, groups


def parse_selectors(items: list[str] | None) -> set[tuple[str, str]]:
    selectors: set[tuple[str, str]] = set()
    for item in items or []:
        chunks = re.split(r"[,;]", item)
        for chunk in chunks:
            chunk = chunk.strip()
            if not chunk:
                continue
            if "=" in chunk:
                sym, orbitals = chunk.split("=", 1)
            elif ":" in chunk:
                sym, orbitals = chunk.split(":", 1)
            else:
                die(f"bad selector {chunk!r}; expected Element=s,p,d, for example Ni=d")
            for orb in orbitals.replace("+", "").replace("/", "").replace(" ", ""):
                orb = orb.lower()
                if orb not in ORBITAL_LABEL_TO_L:
                    die(f"bad orbital label {orb!r}; use s, p, d, f, or g")
                selectors.add((sym.strip(), orb))
    return selectors


def maybe_shift_fermi(x: np.ndarray, fermi: float | None) -> tuple[np.ndarray, str]:
    if fermi is None:
        return x, "Energy (eV)"
    return x - fermi, "Energy - E_F (eV)"


def cmd_plot_dos(args) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = resolve_out_path(args.path)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    kind = args.kind
    if kind == "auto":
        if find_first_file(root, ["PDOS"], ["PDOS*"]) and args.select:
            kind = "pdos"
        elif find_first_file(root, ["LDOS.txt"], ["LDOS*.txt"]):
            kind = "ldos"
        else:
            kind = "dos"

    if kind == "dos":
        dos_file = args.file or find_first_file(root, ["DOS1_smearing.dat", "DOS1"], ["DOS*_smearing.dat", "DOS*"])
        if not dos_file:
            die(f"cannot find DOS file under {root}")
        data = read_numeric_table(dos_file)
        x, xlabel = maybe_shift_fermi(data[:, 0], args.fermi)
        fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=180)
        ax.plot(x, data[:, 1], lw=1.5, label=dos_file.name)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("DOS")
        ax.legend(frameon=False)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(args.out)
        plt.close(fig)
        print(f"wrote DOS plot {args.out}")
        return

    if kind == "pdos":
        pdos_file = args.file or find_first_file(root, ["PDOS"], ["PDOS*"])
        if not pdos_file:
            die(f"cannot find PDOS file under {root}")
        energies, groups = parse_pdos_file(pdos_file)
        selectors = parse_selectors(args.select)
        if selectors:
            groups = {key: val for key, val in groups.items() if key in selectors}
        if not groups:
            die("selected PDOS channels are not present in the PDOS file")
        x, xlabel = maybe_shift_fermi(energies, args.fermi)
        fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)
        for (species, orbital), values in sorted(groups.items()):
            ax.plot(x, values, lw=1.2, label=f"{species}-{orbital}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("PDOS")
        ax.legend(frameon=False, ncol=2)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(args.out)
        plt.close(fig)
        print(f"wrote PDOS plot {args.out}")
        return

    if kind == "ldos":
        ldos_file = args.file or find_first_file(root, ["LDOS.txt"], ["LDOS*.txt", "LDOS_*eV.cube"])
        if not ldos_file:
            die(f"cannot find LDOS.txt or LDOS cube file under {root}")
        fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)
        if ldos_file.suffix.lower() == ".cube":
            values = read_cube_values(ldos_file)
            plane = values[:, :, values.shape[2] // 2]
            im = ax.imshow(plane.T, origin="lower", aspect="auto", cmap="viridis")
            fig.colorbar(im, ax=ax, label="LDOS")
            ax.set_xlabel("grid x")
            ax.set_ylabel("grid y")
            ax.set_title(ldos_file.name)
        else:
            data = read_numeric_table(ldos_file)
            im = ax.imshow(data, origin="lower", aspect="auto", cmap="viridis")
            fig.colorbar(im, ax=ax, label="LDOS")
            ax.set_xlabel("energy grid")
            ax.set_ylabel("line point")
        fig.tight_layout()
        fig.savefig(args.out)
        plt.close(fig)
        print(f"wrote LDOS plot {args.out}")
        return

    die(f"unknown plot kind: {kind}")


def read_cube_values(path: Path) -> np.ndarray:
    _, values = read_cube_grid(path)
    return values


def read_cube_grid(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 6:
        die(f"cube file too short: {path}")
    natoms = abs(int(float(lines[2].split()[0])))
    nx = int(float(lines[3].split()[0]))
    ny = int(float(lines[4].split()[0]))
    nz = int(float(lines[5].split()[0]))
    start = 6 + natoms
    values: list[float] = []
    for line in lines[start:]:
        for part in line.split():
            try:
                values.append(float(part))
            except ValueError:
                pass
    need = nx * ny * nz
    if len(values) < need:
        die(f"cube file has {len(values)} values, expected {need}: {path}")
    header = lines[:start]
    return header, np.array(values[:need], dtype=float).reshape((nx, ny, nz))


def write_cube_grid(path: Path, header: list[str], values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flat = values.reshape(-1)
    lines = list(header)
    for start in range(0, len(flat), 6):
        chunk = flat[start : start + 6]
        lines.append(" ".join(f"{value: .10e}" for value in chunk))
    path.write_text("\n".join(lines) + "\n")


def cube_midplane(values: np.ndarray, axis: str, index: int | None) -> tuple[np.ndarray, int]:
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    size = values.shape[axis_index]
    if index is None:
        index = size // 2
    if not 0 <= index < size:
        die(f"{axis}-axis slice index {index} is outside 0..{size - 1}")
    if axis == "x":
        return values[index, :, :], index
    if axis == "y":
        return values[:, index, :], index
    return values[:, :, index], index


def find_grid_file(root: Path, kind: str) -> Path | None:
    if kind == "elf":
        return find_first_file(root, ["elf.cube"], ["elf*.cube", "*ELF*.cube", "*elf*.cube"])
    if kind == "charge":
        return find_first_file(root, ["SPIN1_CHG.cube", "CHG.cube"], ["*CHG*.cube", "*chg*.cube"])
    return find_first_file(root, [], ["*.cube"])


def plot_grid_slice(values: np.ndarray, out: Path, label: str, axis: str, index: int | None, cmap: str) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plane, used_index = cube_midplane(values, axis, index)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=180)
    im = ax.imshow(plane.T, origin="lower", aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=label)
    ax.set_xlabel("grid")
    ax.set_ylabel("grid")
    ax.set_title(f"{label}, {axis} slice {used_index}")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def cmd_plot_grid(args) -> None:
    root = resolve_out_path(args.path)
    kind = args.kind
    if kind == "auto":
        if args.minus:
            kind = "diff"
        elif find_grid_file(root, "elf"):
            kind = "elf"
        elif find_grid_file(root, "charge"):
            kind = "charge"
        else:
            kind = "cube"

    if kind == "diff":
        plus_file = args.file or (root if root.is_file() else None) or find_grid_file(root, "charge")
        if not plus_file:
            die(f"cannot find charge-density cube under {root}")
        minus_root = resolve_out_path(args.minus) if args.minus else None
        minus_file = args.minus_file or (
            minus_root if minus_root and minus_root.is_file() else find_grid_file(minus_root, "charge") if minus_root else None
        )
        if not minus_file:
            die("charge-density difference needs --minus or --minus-file")
        header, plus = read_cube_grid(plus_file)
        _, minus = read_cube_grid(minus_file)
        if plus.shape != minus.shape:
            die(f"cube grid shapes differ: {plus_file} {plus.shape} vs {minus_file} {minus.shape}")
        values = plus - minus
        label = "Charge density difference"
        cmap = args.cmap or "RdBu_r"
        if args.cube_out:
            write_cube_grid(args.cube_out, header, values)
            print(f"wrote charge-density difference cube {args.cube_out}")
    else:
        if args.file:
            cube_file = args.file
        elif root.is_file():
            cube_file = root
        else:
            cube_file = find_grid_file(root, "elf" if kind == "elf" else "charge" if kind == "charge" else "cube")
        if not cube_file:
            die(f"cannot find {kind} cube under {root}")
        _, values = read_cube_grid(cube_file)
        label = "ELF" if kind == "elf" else "Charge density" if kind == "charge" else cube_file.name
        cmap = args.cmap or ("viridis" if kind == "elf" else "magma")

    plot_grid_slice(values, args.out, label, args.axis, args.index, cmap)
    print(f"wrote {label} plot {args.out}")


def infer_type_map_from_data(paths: list[Path]) -> list[str]:
    names: list[str] = []
    for root in paths:
        type_map = root / "type_map.raw"
        if type_map.is_file():
            for sym in type_map.read_text().split():
                if sym not in names:
                    names.append(sym)
    return names


def deepmd_input(type_map: list[str], train_systems: list[str], valid_systems: list[str], steps: int):
    return {
        "model": {
            "type": "standard",
            "type_map": type_map,
            "descriptor": {
                "type": "se_e2_a",
                "sel": "auto",
                "rcut_smth": 0.5,
                "rcut": 6.0,
                "neuron": [25, 50, 100],
                "resnet_dt": False,
                "axis_neuron": 16,
                "seed": 1,
            },
            "fitting_net": {"type": "ener", "neuron": [240, 240, 240], "resnet_dt": True, "seed": 1},
        },
        "learning_rate": {"type": "exp", "start_lr": 1.0e-3, "stop_lr": 3.51e-8, "decay_steps": 5000},
        "loss": {
            "type": "ener",
            "start_pref_e": 0.02,
            "limit_pref_e": 1.0,
            "start_pref_f": 1000.0,
            "limit_pref_f": 1.0,
            "start_pref_v": 0.0,
            "limit_pref_v": 0.0,
        },
        "training": {
            "training_data": {"systems": train_systems, "batch_size": 1, "auto_prob": "prob_uniform"},
            "validation_data": {
                "systems": valid_systems or train_systems,
                "batch_size": 1,
                "numb_batch": 10,
            },
            "numb_steps": steps,
            "seed": 1,
            "disp_file": "lcurve.out",
            "disp_freq": 1000,
            "save_freq": 10000,
            "save_ckpt": "model.ckpt",
        },
    }


def cmd_make_train(args) -> None:
    train_systems = [str(p.resolve()) for p in args.train_systems]
    valid_systems = [str(p.resolve()) for p in args.valid_systems]
    type_map = args.type_map or infer_type_map_from_data(args.train_systems + args.valid_systems)
    if not type_map:
        die("cannot infer type_map; pass --type-map C H O ...")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(deepmd_input(type_map, train_systems, valid_systems, args.steps), indent=2) + "\n")
    run = args.out.parent / "run_deepmd.sh"
    run.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

DP="${{DP:-{DEFAULT_DP}}}"
"$DP" train "{args.out.name}"
"$DP" freeze -o frozen_model.pb
"""
    )
    run.chmod(0o755)
    print(f"wrote {args.out} and {run}")


def cmd_init_workflow(args) -> None:
    root = args.out
    for sub in ["00_cif", "01_candidates", "02_abacus_sp", "03_deepmd_data", "04_train"]:
        (root / sub).mkdir(parents=True, exist_ok=True)
    guide = f"""# ABACUS -> DeepMD workflow

1. Put starting CIF files in `00_cif/`.
2. Generate perturbed candidate structures:
   `abacuskit make-candidates 00_cif/your.cif --out 01_candidates --count 50 --rattle 0.03 --strain 0.01`
3. Prepare ABACUS labeling jobs:
   `abacuskit prepare-abacus 01_candidates --out 02_abacus_sp --cal-stress --nspin 2`
   Default ABACUS: `{DEFAULT_ABACUS_LABEL}`
   To refresh launch scripts for existing jobs:
   `abacuskit launch-script 02_abacus_sp --array-script 02_abacus_sp/run_all_abacus.sh`
4. Run ABACUS:
   `bash 02_abacus_sp/run_all_abacus.sh`
5. Check and summarize ABACUS jobs:
   `abacuskit check-abacus 02_abacus_sp --json check_report.json --csv check_report.csv`
   `abacuskit collect-metrics 02_abacus_sp --json metrics.json --csv metrics.csv`
   `abacuskit report-metrics --metrics metrics.json --out abacuskit_report.html`
6. Optional convergence test from one prepared job:
   `abacuskit conv-test 02_abacus_sp/000000 --key ecutwfc --values 80 100 120 --out conv_ecutwfc`
7. Convert finished ABACUS outputs to DeepMD npy:
   `abacuskit collect-deepmd 02_abacus_sp --out 03_deepmd_data --split-ratio 0.1`
8. Create and run DeepMD training:
   `abacuskit make-train 03_deepmd_data/train/* --valid-systems 03_deepmd_data/valid/* --out 04_train/input.json`
   `cd 04_train && bash run_deepmd.sh`
"""
    (root / "README_workflow.md").write_text(guide)
    print(f"initialized workflow under {root}")


class MenuExit(Exception):
    pass


class ProgramExit(Exception):
    pass


def prompt_text(label: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    try:
        value = input(f"{label}{suffix}: ").strip()
    except EOFError as exc:
        print()
        raise MenuExit from exc
    if not value and default is not None:
        return default
    return value


def prompt_path(label: str, default: str | None = None) -> Path:
    return Path(prompt_text(label, default)).expanduser()


def prompt_int(label: str, default: int) -> int:
    while True:
        value = prompt_text(label, str(default))
        try:
            return int(value)
        except ValueError:
            print("Please enter an integer.")


def prompt_float(label: str, default: float) -> float:
    while True:
        value = prompt_text(label, str(default))
        try:
            return float(value)
        except ValueError:
            print("Please enter a number.")


def prompt_yes_no(label: str, default: bool = False) -> bool:
    default_text = "y" if default else "n"
    while True:
        value = prompt_text(label, default_text).lower()
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please enter y or n.")


def prompt_choice(label: str, choices: list[str], default: str) -> str:
    choice_text = "/".join(choices)
    while True:
        value = prompt_text(f"{label} ({choice_text})", default).lower()
        if value in choices:
            return value
        print(f"Please choose one of: {choice_text}")


def prompt_multi(label: str) -> list[str]:
    value = prompt_text(label, "")
    return value.split() if value else []


def choose_cif_from_current_dir() -> Path | None:
    cifs = sorted(Path.cwd().glob("*.cif"), key=lambda p: natural_key(p.name))
    if not cifs:
        return None
    print("\nCIF files in current directory:")
    for i, path in enumerate(cifs, start=1):
        print(f"  {i:2d}) {path.name}")
    value = prompt_text("Choose a CIF number, or press Enter to type a path", "")
    if not value:
        return None
    try:
        idx = int(value)
    except ValueError:
        return Path(value).expanduser()
    if 1 <= idx <= len(cifs):
        return cifs[idx - 1]
    print("Invalid number; please type the CIF path manually.")
    return None


def first_cif_if_unique() -> Path | None:
    cifs = sorted(Path.cwd().glob("*.cif"), key=lambda p: natural_key(p.name))
    return cifs[0] if len(cifs) == 1 else None


def interactive_structure_options() -> dict:
    orbital_quality = prompt_choice("Orbital library quality", ["efficiency", "precision"], "efficiency")
    basis_type = prompt_choice("Basis type", ["lcao", "pw"], "lcao")
    supercell_text = prompt_text("Supercell nx ny nz", "1 1 1").split()
    if len(supercell_text) != 3:
        print("Bad supercell input; using 1 1 1.")
        supercell = (1, 1, 1)
    else:
        supercell = tuple(int(x) for x in supercell_text)

    mag = None
    fix_element = None
    fix_index = None
    fix_below = None
    fix_above = None
    element_orbital_quality = None
    pseudo_file = None
    orbital_file = None
    pseudo_dir = DEFAULT_PSEUDO_DIR
    orbital_dir = None

    if prompt_yes_no("Use advanced structure options?", False):
        pseudo_dir = prompt_path("Pseudo directory", str(DEFAULT_PSEUDO_DIR))
        if basis_type == "lcao":
            orbital_dir_text = prompt_text("Orbital directory, empty means default library", "")
            orbital_dir = Path(orbital_dir_text).expanduser() if orbital_dir_text else None
        mag = prompt_multi("Magnetic moments, e.g. Ni=2 Fe=3, empty for none") or None
        fix_element = prompt_multi("Fix elements, e.g. Ni=z C=xy, empty for none") or None
        fix_index = prompt_multi("Fix atom indices, e.g. 1-10=z 15=xy, empty for none") or None
        fix_below = prompt_multi("Fix below cutoff, e.g. z=5.0:xyz, empty for none") or None
        fix_above = prompt_multi("Fix above cutoff, e.g. z=25.0:z, empty for none") or None
        element_orbital_quality = (
            prompt_multi("Element orbital quality, e.g. C=precision Ni=efficiency, empty for none")
            or None
        )
        pseudo_file = prompt_multi("Pseudo overrides, e.g. Ni=/path/Ni.upf, empty for none") or None
        orbital_file = prompt_multi("Orbital overrides, e.g. Ni=/path/Ni.orb, empty for none") or None

    return {
        "pseudo_dir": pseudo_dir,
        "orbital_dir": orbital_dir,
        "orbital_quality": orbital_quality,
        "element_orbital_quality": element_orbital_quality,
        "orbital_file": orbital_file,
        "pseudo_file": pseudo_file,
        "basis_type": basis_type,
        "supercell": supercell,
        "mag": mag,
        "fix_element": fix_element,
        "fix_index": fix_index,
        "fix_below": fix_below,
        "fix_above": fix_above,
    }


def default_cif2stru_state() -> dict:
    return {
        "cif": first_cif_if_unique(),
        "output": Path("STRU"),
        "pseudo_dir": DEFAULT_PSEUDO_DIR,
        "orbital_dir": None,
        "orbital_quality": "efficiency",
        "element_orbital_quality": None,
        "orbital_file": None,
        "pseudo_file": None,
        "basis_type": "lcao",
        "supercell": (1, 1, 1),
        "mag": None,
        "fix_element": None,
        "fix_index": None,
        "fix_below": None,
        "fix_above": None,
    }


def format_menu_value(value) -> str:
    if value is None or value == []:
        return "none"
    if isinstance(value, tuple):
        return " ".join(str(x) for x in value)
    if isinstance(value, list):
        return " ".join(str(x) for x in value) if value else "none"
    return str(value)


def print_cif2stru_menu(state: dict) -> None:
    print(
        f"""
---------- 10x: CIF -> ABACUS STRU ----------
Current settings:
  CIF file        : {format_menu_value(state["cif"])}
  Output          : {state["output"]}
  Orbital basis   : {state["orbital_quality"]} LCAO
  Supercell       : {format_menu_value(state["supercell"])}
  Magnetic moments: {format_menu_value(state["mag"])}
  Fixed elements  : {format_menu_value(state["fix_element"])}
  Fixed indices   : {format_menu_value(state["fix_index"])}
  Fixed cutoff    : below={format_menu_value(state["fix_below"])} above={format_menu_value(state["fix_above"])}

  101) Generate STRU now using current settings
  102) Set orbital basis to precision
  103) Set magnetic moment, e.g. Ni 2
  104) Set supercell, e.g. 2 2 1
  105) Fix by element, e.g. Ni z
  106) Fix by atom index, e.g. 1-10 z
  107) Fix by coordinate cutoff, e.g. below z 5.0 xyz
  108) Select/change CIF file
  109) Change output STRU path
  110) Reset to defaults
  111) Set orbital basis to efficiency
  112) Clear magnetic moments
  113) Clear fixed atom settings
  0) Back to previous menu
  q) Quit abacuskit
"""
    )


def append_option(state: dict, key: str, value: str) -> None:
    items = list(state.get(key) or [])
    items.append(value)
    state[key] = items


def set_extra_setting(state: dict, key: str, value: object) -> None:
    items = []
    for item in state.get("set") or []:
        old_key = item.split("=", 1)[0].strip() if "=" in item else item.split(None, 1)[0].strip()
        if old_key != key:
            items.append(item)
    items.append(f"{key}={value}")
    state["set"] = items


def remove_extra_setting(state: dict, *keys: str) -> None:
    key_set = set(keys)
    items = []
    for item in state.get("set") or []:
        old_key = item.split("=", 1)[0].strip() if "=" in item else item.split(None, 1)[0].strip()
        if old_key not in key_set:
            items.append(item)
    state["set"] = items or None


def get_extra_setting(state: dict, key: str) -> str | None:
    for item in state.get("set") or []:
        if "=" not in item:
            continue
        old_key, value = item.split("=", 1)
        if old_key.strip() == key:
            return value.strip()
    return None


def format_dftu_status(state: dict) -> str:
    mode = get_extra_setting(state, "dft_plus_u")
    if not mode or mode == "0":
        return "off"
    corr = get_extra_setting(state, "orbital_corr") or "-1"
    hubbard_u = get_extra_setting(state, "hubbard_u") or "0.0"
    return f"mode {mode}, orbital_corr={corr}, U={hubbard_u}"


def apply_dftu_settings(state: dict) -> None:
    mode = prompt_choice("DFT+U mode: 1 radius-adjustable projection, 2 first-zeta projection", ["1", "2"], "1")
    orbital_corr = prompt_text(
        "orbital_corr list, e.g. -1 2 -1 (1=p, 2=d, 3=f)",
        get_extra_setting(state, "orbital_corr") or "-1",
    )
    hubbard_u = prompt_text(
        "hubbard_u list in eV, e.g. 0 4.0 0",
        get_extra_setting(state, "hubbard_u") or "0.0",
    )
    set_extra_setting(state, "dft_plus_u", mode)
    set_extra_setting(state, "orbital_corr", orbital_corr)
    set_extra_setting(state, "hubbard_u", hubbard_u)
    set_extra_setting(state, "yukawa_potential", "false")
    set_extra_setting(state, "omc", 0)
    if mode == "1":
        set_extra_setting(state, "onsite_radius", get_extra_setting(state, "onsite_radius") or 3.0)
    else:
        remove_extra_setting(state, "onsite_radius")
    print("DFT+U enabled. Check atom-type order against STRU before running ABACUS.")


def apply_dftu_mixing_aid(state: dict) -> None:
    mixing_restart = prompt_text("mixing_restart for DFT+U restart", "5e-4")
    uramping = prompt_text("uramping in eV, -1 disables U-ramping", "-1.0")
    set_extra_setting(state, "mixing_restart", mixing_restart)
    set_extra_setting(state, "mixing_dmr", "true")
    set_extra_setting(state, "uramping", uramping)
    print("DFT+U convergence aid applied: mixing_restart, mixing_dmr, and uramping are explicit.")


def clear_dftu_settings(state: dict) -> None:
    remove_extra_setting(state, *DFTU_KEYS)
    print("DFT+U settings disabled.")


def clear_dftu_mixing_aid(state: dict) -> None:
    remove_extra_setting(state, *DFTU_MIXING_KEYS)
    print("DFT+U convergence-aid settings cleared.")


def ensure_cif_selected(state: dict) -> None:
    if state["cif"] and Path(state["cif"]).expanduser().is_file():
        state["cif"] = Path(state["cif"]).expanduser()
        return
    unique = first_cif_if_unique()
    if unique:
        state["cif"] = unique
        print(f"Using CIF file: {unique}")
        return
    cif = choose_cif_from_current_dir()
    if cif is None:
        cif = prompt_path("CIF file")
    state["cif"] = cif


def run_cif2stru_from_state(state: dict) -> None:
    ensure_cif_selected(state)
    args = argparse.Namespace(**state)
    cmd_cif2stru(args)


def parse_two_fields(text: str, example: str) -> tuple[str, str]:
    parts = text.replace("=", " ").split()
    if len(parts) != 2:
        die(f"expected two fields, for example: {example}")
    return parts[0], parts[1]


def interactive_cif2stru() -> None:
    state = default_cif2stru_state()
    while True:
        print_cif2stru_menu(state)
        choice = prompt_text("Enter 10x option", "101").lower()
        try:
            if choice in {"q", "quit", "exit"}:
                raise ProgramExit
            if choice in {"100", "0"}:
                return
            if choice == "101":
                run_cif2stru_from_state(state)
                print("STRU generation finished. Exiting abacuskit.")
                raise ProgramExit
            elif choice == "102":
                state["orbital_quality"] = "precision"
                print("Orbital basis set to precision.")
            elif choice == "103":
                sym, moment = parse_two_fields(
                    prompt_text("Element and initial magnetic moment, e.g. Ni 2"),
                    "Ni 2",
                )
                append_option(state, "mag", f"{sym}={moment}")
            elif choice == "104":
                values = prompt_text("Supercell nx ny nz, e.g. 2 2 1").split()
                if len(values) != 3:
                    die("expected three integers, for example: 2 2 1")
                state["supercell"] = tuple(int(x) for x in values)
            elif choice == "105":
                sym, dirs = parse_two_fields(
                    prompt_text("Element and fixed directions, e.g. Ni z"),
                    "Ni z",
                )
                append_option(state, "fix_element", f"{sym}={dirs}")
            elif choice == "106":
                selector, dirs = parse_two_fields(
                    prompt_text("Atom indices and fixed directions, e.g. 1-10 z"),
                    "1-10 z",
                )
                append_option(state, "fix_index", f"{selector}={dirs}")
            elif choice == "107":
                values = prompt_text("Cutoff rule: below/above axis value dirs, e.g. below z 5.0 xyz").split()
                if len(values) != 4 or values[0] not in {"below", "above"}:
                    die("expected: below z 5.0 xyz or above z 25.0 z")
                relation, axis, cutoff, dirs = values
                key = "fix_below" if relation == "below" else "fix_above"
                append_option(state, key, f"{axis}={cutoff}:{dirs}")
            elif choice == "108":
                cif = choose_cif_from_current_dir()
                state["cif"] = cif or prompt_path("CIF file")
            elif choice == "109":
                state["output"] = prompt_path("Output STRU file", "STRU")
            elif choice == "110":
                state.clear()
                state.update(default_cif2stru_state())
                print("CIF -> STRU settings reset to defaults.")
            elif choice == "111":
                state["orbital_quality"] = "efficiency"
                print("Orbital basis set to efficiency.")
            elif choice == "112":
                state["mag"] = None
                print("Magnetic moments cleared.")
            elif choice == "113":
                state["fix_element"] = None
                state["fix_index"] = None
                state["fix_below"] = None
                state["fix_above"] = None
                print("Fixed atom settings cleared.")
            else:
                print("Unknown 10x option.")
        except SystemExit as exc:
            print(exc)


def interactive_make_candidates() -> None:
    print("\n[2] Make candidate CIFs\n")
    cif = choose_cif_from_current_dir()
    if cif is None:
        cif = prompt_path("Seed CIF file")
    supercell_text = prompt_text("Supercell nx ny nz", "1 1 1").split()
    supercell = tuple(int(x) for x in supercell_text) if len(supercell_text) == 3 else (1, 1, 1)
    args = argparse.Namespace(
        cif=cif,
        out=prompt_path("Output directory", "01_candidates"),
        count=prompt_int("Number of structures", 20),
        rattle=prompt_float("Position noise in Angstrom", 0.03),
        strain=prompt_float("Cell strain amplitude", 0.0),
        supercell=supercell,
        seed=prompt_int("Random seed", 1),
    )
    cmd_make_candidates(args)


def interactive_prepare_abacus() -> None:
    print("\n[3] Prepare ABACUS jobs\n")
    cifs = [Path(x).expanduser() for x in prompt_multi("CIF file or directory paths")]
    if not cifs:
        cifs = [prompt_path("CIF file or directory path", "01_candidates")]
    opts = interactive_structure_options()
    args = argparse.Namespace(
        cifs=cifs,
        out=prompt_path("Output job directory", "02_abacus_sp"),
        suffix=prompt_text("ABACUS output suffix", "ABACUS"),
        calculation=prompt_choice("Calculation type", ["scf", "relax", "md"], "scf"),
        device=prompt_choice("Device", ["gpu", "cpu"], "gpu"),
        ks_solver=prompt_text("KS solver", "cusolver"),
        kspacing=prompt_float("K spacing", 0.14),
        ecutwfc=prompt_float("ecutwfc", 100),
        nspin=prompt_int("nspin", 1),
        cal_stress=prompt_yes_no("Calculate stress?", False),
        set=prompt_multi("Extra INPUT key=value items, empty for none") or None,
        abacus_env=prompt_path("ABACUS env script", str(DEFAULT_ABACUS_ENV)),
        mpi_np=prompt_int("MPI processes", 1),
        gpu_ids=prompt_text("GPU ids, empty for default", "") or None,
        omp_threads=prompt_int("OMP threads", 12),
        no_numactl=not prompt_yes_no("Use numactl CPU binding?", True),
        cpu_bind=prompt_text("CPU bind range, empty for auto", "") or None,
        mem_bind=prompt_text("NUMA memory node", "0"),
        **opts,
    )
    cmd_prepare_abacus(args)


def default_input_state() -> dict:
    return {
        "kind": "scf",
        "out": Path("INPUT"),
        "suffix": "ABACUS",
        "pseudo_dir": DEFAULT_PSEUDO_DIR,
        "orbital_dir": None,
        "orbital_quality": "efficiency",
        "basis_type": "lcao",
        "device": "gpu",
        "ks_solver": "cusolver",
        "kspacing": 0.14,
        "ecutwfc": 100,
        "nspin": 1,
        "cal_stress": False,
        "relax_nmax": 100,
        "force_thr_ev": 0.04,
        "stress_thr": 1.0,
        "dos": False,
        "set": None,
        "no_comments": False,
    }


def print_input_menu(state: dict) -> None:
    vdw = get_extra_setting(state, "vdw_method") or "off"
    dipole = get_extra_setting(state, "dip_cor_flag") or "off"
    dftu = format_dftu_status(state)
    print(
        f"""
---------- 30x: Generate ABACUS INPUT ----------
Current settings:
  Output          : {state["out"]}
  Calculation     : {state["kind"]}
  Suffix          : {state["suffix"]}
  Orbital basis   : {state["orbital_quality"]} LCAO
  Device / solver : {state["device"]} / {state["ks_solver"]}
  kspacing        : {state["kspacing"]}
  ecutwfc         : {state["ecutwfc"]}
  nspin           : {state["nspin"]}
  cal_stress      : {state["cal_stress"]}
  DOS/PDOS        : {state["dos"]}
  VDW             : {vdw}
  Dipole corr.    : {dipole}
  DFT+U           : {dftu}
  Extra INPUT     : {format_menu_value(state["set"])}

  301) Generate INPUT now using current settings
  302) Set calculation to scf
  303) Set calculation to relax
  304) Set orbital basis to precision
  305) Set orbital basis to efficiency
  306) Set nspin
  307) Set ecutwfc
  308) Set kspacing
  309) Set device, gpu or cpu
  310) Set ks_solver
  311) Toggle cal_stress
  312) Toggle DOS/PDOS output
  313) Add extra INPUT key=value
  314) Change output INPUT path
  315) Change suffix
  316) Set relax parameters
  317) Clear extra INPUT settings
  318) Reset to defaults
  319) Set VDW correction, e.g. d3_bj
  320) Toggle dipole correction, default Z axis
  321) Apply DOS target template
  322) Apply PDOS target template
  323) Apply band structure target template
  324) Apply COHP matrix-output template
  325) Apply work-function/potential template
  326) Enable/edit DFT+U
  327) Apply DFT+U convergence-aid template
  328) Disable DFT+U
  329) Clear DFT+U convergence-aid settings
  330) Apply ELF cube-output template
  331) Apply charge-density cube-output template
  0) Back to previous menu
  q) Quit abacuskit
"""
    )


def run_input_from_state(state: dict) -> None:
    args = argparse.Namespace(**state)
    cmd_input_template(args)


def apply_input_target_template(state: dict, target: str) -> None:
    if target == "dos":
        state["kind"] = "nscf"
        state["dos"] = True
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_dos", 1)
        set_extra_setting(state, "dos_sigma", 0.07)
        set_extra_setting(state, "dos_edelta_ev", 0.01)
        print("DOS target template applied. Prepare a dense KPT mesh and previous charge density.")
    elif target == "pdos":
        state["kind"] = "nscf"
        state["dos"] = True
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_dos", 2)
        set_extra_setting(state, "dos_sigma", 0.07)
        set_extra_setting(state, "dos_edelta_ev", 0.01)
        print("PDOS target template applied. LCAO basis will output the PDOS file.")
    elif target == "band":
        state["kind"] = "nscf"
        state["dos"] = False
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_band", 1)
        set_extra_setting(state, "out_proj_band", 1)
        set_extra_setting(state, "smearing_method", "gaussian")
        set_extra_setting(state, "smearing_sigma", 0.02)
        print("Band target template applied. Prepare a line-mode KPT file before running ABACUS.")
    elif target == "cohp":
        state["kind"] = "scf"
        set_extra_setting(state, "out_mat_hs2", 1)
        set_extra_setting(state, "out_mat_hs", 1)
        set_extra_setting(state, "out_app_flag", "false")
        print("COHP matrix-output template applied. Use generated H/S matrices for COHP post-processing.")
    elif target == "workfunc":
        state["kind"] = "scf"
        set_extra_setting(state, "out_pot", 2)
        set_extra_setting(state, "efield_flag", "true")
        set_extra_setting(state, "dip_cor_flag", "true")
        set_extra_setting(state, "efield_dir", 2)
        set_extra_setting(state, "efield_amp", 0)
        print("Work-function/potential template applied. Default dipole correction direction is Z.")
    elif target == "elf":
        state["kind"] = "scf"
        set_extra_setting(state, "out_elf", "1 3")
        print("ELF cube-output template applied. ABACUS will write elf.cube under OUT.<suffix>.")
    elif target == "charge":
        state["kind"] = "scf"
        set_extra_setting(state, "out_chg", "1 3")
        print("Charge-density cube-output template applied. ABACUS will write charge-density cube files.")
    else:
        die(f"unknown INPUT target template: {target}")


def interactive_input_template() -> None:
    state = default_input_state()
    while True:
        print_input_menu(state)
        choice = prompt_text("Enter 30x option", "301").lower()
        try:
            if choice in {"q", "quit", "exit"}:
                raise ProgramExit
            if choice in {"0", "300"}:
                return
            if choice == "301":
                run_input_from_state(state)
                print("INPUT generation finished. Exiting abacuskit.")
                raise ProgramExit
            if choice == "302":
                state["kind"] = "scf"
                print("Calculation set to scf.")
            elif choice == "303":
                state["kind"] = "relax"
                print("Calculation set to relax.")
            elif choice == "304":
                state["orbital_quality"] = "precision"
                print("Orbital basis set to precision.")
            elif choice == "305":
                state["orbital_quality"] = "efficiency"
                print("Orbital basis set to efficiency.")
            elif choice == "306":
                state["nspin"] = prompt_int("nspin", state["nspin"])
            elif choice == "307":
                state["ecutwfc"] = prompt_float("ecutwfc", state["ecutwfc"])
            elif choice == "308":
                state["kspacing"] = prompt_float("kspacing", state["kspacing"])
            elif choice == "309":
                state["device"] = prompt_choice("Device", ["gpu", "cpu"], state["device"])
            elif choice == "310":
                state["ks_solver"] = prompt_text("KS solver", state["ks_solver"])
            elif choice == "311":
                state["cal_stress"] = not state["cal_stress"]
                print(f"cal_stress set to {state['cal_stress']}.")
            elif choice == "312":
                state["dos"] = not state["dos"]
                print(f"DOS/PDOS output set to {state['dos']}.")
            elif choice == "313":
                item = prompt_text("Extra INPUT key=value, e.g. mixing_beta=0.05")
                if "=" not in item:
                    die("expected key=value, for example: mixing_beta=0.05")
                append_option(state, "set", item)
            elif choice == "314":
                state["out"] = prompt_path("Output INPUT file", "INPUT")
            elif choice == "315":
                state["suffix"] = prompt_text("ABACUS output suffix", state["suffix"])
            elif choice == "316":
                state["kind"] = "relax"
                state["relax_nmax"] = prompt_int("relax_nmax", state["relax_nmax"])
                state["force_thr_ev"] = prompt_float("force_thr_ev", state["force_thr_ev"])
                state["stress_thr"] = prompt_float("stress_thr", state["stress_thr"])
            elif choice == "317":
                state["set"] = None
                print("Extra INPUT settings cleared.")
            elif choice == "318":
                state.clear()
                state.update(default_input_state())
                print("INPUT settings reset to defaults.")
            elif choice == "319":
                method = prompt_choice("VDW method", ["d3_bj", "d3_0", "d2", "none"], "d3_bj")
                if method == "none":
                    remove_extra_setting(state, "vdw_method")
                    print("VDW correction disabled.")
                else:
                    set_extra_setting(state, "vdw_method", method)
                    print(f"VDW correction set to {method}.")
            elif choice == "320":
                if get_extra_setting(state, "dip_cor_flag"):
                    remove_extra_setting(
                        state,
                        "efield_flag",
                        "dip_cor_flag",
                        "efield_dir",
                        "efield_amp",
                        "efield_pos_max",
                        "efield_pos_dec",
                    )
                    print("Dipole correction disabled.")
                else:
                    axis = prompt_choice("Dipole correction axis", ["z", "x", "y"], "z")
                    axis_to_dir = {"x": 0, "y": 1, "z": 2}
                    set_extra_setting(state, "efield_flag", "true")
                    set_extra_setting(state, "dip_cor_flag", "true")
                    set_extra_setting(state, "efield_dir", axis_to_dir[axis])
                    set_extra_setting(state, "efield_amp", 0)
                    print(f"Dipole correction enabled along {axis.upper()} with efield_amp=0.")
            elif choice == "321":
                apply_input_target_template(state, "dos")
            elif choice == "322":
                apply_input_target_template(state, "pdos")
            elif choice == "323":
                apply_input_target_template(state, "band")
            elif choice == "324":
                apply_input_target_template(state, "cohp")
            elif choice == "325":
                apply_input_target_template(state, "workfunc")
            elif choice == "326":
                apply_dftu_settings(state)
            elif choice == "327":
                apply_dftu_mixing_aid(state)
            elif choice == "328":
                clear_dftu_settings(state)
            elif choice == "329":
                clear_dftu_mixing_aid(state)
            elif choice == "330":
                apply_input_target_template(state, "elf")
            elif choice == "331":
                apply_input_target_template(state, "charge")
            else:
                print("Unknown 30x option.")
        except SystemExit as exc:
            print(exc)


def interactive_check_abacus() -> None:
    print("\n[5] Check ABACUS job status in current directory\n")
    args = argparse.Namespace(jobs=[Path(".")], json=None, csv=None)
    cmd_check_abacus(args)
    print("ABACUS job check finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_launch_script() -> None:
    print("\n[13] Create ABACUS launch scripts\n")
    jobs = [Path(x).expanduser() for x in prompt_multi("ABACUS job directory paths")]
    if not jobs:
        jobs = [prompt_path("ABACUS job directory path", ".")]
    args = argparse.Namespace(
        jobs=jobs,
        abacus_env=prompt_path("ABACUS env script", str(DEFAULT_ABACUS_ENV)),
        mpi_np=prompt_int("MPI processes", 1),
        gpu_ids=prompt_text("GPU ids", "0") or "0",
        omp_threads=prompt_int("OMP threads", 12),
        no_numactl=not prompt_yes_no("Use numactl CPU binding?", True),
        cpu_bind=prompt_text("CPU bind range, empty for auto", "") or None,
        mem_bind=prompt_text("NUMA memory node", "0"),
        script_name=prompt_text("Per-job script name", "run_abacus.sh"),
        array_script=prompt_path("Launcher script path", "run_all_abacus.sh"),
    )
    cmd_launch_script(args)


def interactive_kpt() -> None:
    print("\n[10] Generate ABACUS KPT\n")
    mesh_text = prompt_text("K mesh nx ny nz", "3 3 1").replace(",", " ").split()
    if len(mesh_text) != 3:
        die("expected three integers, for example: 3 3 1")
    shift_text = prompt_text("K shift sx sy sz", "0 0 0").replace(",", " ").split()
    if len(shift_text) != 3:
        die("expected three integers, for example: 0 0 0")
    args = argparse.Namespace(
        mesh=[int(x) for x in mesh_text],
        shift=[int(x) for x in shift_text],
        model=prompt_choice("KPT model", ["gamma", "mp"], "gamma"),
        out=prompt_path("Output KPT file", "KPT"),
    )
    cmd_kpt(args)


def interactive_conv_test() -> None:
    print("\n[11] Prepare convergence-test jobs\n")
    jobs = [Path(x).expanduser() for x in prompt_multi("Template ABACUS job paths")]
    if not jobs:
        jobs = [prompt_path("Template ABACUS job path", ".")]
    key = prompt_text("INPUT key to sweep, e.g. ecutwfc/kspacing/kpt", "ecutwfc")
    values = prompt_multi("Values, e.g. 80 100 120 or '2 2 1','3 3 1'")
    if not values:
        die("at least one value is required")
    args = argparse.Namespace(
        jobs=jobs,
        key=key,
        values=values,
        out=prompt_path("Output convergence directory", f"conv_{key}"),
        kpt_model=prompt_choice("KPT model if key=kpt", ["gamma", "mp"], "gamma"),
        force=prompt_yes_no("Overwrite existing target jobs?", False),
    )
    cmd_conv_test(args)


def interactive_collect_report() -> None:
    print("\n[12] Collect ABACUS metrics / report\n")
    jobs = [Path(x).expanduser() for x in prompt_multi("ABACUS job directory paths")]
    if not jobs:
        jobs = [prompt_path("ABACUS job directory path", "02_abacus_sp")]
    json_path = prompt_path("Metrics JSON output", "metrics.json")
    csv_path = prompt_path("Metrics CSV output", "metrics.csv")
    args = argparse.Namespace(jobs=jobs, json=json_path, csv=csv_path, metrics=None)
    cmd_collect_metrics(args)
    if prompt_yes_no("Generate HTML report?", True):
        report_args = argparse.Namespace(metrics=json_path, out=prompt_path("HTML report output", "abacuskit_report.html"), keys=None)
        cmd_report_metrics(report_args)


def interactive_plot_dos() -> None:
    print("\n[6] Plot DOS / PDOS / LDOS\n")
    args = argparse.Namespace(
        path=prompt_path("ABACUS job or OUT.* directory"),
        kind=prompt_choice("Plot kind", ["auto", "dos", "pdos", "ldos"], "auto"),
        file=None,
        select=prompt_multi("PDOS selectors, e.g. C=p H=s Ni=d, empty for none") or None,
        fermi=None,
        out=prompt_path("Output image", "dos.png"),
    )
    fermi_text = prompt_text("Fermi energy shift in eV, empty for none", "")
    args.fermi = float(fermi_text) if fermi_text else None
    file_text = prompt_text("Explicit data file, empty for auto", "")
    args.file = Path(file_text).expanduser() if file_text else None
    cmd_plot_dos(args)


def interactive_plot_grid() -> None:
    print("\n[14] Plot ELF / charge density / charge-density difference\n")
    kind = prompt_choice("Plot kind", ["auto", "elf", "charge", "diff", "cube"], "auto")
    args = argparse.Namespace(
        path=prompt_path("ABACUS job, OUT.* directory, or cube file"),
        kind=kind,
        file=None,
        minus=None,
        minus_file=None,
        cube_out=None,
        axis=prompt_choice("Slice axis", ["z", "x", "y"], "z"),
        index=None,
        cmap=None,
        out=None,
    )
    index_text = prompt_text("Slice index, empty for middle", "")
    args.index = int(index_text) if index_text else None
    file_text = prompt_text("Explicit cube file, empty for auto", "")
    args.file = Path(file_text).expanduser() if file_text else None
    if kind == "diff":
        args.minus = prompt_path("Subtracted ABACUS job, OUT.* directory, or cube file")
        minus_file_text = prompt_text("Explicit subtracted cube file, empty for auto", "")
        args.minus_file = Path(minus_file_text).expanduser() if minus_file_text else None
        cube_out_text = prompt_text("Output difference cube", "charge_diff.cube")
        args.cube_out = Path(cube_out_text).expanduser() if cube_out_text else None
    default_out = "charge_diff.png" if kind == "diff" else f"{kind if kind != 'auto' else 'grid'}.png"
    args.out = prompt_path("Output image", default_out)
    cmd_plot_grid(args)


def interactive_collect_deepmd() -> None:
    print("\n[7] Collect ABACUS outputs to DeepMD data\n")
    jobs = [Path(x).expanduser() for x in prompt_multi("ABACUS job directory paths")]
    if not jobs:
        jobs = [prompt_path("ABACUS job directory path", "02_abacus_sp")]
    fmt = prompt_choice("dpdata format", ["auto", "abacus/scf", "abacus/md", "abacus/relax"], "auto")
    args = argparse.Namespace(
        jobs=jobs,
        out=prompt_path("Output DeepMD data directory", "03_deepmd_data"),
        fmt=None if fmt == "auto" else fmt,
        set_size=prompt_int("DeepMD set size", 5000),
        split_ratio=prompt_float("Validation split ratio", 0.1),
        seed=prompt_int("Random seed", 1),
        python=DEFAULT_DEEPMD_PYTHON,
    )
    cmd_collect_deepmd(args)


def interactive_make_train() -> None:
    print("\n[8] Make DeepMD training input\n")
    train_systems = [Path(x).expanduser() for x in prompt_multi("Training system paths")]
    valid_systems = [Path(x).expanduser() for x in prompt_multi("Validation system paths, empty for none")]
    type_map = prompt_multi("Type map, e.g. C H O Ni, empty to infer") or None
    args = argparse.Namespace(
        train_systems=train_systems,
        valid_systems=valid_systems,
        type_map=type_map,
        steps=prompt_int("Training steps", 100000),
        out=prompt_path("Output input.json", "04_train/input.json"),
    )
    cmd_make_train(args)


def interactive_init_workflow() -> None:
    print("\n[9] Init workflow skeleton\n")
    args = argparse.Namespace(out=prompt_path("Workflow root directory", "abacus_deepmd_project"))
    cmd_init_workflow(args)


def print_interactive_menu() -> None:
    print_terminal_logo()
    print(
        f"""
============== abacuskit {__version__} ==============
Author: {__author__}
Affiliation: {__affiliation__}

  1) CIF -> ABACUS STRU
  2) Make candidate CIFs
  3) Generate ABACUS INPUT
  4) Prepare ABACUS jobs
  5) Check ABACUS job status
  6) Plot DOS / PDOS / LDOS
  7) Collect ABACUS outputs to DeepMD data
  8) Make DeepMD training input
  9) Init workflow skeleton
  10) Generate ABACUS KPT
  11) Prepare convergence-test jobs
  12) Collect ABACUS metrics / report
  13) Create ABACUS launch scripts
  14) Plot ELF / charge density / charge-density difference
  h) Show command-line help
  q) Quit abacuskit
  0) Exit
"""
    )


def interactive_menu() -> None:
    actions = {
        "1": interactive_cif2stru,
        "2": interactive_make_candidates,
        "3": interactive_input_template,
        "4": interactive_prepare_abacus,
        "5": interactive_check_abacus,
        "6": interactive_plot_dos,
        "7": interactive_collect_deepmd,
        "8": interactive_make_train,
        "9": interactive_init_workflow,
        "10": interactive_kpt,
        "11": interactive_conv_test,
        "12": interactive_collect_report,
        "13": interactive_launch_script,
        "14": interactive_plot_grid,
    }
    while True:
        print_interactive_menu()
        choice = prompt_text("Enter option number", "1").lower()
        if choice in {"q", "quit", "exit", "0"}:
            print("Bye.")
            return
        if choice in {"h", "help"}:
            build_parser().print_help()
        elif choice in actions:
            try:
                actions[choice]()
            except MenuExit:
                return
            except ProgramExit:
                print("Bye.")
                return
            except SystemExit as exc:
                print(exc)
            except KeyboardInterrupt:
                print("\nCancelled.")
        else:
            print("Unknown option.")


def cmd_menu(args) -> None:
    interactive_menu()


def add_common_structure_args(parser) -> None:
    parser.add_argument("--pseudo-dir", type=Path, default=DEFAULT_PSEUDO_DIR)
    parser.add_argument("--orbital-dir", type=Path)
    parser.add_argument("--orbital-quality", choices=sorted(DEFAULT_ORBITAL_DIRS), default="efficiency")
    parser.add_argument(
        "--element-orbital-quality",
        action="append",
        help="per-element APNS orbital library choice, e.g. --element-orbital-quality C=precision",
    )
    parser.add_argument("--orbital-file", action="append", help="exact orbital override, e.g. Ni=/path/Ni.orb")
    parser.add_argument("--pseudo-file", action="append", help="exact pseudopotential override, e.g. Ni=/path/Ni.upf")
    parser.add_argument("--basis-type", choices=["lcao", "pw"], default="lcao")
    parser.add_argument("--supercell", type=int, nargs=3, default=(1, 1, 1))
    parser.add_argument("--mag", action="append", help="initial magnetic moment, e.g. --mag Fe=2 --mag Ni=2")
    parser.add_argument("--fix-element", action="append", help="fix all atoms of an element along directions, e.g. C=z")
    parser.add_argument("--fix-index", action="append", help="fix 1-based atom indices, e.g. 1-10,15=xy")
    parser.add_argument("--fix-below", action="append", help="fix atoms below coordinate cutoff, e.g. z=3.0:xy")
    parser.add_argument("--fix-above", action="append", help="fix atoms above coordinate cutoff, e.g. z=20.0:z")


def build_parser() -> argparse.ArgumentParser:
    epilog = f"Version: {__version__}\nAuthor: {__author__}, {__affiliation__}"
    parser = argparse.ArgumentParser(
        prog="abacuskit",
        description=__doc__,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"abacuskit {__version__} ({__author__}, {__affiliation__})",
    )
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("menu", help="start the interactive numbered menu")
    p.set_defaults(func=cmd_menu)

    p = sub.add_parser("abacus-versions", help="list local ABACUS installations")
    p.add_argument("--root", type=Path, default=DEFAULT_ABACUS_ROOT)
    p.add_argument("--json", type=Path, help="write version list to JSON")
    p.set_defaults(func=cmd_abacus_versions)

    p = sub.add_parser("cif2stru", help="convert a CIF file to ABACUS STRU")
    p.add_argument("cif", type=Path)
    p.add_argument("-o", "--output", type=Path, default=Path("STRU"))
    add_common_structure_args(p)
    p.set_defaults(func=cmd_cif2stru)

    p = sub.add_parser("make-candidates", help="make randomly displaced/strained CIFs")
    p.add_argument("cif", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--count", type=int, default=20)
    p.add_argument("--rattle", type=float, default=0.03, help="Gaussian position noise in Angstrom")
    p.add_argument("--strain", type=float, default=0.0, help="symmetric random cell strain amplitude")
    p.add_argument("--supercell", type=int, nargs=3, default=(1, 1, 1))
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(func=cmd_make_candidates)

    p = sub.add_parser("prepare-abacus", help="write STRU/INPUT/run scripts for CIF files")
    p.add_argument("cifs", type=Path, nargs="+")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--suffix", default="ABACUS")
    p.add_argument("--calculation", choices=["scf", "relax", "md"], default="scf")
    p.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    p.add_argument("--ks-solver", default="cusolver")
    p.add_argument("--kspacing", type=float, default=0.14)
    p.add_argument("--ecutwfc", type=float, default=100)
    p.add_argument("--nspin", type=int, default=1)
    p.add_argument("--cal-stress", action="store_true")
    p.add_argument("--set", action="append", help="extra INPUT key=value; can be repeated")
    p.add_argument("--abacus-env", type=Path, default=DEFAULT_ABACUS_ENV)
    p.add_argument("--mpi-np", type=int, default=1)
    p.add_argument("--gpu-ids", default=None, help="for example 0 or 0,1")
    p.add_argument("--omp-threads", type=int, default=12)
    p.add_argument("--no-numactl", action="store_true", help="do not wrap mpirun with numactl")
    p.add_argument("--cpu-bind", help="numactl CPU range, default depends on MPI ranks")
    p.add_argument("--mem-bind", default="0", help="numactl memory node, default 0")
    add_common_structure_args(p)
    p.set_defaults(func=cmd_prepare_abacus)

    p = sub.add_parser("input-template", help="write an ABACUS INPUT template for scf or relax")
    p.add_argument("--kind", choices=["scf", "relax"], required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--suffix", default="ABACUS")
    p.add_argument("--pseudo-dir", type=Path, default=DEFAULT_PSEUDO_DIR)
    p.add_argument("--orbital-dir", type=Path)
    p.add_argument("--orbital-quality", choices=sorted(DEFAULT_ORBITAL_DIRS), default="efficiency")
    p.add_argument("--basis-type", choices=["lcao", "pw"], default="lcao")
    p.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    p.add_argument("--ks-solver", default="cusolver")
    p.add_argument("--kspacing", type=float, default=0.14)
    p.add_argument("--ecutwfc", type=float, default=100)
    p.add_argument("--nspin", type=int, default=1)
    p.add_argument("--cal-stress", action="store_true")
    p.add_argument("--relax-nmax", type=int, default=100)
    p.add_argument("--force-thr-ev", type=float, default=0.04)
    p.add_argument("--stress-thr", type=float, default=1.0)
    p.add_argument("--dos", action="store_true", help="include DOS/PDOS output parameters")
    p.add_argument("--set", action="append", help="extra INPUT key=value; can be repeated")
    p.add_argument("--no-comments", action="store_true")
    p.set_defaults(func=cmd_input_template)

    p = sub.add_parser("kpt", help="write an ABACUS KPT file")
    p.add_argument("--mesh", type=int, nargs=3, required=True, metavar=("NX", "NY", "NZ"))
    p.add_argument("--shift", type=int, nargs=3, default=(0, 0, 0), metavar=("SX", "SY", "SZ"))
    p.add_argument("--model", choices=["gamma", "mp"], default="gamma")
    p.add_argument("--out", type=Path, default=Path("KPT"))
    p.set_defaults(func=cmd_kpt)

    p = sub.add_parser("check-abacus", help="check whether ABACUS jobs finished and converged")
    p.add_argument("jobs", type=Path, nargs="+")
    p.add_argument("--json", type=Path, help="write full JSON report")
    p.add_argument("--csv", type=Path, help="write CSV report")
    p.set_defaults(func=cmd_check_abacus)

    p = sub.add_parser("launch-script", help="create GPU ABACUS launch scripts for existing jobs")
    p.add_argument("jobs", type=Path, nargs="+", help="ABACUS job directories or a root containing jobs")
    p.add_argument("--abacus-env", type=Path, default=DEFAULT_ABACUS_ENV)
    p.add_argument("--mpi-np", type=int, default=1)
    p.add_argument("--gpu-ids", default="0", help="CUDA_VISIBLE_DEVICES, default 0")
    p.add_argument("--omp-threads", type=int, default=12)
    p.add_argument("--no-numactl", action="store_true", help="do not wrap mpirun with numactl")
    p.add_argument("--cpu-bind", help="numactl CPU range, default depends on MPI ranks")
    p.add_argument("--mem-bind", default="0", help="numactl memory node, default 0")
    p.add_argument("--script-name", default="run_abacus.sh")
    p.add_argument("--array-script", type=Path, default=Path("run_all_abacus.sh"))
    p.set_defaults(func=cmd_launch_script)

    p = sub.add_parser("collect-metrics", help="collect ABACUS metrics similar to abacustest collectdata")
    p.add_argument("jobs", type=Path, nargs="+")
    p.add_argument("--json", type=Path, default=Path("metrics.json"))
    p.add_argument("--csv", type=Path)
    p.add_argument("--metrics", nargs="+", help="columns to print/write to CSV")
    p.set_defaults(func=cmd_collect_metrics)

    p = sub.add_parser("report-metrics", help="write a simple HTML report from collect-metrics JSON")
    p.add_argument("--metrics", type=Path, required=True, help="metrics JSON from collect-metrics")
    p.add_argument("--out", type=Path, default=Path("abacuskit_report.html"))
    p.add_argument("--keys", nargs="+", help="columns to include in the HTML table")
    p.set_defaults(func=cmd_report_metrics)

    p = sub.add_parser("conv-test", help="prepare convergence-test jobs by sweeping one INPUT key")
    p.add_argument("jobs", type=Path, nargs="+", help="template ABACUS job directories")
    p.add_argument("--key", required=True, help="INPUT key to sweep, or kpt to sweep KPT mesh")
    p.add_argument("--values", nargs="+", required=True, help="values to sweep; quote KPT meshes like '3 3 1'")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--kpt-model", choices=["gamma", "mp"], default="gamma")
    p.add_argument("--force", action="store_true", help="overwrite existing generated jobs")
    p.set_defaults(func=cmd_conv_test)

    p = sub.add_parser("plot-dos", help="plot ABACUS DOS, PDOS, or LDOS output")
    p.add_argument("path", type=Path, help="ABACUS job directory, OUT.* directory, or data file")
    p.add_argument("--kind", choices=["auto", "dos", "pdos", "ldos"], default="auto")
    p.add_argument("--file", type=Path, help="explicit DOS/PDOS/LDOS file")
    p.add_argument("--select", action="append", help="PDOS selector, e.g. C=p --select H=s --select O=p --select Ni=d")
    p.add_argument("--fermi", type=float, help="shift energy by this Fermi energy in eV")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_plot_dos)

    p = sub.add_parser("plot-grid", help="plot ELF/charge-density cube files or charge-density differences")
    p.add_argument("path", type=Path, help="ABACUS job directory, OUT.* directory, or cube file")
    p.add_argument("--kind", choices=["auto", "elf", "charge", "diff", "cube"], default="auto")
    p.add_argument("--file", type=Path, help="explicit cube file for ELF/charge or positive diff term")
    p.add_argument("--minus", type=Path, help="ABACUS job directory, OUT.* directory, or cube file to subtract")
    p.add_argument("--minus-file", type=Path, help="explicit cube file to subtract")
    p.add_argument("--cube-out", type=Path, help="write charge-density difference cube")
    p.add_argument("--axis", choices=["x", "y", "z"], default="z")
    p.add_argument("--index", type=int, help="slice index; default is the middle slice")
    p.add_argument("--cmap", help="matplotlib colormap name")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_plot_grid)

    p = sub.add_parser("collect-deepmd", help="convert ABACUS outputs to DeepMD npy via dpdata")
    p.add_argument("jobs", type=Path, nargs="+")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--fmt", choices=["abacus/scf", "abacus/md", "abacus/relax"])
    p.add_argument("--set-size", type=int, default=5000)
    p.add_argument("--split-ratio", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--python", type=Path, default=DEFAULT_DEEPMD_PYTHON)
    p.set_defaults(func=cmd_collect_deepmd)

    p = sub.add_parser("make-train", help="write a DeepMD input.json and run script")
    p.add_argument("train_systems", type=Path, nargs="+")
    p.add_argument("--valid-systems", type=Path, nargs="*", default=[])
    p.add_argument("--type-map", nargs="+")
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_make_train)

    p = sub.add_parser("init-workflow", help="create a standard workflow directory skeleton")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_init_workflow)
    return parser


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        interactive_menu()
        return 0
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
