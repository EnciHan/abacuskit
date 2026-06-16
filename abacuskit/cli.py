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
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

try:
    from . import __affiliation__, __author__, __version__
    from .bader import run_bader_analysis, write_bader_csv, write_bader_json
    from .cohp import build_orbital_map, format_orbital_map, resolve_orbital_arguments, run_cohp
except ImportError:
    __version__ = "v1.2.4"
    __author__ = "Han Enci, Zhong Lisheng, Yu Yutong, Xu Mengting, Chen Jingyuan"
    __affiliation__ = "Xi'an University of Technology"
    from bader import run_bader_analysis, write_bader_csv, write_bader_json
    from cohp import build_orbital_map, format_orbital_map, resolve_orbital_arguments, run_cohp

BOHR_PER_ANGSTROM = 1.88972612546
USER_CONFIG_PATH = Path.home() / ".abacuskit" / "config.json"


def env_path(name: str, default: str | Path) -> Path:
    return Path(os.environ.get(name, str(default))).expanduser()


APNS_RESOURCE_NAMES = {
    "pseudo": "apns-pseudopotentials-v1",
    "orbital_efficiency": "apns-orbitals-efficiency-v1",
    "orbital_precision": "apns-orbitals-precision-v1",
}
RESOURCE_CONFIG_KEYS = {
    "ABACUSKIT_PSEUDO_DIR": "pseudo_dir",
    "ABACUSKIT_ORBITAL_EFFICIENCY_DIR": "orbital_efficiency_dir",
    "ABACUSKIT_ORBITAL_PRECISION_DIR": "orbital_precision_dir",
}


def load_user_config() -> dict[str, object]:
    if not USER_CONFIG_PATH.is_file():
        return {}
    try:
        data = json.loads(USER_CONFIG_PATH.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


USER_CONFIG = load_user_config()


def unique_existing_dirs(paths: Iterable[Path]) -> list[Path]:
    result: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        expanded = path.expanduser()
        try:
            resolved = expanded.resolve()
        except OSError:
            resolved = expanded.absolute()
        if resolved in seen or not expanded.is_dir():
            continue
        seen.add(resolved)
        result.append(expanded)
    return result


def apns_search_roots() -> list[Path]:
    env_roots = os.environ.get("ABACUSKIT_APNS_SEARCH_ROOTS") or os.environ.get("ABACUSKIT_APNS_ROOT")
    roots: list[Path] = []
    if env_roots:
        roots.extend(Path(item) for item in env_roots.split(os.pathsep) if item)
    home = Path.home()
    cwd = Path.cwd()
    roots.extend(
        [
            cwd,
            cwd.parent,
            home / "data" / "abacus-lib",
            home / "data",
            home / "apps" / "abacus-lib",
            home / "apps",
            home,
            Path("/opt"),
            Path("/usr/local/share"),
        ]
    )
    return unique_existing_dirs(roots)


def direct_apns_candidates(name: str, roots: list[Path]) -> list[Path]:
    candidates: list[Path] = []
    for root in roots:
        candidates.extend(
            [
                root / name,
                root / "abacus-lib" / name,
                root / "data" / "abacus-lib" / name,
                root / "apps" / "abacus-lib" / name,
            ]
        )
    return candidates


def walk_find_apns_dir(name: str, roots: list[Path], max_depth: int = 5) -> Path | None:
    skip_names = {
        ".cache",
        ".conda",
        ".git",
        ".local",
        "__pycache__",
        "node_modules",
    }
    for root in roots:
        base_depth = len(root.parts)
        try:
            walker = os.walk(root)
            for current, dirs, _ in walker:
                current_path = Path(current)
                depth = len(current_path.parts) - base_depth
                if depth >= max_depth:
                    dirs[:] = []
                else:
                    dirs[:] = [item for item in dirs if item not in skip_names and not item.startswith(".")]
                if current_path.name == name:
                    return current_path
        except OSError:
            continue
    return None


def find_apns_dir(name: str) -> Path | None:
    roots = apns_search_roots()
    for candidate in direct_apns_candidates(name, roots):
        if candidate.is_dir():
            return candidate
    return walk_find_apns_dir(name, roots)


def resource_path(env_name: str, apns_name: str, fallback: str | Path) -> Path:
    if env_name in os.environ:
        return Path(os.environ[env_name]).expanduser()
    config_key = RESOURCE_CONFIG_KEYS.get(env_name)
    if config_key:
        configured = USER_CONFIG.get(config_key)
        if isinstance(configured, str) and configured:
            return Path(configured).expanduser()
    return find_apns_dir(apns_name) or Path(fallback).expanduser()


DEFAULT_PSEUDO_DIR = resource_path("ABACUSKIT_PSEUDO_DIR", APNS_RESOURCE_NAMES["pseudo"], "pseudopotentials")
DEFAULT_ORBITAL_DIRS = {
    "efficiency": resource_path(
        "ABACUSKIT_ORBITAL_EFFICIENCY_DIR",
        APNS_RESOURCE_NAMES["orbital_efficiency"],
        "orbitals/efficiency",
    ),
    "precision": resource_path(
        "ABACUSKIT_ORBITAL_PRECISION_DIR",
        APNS_RESOURCE_NAMES["orbital_precision"],
        "orbitals/precision",
    ),
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
ORBITAL_COLORS = {
    "s": "#1f77b4",
    "p": "#d62728",
    "d": "#2ca02c",
    "f": "#9467bd",
    "g": "#8c564b",
}
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
PRECISION_TARGET_KEYS = {"out_band", "out_proj_band", "out_dos"}
HYBRID_FUNCTIONALS = {"hse", "hse06", "pbe0", "hf"}
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


def choose_orbitals_for_symbols(symbols: Iterable[str], orbital_dir: Path) -> dict[str, Path]:
    orbital_index = index_library(orbital_dir, (".orb",))
    orbitals: dict[str, Path] = {}
    missing = []
    for sym in symbols:
        candidates = [p for p in orbital_index.get(sym, []) if p.name.startswith(f"{sym}_")]
        if not candidates:
            missing.append(f"orbital for {sym} in {orbital_dir}")
        else:
            orbitals[sym] = sorted(candidates, key=lambda p: score_orbital(p, sym))[0]
    if missing:
        die("missing library files: " + ", ".join(missing))
    return orbitals


def read_stru_species_symbols(stru: Path) -> list[str]:
    lines = stru.read_text(errors="ignore").splitlines()
    symbols: list[str] = []
    in_species = False
    for raw in lines:
        body = raw.split("#", 1)[0].strip()
        if not body:
            continue
        key = body.upper()
        if key == "ATOMIC_SPECIES":
            in_species = True
            continue
        if in_species and key in {
            "NUMERICAL_ORBITAL",
            "LATTICE_CONSTANT",
            "LATTICE_VECTORS",
            "ATOMIC_POSITIONS",
        }:
            break
        if in_species:
            symbols.append(body.split()[0])
    return symbols


def sync_stru_orbitals(stru: Path, orbital_dir: Path, backup: bool = True) -> bool:
    if not stru.is_file():
        return False
    lines = stru.read_text(errors="ignore").splitlines()
    symbols = read_stru_species_symbols(stru)
    if not symbols:
        return False
    orbitals = choose_orbitals_for_symbols(symbols, orbital_dir)
    start = None
    for idx, raw in enumerate(lines):
        if raw.split("#", 1)[0].strip().upper() == "NUMERICAL_ORBITAL":
            start = idx
            break
    new_orbital_lines = [orbitals[sym].name for sym in symbols]
    if start is None:
        insert = None
        for idx, raw in enumerate(lines):
            if raw.split("#", 1)[0].strip().upper() == "LATTICE_CONSTANT":
                insert = idx
                break
        if insert is None:
            die(f"cannot find insertion point for NUMERICAL_ORBITAL in {stru}")
        lines[insert:insert] = ["", "NUMERICAL_ORBITAL", *new_orbital_lines]
    else:
        line_idx = start + 1
        replaced = 0
        while line_idx < len(lines) and replaced < len(symbols):
            body = lines[line_idx].split("#", 1)[0].strip()
            if body:
                lines[line_idx] = new_orbital_lines[replaced]
                replaced += 1
            line_idx += 1
        if replaced < len(symbols):
            lines[start + 1 : line_idx] = new_orbital_lines
    new_text = "\n".join(lines) + "\n"
    old_text = stru.read_text(errors="ignore")
    if new_text == old_text:
        return False
    if backup:
        backup_path = stru.with_suffix(stru.suffix + ".orbital.bak") if stru.suffix else stru.with_name(stru.name + ".orbital.bak")
        backup_path.write_text(old_text)
    stru.write_text(new_text)
    return True


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


def clean_stru_lines(path: Path) -> list[str]:
    lines: list[str] = []
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            lines.append(line)
    return lines


def find_stru_section(lines: list[str], name: str) -> int | None:
    target = name.upper()
    for idx, line in enumerate(lines):
        if line.upper() == target:
            return idx
    return None


def parse_stru_atoms(path: Path):
    lines = clean_stru_lines(path)
    lattice_constant = 1.0
    idx = find_stru_section(lines, "LATTICE_CONSTANT")
    if idx is not None and idx + 1 < len(lines):
        lattice_constant = float(lines[idx + 1].split()[0])

    idx = find_stru_section(lines, "LATTICE_VECTORS")
    if idx is None or idx + 3 >= len(lines):
        die(f"cannot find LATTICE_VECTORS in {path}")
    cell_bohr = np.array([[float(x) for x in lines[idx + offset].split()[:3]] for offset in (1, 2, 3)])
    cell = cell_bohr * lattice_constant / BOHR_PER_ANGSTROM

    idx = find_stru_section(lines, "ATOMIC_POSITIONS")
    if idx is None or idx + 1 >= len(lines):
        die(f"cannot find ATOMIC_POSITIONS in {path}")
    coord_type = lines[idx + 1].strip().lower()
    scaled = coord_type.startswith(("direct", "crystal"))
    cartesian_angstrom = "angstrom" in coord_type

    symbols: list[str] = []
    positions: list[list[float]] = []
    cursor = idx + 2
    while cursor < len(lines):
        symbol = lines[cursor].split()[0]
        if symbol.upper() in {"ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS"}:
            break
        if symbol not in atomic_numbers:
            cursor += 1
            continue
        if cursor + 2 >= len(lines):
            die(f"incomplete ATOMIC_POSITIONS block for {symbol} in {path}")
        try:
            count = int(float(lines[cursor + 2].split()[0]))
        except ValueError:
            die(f"bad atom count for {symbol} in {path}: {lines[cursor + 2]}")
        start = cursor + 3
        for line in lines[start : start + count]:
            parts = line.split()
            if len(parts) < 3:
                die(f"bad atomic position line in {path}: {line}")
            symbols.append(symbol)
            positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
        cursor = start + count

    if not symbols:
        die(f"cannot parse atoms from {path}")
    if scaled:
        atoms = Atoms(symbols=symbols, scaled_positions=positions, cell=cell, pbc=True)
    else:
        factor = 1.0 if cartesian_angstrom else lattice_constant / BOHR_PER_ANGSTROM
        atoms = Atoms(symbols=symbols, positions=np.array(positions) * factor, cell=cell, pbc=True)
    return atoms


def canonicalize_axis_aligned_atoms(atoms):
    cell = np.array(atoms.cell.array, dtype=float)
    if cell.shape != (3, 3):
        return atoms

    tol = max(1.0e-8, float(np.max(np.abs(cell))) * 1.0e-8)
    axis_for_row: list[int] = []
    lengths = np.zeros(3)
    for row in cell:
        nonzero = [idx for idx, value in enumerate(row) if abs(value) > tol]
        if len(nonzero) != 1:
            return atoms
        axis = nonzero[0]
        if axis in axis_for_row:
            return atoms
        axis_for_row.append(axis)
        lengths[axis] = abs(float(row[axis]))

    if np.any(lengths <= tol):
        return atoms
    return Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=atoms.get_positions(),
        cell=np.diag(lengths),
        pbc=atoms.pbc,
    )


def fix_stru_range(
    stru: Path,
    out: Path,
    axis: str,
    lower: float,
    upper: float,
    backup: bool = True,
) -> int:
    if lower > upper:
        lower, upper = upper, lower
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    raw_lines = stru.read_text(errors="ignore").splitlines()
    original_lines = list(raw_lines)
    entries: list[tuple[int, str]] = []
    for raw_idx, raw in enumerate(raw_lines):
        line = raw.split("#", 1)[0].strip()
        if line:
            entries.append((raw_idx, line))
    cleaned = [line for _, line in entries]
    idx = find_stru_section(cleaned, "ATOMIC_POSITIONS")
    if idx is None or idx + 1 >= len(cleaned):
        die(f"cannot find ATOMIC_POSITIONS in {stru}")
    coord_type = cleaned[idx + 1].strip().lower()
    if not coord_type.startswith("cartesian"):
        die("fix-stru-range currently supports Cartesian ATOMIC_POSITIONS only")

    fixed = 0
    cursor = idx + 2
    while cursor < len(cleaned):
        symbol = cleaned[cursor].split()[0]
        if symbol.upper() in {"ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS"}:
            break
        if symbol not in atomic_numbers:
            cursor += 1
            continue
        if cursor + 2 >= len(cleaned):
            die(f"incomplete ATOMIC_POSITIONS block for {symbol} in {stru}")
        count = int(float(cleaned[cursor + 2].split()[0]))
        start = cursor + 3
        for line_idx in range(start, start + count):
            if line_idx >= len(cleaned):
                die(f"incomplete coordinate list for {symbol} in {stru}")
            parts = cleaned[line_idx].split()
            if len(parts) < 3:
                die(f"bad atomic position line in {stru}: {cleaned[line_idx]}")
            coord = float(parts[axis_index])
            if lower <= coord <= upper:
                coords = [float(parts[0]), float(parts[1]), float(parts[2])]
                raw_idx = entries[line_idx][0]
                raw_lines[raw_idx] = f"{coords[0]:18.10f} {coords[1]:18.10f} {coords[2]:18.10f} 0 0 0"
                fixed += 1
        cursor = start + count

    if fixed == 0:
        print(f"No atoms matched {axis} in [{lower:g}, {upper:g}].")
    if out.resolve() == stru.resolve() and backup:
        backup_path = stru.with_suffix(stru.suffix + ".bak") if stru.suffix else stru.with_name(stru.name + ".bak")
        backup_path.write_text("\n".join(original_lines) + "\n")
        print(f"backed up original STRU to {backup_path}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(raw_lines) + "\n")
    return fixed


def swap_yz_values(values: list[float]) -> list[float]:
    return [values[0], values[2], values[1]]


def rotate_stru_vacuum_z_to_y(
    stru: Path,
    out: Path,
    backup: bool = True,
) -> None:
    raw_lines = stru.read_text(errors="ignore").splitlines()
    original_lines = list(raw_lines)
    entries: list[tuple[int, str]] = []
    for raw_idx, raw in enumerate(raw_lines):
        line = raw.split("#", 1)[0].strip()
        if line:
            entries.append((raw_idx, line))
    cleaned = [line for _, line in entries]

    idx = find_stru_section(cleaned, "LATTICE_VECTORS")
    if idx is None or idx + 3 >= len(cleaned):
        die(f"cannot find LATTICE_VECTORS in {stru}")
    for line_idx in range(idx + 1, idx + 4):
        parts = cleaned[line_idx].split()
        if len(parts) < 3:
            die(f"bad lattice vector line in {stru}: {cleaned[line_idx]}")
        vector = swap_yz_values([float(parts[0]), float(parts[1]), float(parts[2])])
        raw_idx = entries[line_idx][0]
        raw_lines[raw_idx] = f"{vector[0]:18.10f} {vector[1]:18.10f} {vector[2]:18.10f}"

    idx = find_stru_section(cleaned, "ATOMIC_POSITIONS")
    if idx is None or idx + 1 >= len(cleaned):
        die(f"cannot find ATOMIC_POSITIONS in {stru}")
    coord_type = cleaned[idx + 1].strip().lower()
    if not coord_type.startswith("cartesian"):
        die("rotate-vacuum-z-to-y currently supports Cartesian ATOMIC_POSITIONS only")

    cursor = idx + 2
    while cursor < len(cleaned):
        symbol = cleaned[cursor].split()[0]
        if symbol.upper() in {"ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS"}:
            break
        if symbol not in atomic_numbers:
            cursor += 1
            continue
        if cursor + 2 >= len(cleaned):
            die(f"incomplete ATOMIC_POSITIONS block for {symbol} in {stru}")
        count = int(float(cleaned[cursor + 2].split()[0]))
        start = cursor + 3
        for line_idx in range(start, start + count):
            if line_idx >= len(cleaned):
                die(f"incomplete coordinate list for {symbol} in {stru}")
            parts = cleaned[line_idx].split()
            if len(parts) < 3:
                die(f"bad atomic position line in {stru}: {cleaned[line_idx]}")
            coords = swap_yz_values([float(parts[0]), float(parts[1]), float(parts[2])])
            flags = parts[3:6] if len(parts) >= 6 else ["1", "1", "1"]
            flags = [flags[0], flags[2], flags[1]]
            raw_idx = entries[line_idx][0]
            raw_lines[raw_idx] = (
                f"{coords[0]:18.10f} {coords[1]:18.10f} {coords[2]:18.10f} "
                f"{flags[0]} {flags[1]} {flags[2]}"
            )
        cursor = start + count

    if out.resolve() == stru.resolve() and backup:
        backup_path = stru.with_suffix(stru.suffix + ".bak") if stru.suffix else stru.with_name(stru.name + ".bak")
        backup_path.write_text("\n".join(original_lines) + "\n")
        print(f"backed up original STRU to {backup_path}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(raw_lines) + "\n")


def read_atoms(path: Path, supercell: tuple[int, int, int] = (1, 1, 1)):
    if path.name.upper() == "STRU":
        atoms = parse_stru_atoms(path)
    else:
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
    def add_default(key: str, value: object) -> None:
        if key not in extra:
            params.append((key, value))

    params: list[tuple[str, object]] = [
        ("suffix", suffix),
        ("calculation", calculation),
        ("stru_file", "STRU"),
        ("pseudo_dir", pseudo_dir),
    ]
    if basis_type == "lcao" and "orbital_dir" not in extra:
        params.append(("orbital_dir", orbital_dir))
    add_default("basis_type", basis_type)
    add_default("ks_solver", ks_solver)
    add_default("device", device)
    add_default("symmetry", 0)
    add_default("gamma_only", 0)
    add_default("kspacing", kspacing)
    if "dft_functional" not in extra:
        params.append(("dft_functional", "PBE"))
    add_default("ecutwfc", ecutwfc)
    add_default("nspin", nspin)
    add_default("scf_thr", "1e-6")
    add_default("scf_nmax", 300)
    add_default("mixing_type", "broyden")
    add_default("mixing_beta", 0.10)
    add_default("mixing_ndim", 20)
    add_default("cal_force", 1)
    add_default("cal_stress", 1 if cal_stress else 0)
    add_default("out_wfc_lcao", 0)
    add_default("out_chg", 0)
    add_default("smearing_method", "gauss")
    add_default("smearing_sigma", 0.015)
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


def precision_target_requested(args, extra: dict[str, str]) -> bool:
    if getattr(args, "dos", False):
        return True
    dft_functional = extra.get("dft_functional", "").strip().lower()
    if dft_functional in HYBRID_FUNCTIONALS or any(key.startswith("exx_") for key in extra):
        return True
    if (getattr(args, "kind", "") or "").lower() != "nscf":
        return False
    return any(key in extra for key in PRECISION_TARGET_KEYS)


def input_template_orbital_dir(args, extra: dict[str, str]) -> tuple[Path, bool]:
    if "orbital_dir" in extra:
        return Path(extra["orbital_dir"]), False
    if args.orbital_dir:
        return args.orbital_dir, False
    if precision_target_requested(args, extra):
        return DEFAULT_ORBITAL_DIRS["precision"], True
    return DEFAULT_ORBITAL_DIRS[args.orbital_quality], False


def input_template_params(args) -> list[tuple[str, object]]:
    extra = parse_key_values(args.set)
    orbital_dir, _ = input_template_orbital_dir(args, extra)
    def add_default(key: str, value: object) -> None:
        if key not in extra:
            params.append((key, value))

    params: list[tuple[str, object]] = [
        ("suffix", args.suffix),
        ("calculation", args.kind),
        ("stru_file", "STRU"),
        ("pseudo_dir", args.pseudo_dir),
    ]
    if args.basis_type == "lcao" and "orbital_dir" not in extra:
        params.append(("orbital_dir", orbital_dir))
    add_default("basis_type", args.basis_type)
    add_default("ks_solver", args.ks_solver)
    add_default("device", args.device)
    add_default("symmetry", 0)
    add_default("gamma_only", 0)
    if not getattr(args, "no_kspacing", False):
        add_default("kspacing", args.kspacing)
    if "dft_functional" not in extra:
        params.append(("dft_functional", "PBE"))
    add_default("ecutwfc", args.ecutwfc)
    add_default("nspin", args.nspin)
    add_default("scf_thr", "1e-6")
    add_default("scf_nmax", 300)
    add_default("mixing_type", "broyden")
    add_default("mixing_beta", 0.10)
    add_default("mixing_ndim", 20)
    add_default("cal_force", 1)
    add_default("cal_stress", 1 if args.kind == "relax" or args.cal_stress else 0)
    add_default("out_wfc_lcao", 0)
    add_default("out_chg", 0)
    add_default("smearing_method", "gauss")
    add_default("smearing_sigma", 0.015)
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
    if getattr(args, "hybrid_hse_scf", False) and getattr(args, "hybrid_hse_band", False):
        die("--hybrid-hse-scf and --hybrid-hse-band are mutually exclusive")
    if getattr(args, "hybrid_hse_scf", False) or getattr(args, "hybrid_hse_band", False):
        user_settings = list(args.set or [])
        template_settings: list[str]
        args.basis_type = "lcao"
        args.device = "cpu"
        args.ks_solver = "genelpa"
        args.orbital_quality = "precision"
        args.orbital_dir = None
        if args.hybrid_hse_scf:
            args.kind = "scf"
            args.no_kspacing = False
            template_settings = ["dft_functional=hse", "exx_hybrid_step=100", "out_chg=1"]
        else:
            args.kind = "nscf"
            args.no_kspacing = True
            template_settings = [
                "dft_functional=hse",
                "exx_hybrid_step=100",
                "kpoint_file=KPT",
                "init_chg=file",
                "read_file_dir=./",
                "out_band=1",
                "cal_force=0",
                "cal_stress=0",
                "smearing_method=gaussian",
                "smearing_sigma=0.02",
            ]
        args.set = template_settings + user_settings
    extra = parse_key_values(args.set)
    orbital_dir, auto_precision = input_template_orbital_dir(args, extra)
    text = format_input_params(input_template_params(args), with_comments=not args.no_comments)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"wrote {args.kind} INPUT template to {args.out}")
    stru = args.out.parent / "STRU"
    if auto_precision and getattr(args, "basis_type", "lcao") == "lcao" and stru.is_file():
        changed = sync_stru_orbitals(stru, orbital_dir)
        status = "updated" if changed else "already used"
        print(f"{status} precision numerical orbitals in {stru}")


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


def seekpath_label(label: str) -> str:
    clean = label.replace("\\Gamma", "GAMMA").replace("Gamma", "GAMMA").replace("Γ", "GAMMA")
    return clean.replace("$", "").replace("\\", "")


def atoms_to_seekpath_structure(atoms):
    if atoms.cell is None or abs(atoms.cell.volume) < 1.0e-8:
        die("structure has no valid periodic cell")
    numbers = atoms.get_atomic_numbers()
    if not len(numbers):
        die("structure has no atoms")
    return (
        atoms.cell.array,
        atoms.get_scaled_positions(wrap=True),
        numbers,
    )


def get_seekpath_result(args) -> dict:
    try:
        import seekpath
    except ImportError as exc:
        die("high-symmetry KPT generation needs seekpath and spglib; reinstall with: pip install -U abacuskit")
        raise exc

    atoms = read_atoms(args.structure)
    structure = atoms_to_seekpath_structure(atoms)
    return seekpath.get_path(
        structure,
        with_time_reversal=not args.no_time_reversal,
        recipe="hpkot",
        threshold=args.threshold,
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
    )


def kpath_special_points(seek_result: dict, points_per_segment: int) -> list[tuple[str, tuple[float, float, float], int]]:
    coords = seek_result["point_coords"]
    path = seek_result["path"]
    special: list[tuple[str, tuple[float, float, float], int]] = []
    previous_end: str | None = None
    for start, end in path:
        if previous_end != start:
            special.append((start, tuple(coords[start]), points_per_segment))
        elif special:
            special[-1] = (special[-1][0], special[-1][1], points_per_segment)
        special.append((end, tuple(coords[end]), 1))
        previous_end = end
    if not special:
        die("seekpath did not return a high-symmetry k-path")
    return special


def write_kpt_path(
    path: Path,
    special_points: list[tuple[str, tuple[float, float, float], int]],
    comment: str,
) -> None:
    lines = [
        "K_POINTS",
        str(len(special_points)),
        "Line",
    ]
    for label, coords, count in special_points:
        lines.append(
            f"{coords[0]: .10f} {coords[1]: .10f} {coords[2]: .10f} {count:d} # {seekpath_label(label)}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    if comment:
        print(comment)


def write_high_symmetry_points(path: Path, seek_result: dict) -> None:
    lines = [
        "# High-symmetry points generated by SeeK-path/HPKOT",
        f"# bravais_lattice: {seek_result.get('bravais_lattice', 'unknown')}",
        f"# bravais_lattice_extended: {seek_result.get('bravais_lattice_extended', 'unknown')}",
    ]
    for label, coords in sorted(seek_result["point_coords"].items()):
        lines.append(f"{seekpath_label(label):<12} {coords[0]: .10f} {coords[1]: .10f} {coords[2]: .10f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def cmd_kpt(args) -> None:
    write_kpt(args.out, tuple(args.mesh), tuple(args.shift), args.model)
    print(f"wrote {args.model} KPT to {args.out}")


def reciprocal_mesh_from_kspacing(atoms, kspacing: float) -> tuple[int, int, int]:
    if kspacing <= 0:
        die("kspacing must be positive")
    reciprocal_lengths = atoms.cell.reciprocal().lengths() * (2.0 * math.pi)
    return tuple(max(1, int(math.ceil(float(length) / kspacing))) for length in reciprocal_lengths)


def cmd_kpt_path(args) -> None:
    if args.points_per_segment < 1:
        die("--points-per-segment must be at least 1")
    seek_result = get_seekpath_result(args)
    special_points = kpath_special_points(seek_result, args.points_per_segment)
    write_kpt_path(
        args.out,
        special_points,
        (
            "SeeK-path bravais lattice: "
            f"{seek_result.get('bravais_lattice', 'unknown')} "
            f"({seek_result.get('bravais_lattice_extended', 'unknown')})"
        ),
    )
    if args.high_symmetry_points:
        write_high_symmetry_points(args.high_symmetry_points, seek_result)
        print(f"wrote high-symmetry point table {args.high_symmetry_points}")
    print(f"wrote SeeK-path line-mode KPT to {args.out}")


def find_default_structure_file(root: Path = Path(".")) -> Path | None:
    candidates = [
        root / "STRU",
        root / "POSCAR",
        root / "CONTCAR",
    ]
    candidates.extend(sorted(root.glob("*.cif"), key=lambda p: natural_key(p.name)))
    for path in candidates:
        if path.is_file():
            return path
    return None


def cmd_kpt_auto(args) -> None:
    structure = args.structure or find_default_structure_file()
    if not structure:
        die("cannot auto-generate KPT because no STRU/POSCAR/CONTCAR/CIF was found in current directory")

    if args.mode in {"auto", "path"}:
        path_args = argparse.Namespace(
            structure=structure,
            out=args.out,
            high_symmetry_points=args.high_symmetry_points,
            points_per_segment=args.points_per_segment,
            symprec=args.symprec,
            angle_tolerance=args.angle_tolerance,
            threshold=args.threshold,
            no_time_reversal=args.no_time_reversal,
        )
        try:
            cmd_kpt_path(path_args)
            return
        except SystemExit:
            if args.mode == "path":
                raise
            print("High-symmetry KPT generation failed; falling back to regular Gamma mesh.")

    atoms = read_atoms(structure)
    mesh = reciprocal_mesh_from_kspacing(atoms, args.kspacing)
    write_kpt(args.out, mesh, (0, 0, 0), "gamma")
    print(f"wrote Gamma KPT to {args.out}; kspacing={args.kspacing:g}, mesh={mesh[0]} {mesh[1]} {mesh[2]}")


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
        r"scf\s+is\s+converged",
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


def split_band_table_blocks(data: np.ndarray) -> list[np.ndarray]:
    if data.shape[0] < 2 or data.shape[1] < 2:
        return [data]
    resets = np.where((np.diff(data[:, 0]) < 0) | (np.diff(data[:, 1]) < -1.0e-10))[0] + 1
    if resets.size == 0:
        return [data]
    starts = [0, *[int(i) for i in resets], data.shape[0]]
    return [data[starts[i] : starts[i + 1]] for i in range(len(starts) - 1) if starts[i + 1] > starts[i]]


def select_band_table_block(data: np.ndarray) -> tuple[np.ndarray, str | None]:
    blocks = split_band_table_blocks(data)
    if len(blocks) == 1:
        return data, None
    first = blocks[0]
    same_blocks = [
        block
        for block in blocks[1:]
        if block.shape == first.shape and np.allclose(block[:, 1:], first[:, 1:], rtol=0.0, atol=1.0e-8)
    ]
    if len(same_blocks) == len(blocks) - 1:
        return first, f"ignored {len(blocks) - 1} duplicated band block(s)"
    return first, f"found {len(blocks)} band blocks; plotted the first block"


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
        chunks = re.split(r"[;\s]+", item)
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
            orbital_text = orbitals.replace(" ", "")
            if re.search(r"[,/+]", orbital_text):
                orbital_items = [orb for orb in re.split(r"[,/+]+", orbital_text) if orb]
            else:
                orbital_items = list(orbital_text)
            for orb in orbital_items:
                orb = orb.lower()
                if orb not in ORBITAL_LABEL_TO_L:
                    die(f"bad orbital label {orb!r}; use s, p, d, f, or g")
                selectors.add((sym.strip(), orb))
    return selectors


def maybe_shift_fermi(x: np.ndarray, fermi: float | None) -> tuple[np.ndarray, str]:
    if fermi is None:
        return x, "Energy (eV)"
    return x - fermi, "Energy - E_F (eV)"


def parse_fermi_energy(text: str) -> float | None:
    patterns = [
        r"\bEFERMI\s*=\s*([-+0-9.eE]+)\s*eV",
        r"\bE_Fermi(?:_up|_dw)?\s+[-+0-9.eE]+\s+([-+0-9.eE]+)",
        r"\bE_Fermi(?:_up|_dw)?\s*[:=]\s*([-+0-9.eE]+)\s*eV",
    ]
    for pattern in patterns:
        value = parse_last_float(pattern, text)
        if value is not None:
            return value
    return None


def find_fermi_energy(root: Path) -> tuple[float | None, Path | None]:
    primary: list[Path] = []
    fallback: list[Path] = []
    search_roots = [root.parent if root.is_file() else root]
    outdir = find_abacus_outdir(search_roots[0])
    if outdir:
        search_roots.append(outdir)
    for base in search_roots:
        if base.is_file() and base.name.startswith("running_"):
            primary.append(base)
        elif base.is_dir():
            primary.extend(sorted(base.glob("running_*.log"), key=lambda p: (p.stat().st_mtime, natural_key(p.name))))
            outdir = find_abacus_outdir(base)
            if outdir and outdir != base:
                primary.extend(sorted(outdir.glob("running_*.log"), key=lambda p: (p.stat().st_mtime, natural_key(p.name))))
            scf_out = base.parent / "scf" / "OUT.ABACUS"
            if scf_out.is_dir():
                fallback.extend(sorted(scf_out.glob("running_*.log"), key=lambda p: (p.stat().st_mtime, natural_key(p.name))))
    seen: set[Path] = set()
    for group in (primary, fallback):
        unique = []
        for path in group:
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                unique.append(path)
        for log in reversed(unique):
            fermi = parse_fermi_energy(log.read_text(errors="ignore"))
            if fermi is not None:
                return fermi, log
    return None, None


def find_band_file(root: Path, explicit: Path | None = None) -> Path | None:
    if explicit:
        return explicit.expanduser()
    if root.is_file():
        return root
    outdir = resolve_out_path(root)
    return find_first_file(outdir, ["BANDS_1.dat", "BANDS_1", "band.txt"], ["BANDS*.dat", "BANDS*", "band*.txt"])


def kpoint_label(coords: tuple[float, float, float]) -> str:
    common = {
        (0.0, 0.0, 0.0): r"$\Gamma$",
        (0.5, 0.0, 0.0): "X",
        (0.0, 0.5, 0.0): "Y",
        (0.0, 0.0, 0.5): "Z",
        (0.5, 0.5, 0.0): "M",
        (0.5, 0.0, 0.5): "A",
        (0.0, 0.5, 0.5): "B",
        (0.5, 0.5, 0.5): "R",
    }
    rounded = tuple(round(x, 6) for x in coords)
    if rounded in common:
        return common[rounded]
    return "(" + ",".join(f"{x:g}" for x in rounded) + ")"


def normalize_high_symmetry_label(label: str) -> str:
    clean = label.strip()
    if not clean:
        return clean
    upper = clean.upper()
    if upper in {"GAMMA", "G", "Γ"}:
        return r"$\Gamma$"
    return clean.replace("Γ", r"$\Gamma$")


def load_high_symmetry_points(path: Path | None) -> dict[str, tuple[float, float, float]]:
    if not path or not path.is_file():
        return {}
    points: dict[str, tuple[float, float, float]] = {}
    for raw in path.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            coords = tuple(float(v) for v in parts[1:4])
        except ValueError:
            continue
        points[parts[0].strip()] = coords
    return points


def parse_kpoint_label(comment: str, coords: tuple[float, float, float], high_symmetry_points: Path | None = None) -> str:
    label = comment.strip()
    if label:
        return normalize_high_symmetry_label(label.split()[0])
    points = load_high_symmetry_points(high_symmetry_points)
    rounded = tuple(round(x, 10) for x in coords)
    for name, ref in points.items():
        if tuple(round(x, 10) for x in ref) == rounded:
            return normalize_high_symmetry_label(name)
    return normalize_high_symmetry_label(kpoint_label(coords))


def find_kpt_file(root: Path) -> Path | None:
    if root.is_file():
        root = root.parent
    candidates = [root / "KPT", root.parent / "KPT"]
    outdir = find_abacus_outdir(root)
    if outdir:
        candidates += [outdir / "KPT", outdir.parent / "KPT"]
    for path in candidates:
        if path.is_file():
            return path
    return None


def find_high_symmetry_points_file(root: Path, kpt: Path | None = None) -> Path | None:
    bases: list[Path] = []
    if kpt and kpt.is_file():
        bases.extend([kpt.parent, kpt.parent.parent])
    if root.is_file():
        bases.extend([root.parent, root.parent.parent])
    else:
        bases.extend([root, root.parent])
        outdir = find_abacus_outdir(root)
        if outdir:
            bases.extend([outdir, outdir.parent, outdir.parent.parent])
    seen: set[Path] = set()
    for base in bases:
        candidate = base / "HIGH_SYMMETRY_POINTS"
        try:
            resolved = candidate.resolve()
        except OSError:
            resolved = candidate.absolute()
        if resolved in seen:
            continue
        seen.add(resolved)
        if candidate.is_file():
            return candidate
    return None


def parse_line_kpt_ticks(
    kpt: Path | None,
    x: np.ndarray,
    high_symmetry_points: Path | None = None,
) -> tuple[list[float], list[str]]:
    if not kpt or not kpt.is_file():
        return [float(x[0]), float(x[-1])], ["", ""]
    raw_lines = [line.strip() for line in kpt.read_text(errors="ignore").splitlines()]
    entries = []
    for raw in raw_lines:
        body, _, comment = raw.partition("#")
        body = body.strip()
        if body:
            entries.append((body, comment.strip()))
    if len(entries) < 4 or entries[2][0].lower() != "line":
        return [float(x[0]), float(x[-1])], ["", ""]
    try:
        npoint = int(float(entries[1][0].split()[0]))
    except ValueError:
        return [float(x[0]), float(x[-1])], ["", ""]
    points = []
    for body, comment in entries[3 : 3 + npoint]:
        parts = body.split()
        if len(parts) < 4:
            continue
        try:
            coords = tuple(float(v) for v in parts[:3])
            count = int(float(parts[3]))
        except ValueError:
            continue
        label = parse_kpoint_label(comment, coords, high_symmetry_points)
        points.append((coords, count, label))
    if len(points) < 2:
        return [float(x[0]), float(x[-1])], ["", ""]
    indices = [0]
    total = 0
    for _, count, _ in points[:-1]:
        total += max(count, 1)
        indices.append(min(total, len(x) - 1))
    ticks = [float(x[i]) for i in indices]
    labels = [label for _, _, label in points[: len(ticks)]]
    if len(ticks) > 2:
        merged_ticks: list[float] = []
        merged_labels: list[str] = []
        min_sep = max(float(x[-1] - x[0]) * 0.025, 1.0e-8)
        for tick, label in zip(ticks, labels):
            if merged_ticks and abs(tick - merged_ticks[-1]) < min_sep:
                merged_ticks[-1] = 0.5 * (merged_ticks[-1] + tick)
                if label and label not in merged_labels[-1].split("|"):
                    merged_labels[-1] = f"{merged_labels[-1]}|{label}" if merged_labels[-1] else label
            else:
                merged_ticks.append(tick)
                merged_labels.append(label)
        ticks, labels = merged_ticks, merged_labels
    return ticks, labels


def line_kpt_break_indices(kpt: Path | None, nk: int) -> list[int]:
    if not kpt or not kpt.is_file():
        return []
    raw_lines = [line.strip() for line in kpt.read_text(errors="ignore").splitlines()]
    entries = []
    for raw in raw_lines:
        body = raw.partition("#")[0].strip()
        if body:
            entries.append(body)
    if len(entries) < 4 or entries[2].lower() != "line":
        return []
    try:
        npoint = int(float(entries[1].split()[0]))
    except ValueError:
        return []
    counts: list[int] = []
    for line in entries[3 : 3 + npoint]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            counts.append(int(float(parts[3])))
        except ValueError:
            continue
    breaks: list[int] = []
    total = 0
    for count in counts[:-1]:
        total += max(count, 1)
        if count <= 1 and 0 < total < nk:
            breaks.append(total)
    return breaks


def compact_line_kpt_breaks(x: np.ndarray, break_indices: list[int]) -> np.ndarray:
    if not break_indices:
        return x
    compact = np.array(x, dtype=float, copy=True)
    for start in sorted(set(break_indices)):
        if start <= 0 or start >= len(compact):
            continue
        jump = float(compact[start] - compact[start - 1])
        if jump > 1.0e-10:
            compact[start:] -= jump
    return compact


def find_pband_file(root: Path) -> Path | None:
    if root.is_file():
        root = root.parent
    candidates = [
        root / "pbands1.xml",
        root / "PBANDS1.xml",
        root / "pbands.xml",
    ]
    outdir = find_abacus_outdir(root)
    if outdir:
        candidates.extend(
            [
                outdir / "pbands1.xml",
                outdir / "PBANDS1.xml",
                outdir / "pbands.xml",
            ]
        )
    for path in candidates:
        if path.is_file():
            return path
    for pattern in ("pbands*.xml", "PBANDS*.xml"):
        found = sorted([p for p in root.glob(pattern) if p.is_file()], key=lambda p: natural_key(p.name))
        if found:
            return found[0]
        if outdir:
            found = sorted([p for p in outdir.glob(pattern) if p.is_file()], key=lambda p: natural_key(p.name))
            if found:
                return found[0]
    return None


def find_input_info_file(root: Path) -> Path | None:
    if root.is_file():
        root = root.parent
    candidates = [root / "INPUT.info"]
    outdir = find_abacus_outdir(root)
    if outdir:
        candidates.append(outdir / "INPUT.info")
    for path in candidates:
        if path.is_file():
            return path
    return None


def parse_projected_band_file(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], list[dict[str, object]]]:
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        die(f"cannot parse projected band XML {path}: {exc}")
    band_node = root.find("band_structure")
    if band_node is None:
        die(f"cannot find band_structure in {path}")
    band_text = band_node.text or ""
    band_rows = [line.strip() for line in band_text.splitlines() if line.strip()]
    bands = np.array([[float(x) for x in row.split()] for row in band_rows], dtype=float)
    x = np.arange(bands.shape[0], dtype=float)
    channels: dict[str, np.ndarray] = {label: np.zeros_like(bands) for label in ORBITAL_COLORS}
    orbitals: list[dict[str, object]] = []
    for orbital in root.findall("orbital"):
        l_raw = orbital.attrib.get("l", "0").strip()
        try:
            l_value = int(float(l_raw))
        except ValueError:
            continue
        label = L_TO_ORBITAL_LABEL.get(l_value)
        if label is None:
            continue
        data_node = orbital.find("data")
        if data_node is None or data_node.text is None:
            continue
        data_rows = [line.strip() for line in data_node.text.splitlines() if line.strip()]
        if not data_rows:
            continue
        matrix = np.array([[float(x) for x in row.split()] for row in data_rows], dtype=float)
        if matrix.shape[0] != bands.shape[0]:
            continue
        if matrix.shape[1] < bands.shape[1]:
            continue
        channels[label][:, :] += np.clip(matrix[:, : bands.shape[1]], 0.0, None)
        orbitals.append(
            {
                "label": label,
                "species": orbital.attrib.get("species", "").strip(),
                "atom_index": orbital.attrib.get("atom_index", "").strip(),
            }
        )
    return x, bands, channels, orbitals


def project_band_weights(pband_file: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, object]]]:
    try:
        root = ET.parse(pband_file).getroot()
    except ET.ParseError as exc:
        die(f"cannot parse projected band XML {pband_file}: {exc}")
    band_node = root.find("band_structure")
    if band_node is None:
        die(f"cannot find band_structure in {pband_file}")
    band_text = band_node.text or ""
    band_rows = [line.strip() for line in band_text.splitlines() if line.strip()]
    bands = np.array([[float(x) for x in row.split()] for row in band_rows], dtype=float)
    x = np.arange(bands.shape[0], dtype=float)

    orbital_rows: list[np.ndarray] = []
    orbitals: list[dict[str, object]] = []
    for orbital in root.findall("orbital"):
        data_node = orbital.find("data")
        if data_node is None or data_node.text is None:
            continue
        data_rows = [line.strip() for line in data_node.text.splitlines() if line.strip()]
        if not data_rows:
            continue
        matrix = np.array([[float(x) for x in row.split()] for row in data_rows], dtype=float)
        if matrix.shape != bands.shape:
            continue
        orbital_rows.append(np.clip(matrix, 0.0, None))
        try:
            l_value = int(float(orbital.attrib.get("l", "0").strip()))
        except ValueError:
            l_value = -1
        orbitals.append(
            {
                "label": L_TO_ORBITAL_LABEL.get(l_value, ""),
                "species": orbital.attrib.get("species", "").strip(),
                "atom_index": orbital.attrib.get("atom_index", "").strip(),
                "l": l_value,
            }
        )
    if not orbital_rows:
        die(f"no projected orbital data parsed from {pband_file}")
    weights = np.stack(orbital_rows, axis=0)
    total = weights.sum(axis=0)
    total[total <= 0.0] = 1.0
    return x, bands, weights / total, orbitals


def format_band_label(label: str) -> str:
    return normalize_high_symmetry_label(label)


def nearest_k_label(k_value: float, ticks: list[float], labels: list[str]) -> str:
    if not ticks or not labels:
        return f"k={k_value:.6f}"
    idx = min(range(len(ticks)), key=lambda i: abs(ticks[i] - k_value))
    if abs(ticks[idx] - k_value) < 1.0e-8 and idx < len(labels) and labels[idx]:
        return labels[idx]
    return f"k={k_value:.6f}"


def analyze_band_gap(
    shifted_bands: np.ndarray,
    x: np.ndarray,
    ticks: list[float],
    labels: list[str],
    tol: float = 1.0e-5,
) -> dict[str, object]:
    occupied = shifted_bands <= tol
    unoccupied = shifted_bands > tol
    if not occupied.any() or not unoccupied.any():
        return {"gap_ev": 0.0, "kind": "metallic", "message": "metallic or incomplete occupied/unoccupied bands"}

    vbm = float(shifted_bands[occupied].max())
    cbm = float(shifted_bands[unoccupied].min())
    gap = max(0.0, cbm - vbm)
    v_indices = np.argwhere(np.isclose(shifted_bands, vbm, atol=tol))
    c_indices = np.argwhere(np.isclose(shifted_bands, cbm, atol=tol))
    v_k = int(v_indices[0][0])
    c_k = int(c_indices[0][0])
    direct = any(int(v[0]) == int(c[0]) for v in v_indices for c in c_indices)
    kind = "direct" if direct else "indirect"
    return {
        "gap_ev": gap,
        "kind": kind,
        "vbm_ev": vbm,
        "cbm_ev": cbm,
        "vbm_k": float(x[v_k]),
        "cbm_k": float(x[c_k]),
        "vbm_label": nearest_k_label(float(x[v_k]), ticks, labels),
        "cbm_label": nearest_k_label(float(x[c_k]), ticks, labels),
    }


def format_band_gap_message(gap_info: dict[str, object]) -> str:
    if gap_info.get("kind") == "metallic":
        return "Band gap: 0.000000 eV (metallic)"
    return (
        f"Band gap: {gap_info['gap_ev']:.6f} eV ({gap_info['kind']}); "
        f"VBM {gap_info['vbm_ev']:.6f} eV at {gap_info['vbm_label']}, "
        f"CBM {gap_info['cbm_ev']:.6f} eV at {gap_info['cbm_label']}"
    )


def contiguous_k_ranges(x: np.ndarray) -> list[slice]:
    if len(x) < 2:
        return [slice(0, len(x))]
    breaks = np.where(np.diff(x) <= 1.0e-10)[0] + 1
    starts = [0, *[int(i) for i in breaks], len(x)]
    return [slice(starts[i], starts[i + 1]) for i in range(len(starts) - 1) if starts[i + 1] > starts[i]]


def grouped_orbital_weights(weights: np.ndarray, orbitals: list[dict[str, object]]) -> dict[str, np.ndarray]:
    groups: dict[str, np.ndarray] = {}
    for idx, orbital in enumerate(orbitals):
        label = str(orbital.get("label", ""))
        if label not in ORBITAL_COLORS:
            continue
        if label not in groups:
            groups[label] = np.zeros_like(weights[idx])
        groups[label] += weights[idx]
    return groups


def load_band_plot_data(
    root: Path,
    explicit_file: Path | None = None,
    explicit_kpt: Path | None = None,
    fermi: float | None = None,
) -> dict[str, object]:
    band_file = find_band_file(root, explicit_file)
    if not band_file or not band_file.is_file():
        die(f"cannot find BANDS_*.dat or band.txt under {root}")
    data, block_note = select_band_table_block(read_numeric_table(band_file))
    if data.shape[1] < 3:
        die(f"band file needs at least 3 columns: {band_file}")
    fermi_log = None
    if fermi is None:
        fermi, fermi_log = find_fermi_energy(root if not root.is_file() else band_file.parent)
    if fermi is None:
        die("cannot determine Fermi energy from running_*.log; pass --fermi explicitly")
    x = data[:, 1]
    bands = data[:, 2:] - fermi
    kpt_file = explicit_kpt or find_kpt_file(root)
    high_symmetry_points = find_high_symmetry_points_file(root, kpt_file)
    x = compact_line_kpt_breaks(x, line_kpt_break_indices(kpt_file, len(x)))
    ticks, labels = parse_line_kpt_ticks(kpt_file, x, high_symmetry_points)
    return {
        "band_file": band_file,
        "block_note": block_note,
        "fermi": fermi,
        "fermi_log": fermi_log,
        "x": x,
        "bands": bands,
        "ticks": ticks,
        "labels": labels,
    }


def load_projected_band_channels(
    pband_file: Path | None,
    bands_shape: tuple[int, int],
) -> dict[str, np.ndarray]:
    if not pband_file:
        return {}
    _, pband_bands, pband_weights, pband_orbitals = project_band_weights(pband_file)
    if pband_bands.shape != bands_shape:
        print(f"note: projected band shape {pband_bands.shape} does not match band shape {bands_shape}; band plot is monochrome")
        return {}
    return grouped_orbital_weights(pband_weights, pband_orbitals)


def plot_band_axes(
    ax,
    x: np.ndarray,
    bands: np.ndarray,
    ticks: list[float],
    labels: list[str],
    pband_channels: dict[str, np.ndarray] | None = None,
    linewidth: float = 0.8,
    color: str = "C0",
) -> None:
    pband_channels = pband_channels or {}
    ranges = contiguous_k_ranges(x)
    background_color = "0.18" if pband_channels else color
    background_alpha = 0.35 if pband_channels else 1.0
    for i in range(bands.shape[1]):
        for band_range in ranges:
            if band_range.stop - band_range.start < 2:
                continue
            ax.plot(
                x[band_range],
                bands[band_range, i],
                lw=linewidth,
                color=background_color,
                alpha=background_alpha,
                zorder=1,
            )
    if pband_channels:
        x_points = np.repeat(x, bands.shape[1])
        y_points = bands.reshape(-1)
        for orbital in ("s", "p", "d", "f", "g"):
            values = pband_channels.get(orbital)
            if values is None:
                continue
            flat = np.clip(values.reshape(-1), 0.0, 1.0)
            mask = flat > 0.025
            if not mask.any():
                continue
            sizes = 2.0 + 20.0 * flat[mask]
            ax.scatter(
                x_points[mask],
                y_points[mask],
                s=sizes,
                c=ORBITAL_COLORS[orbital],
                alpha=0.42,
                marker="o",
                linewidths=0,
                zorder=2,
            )
    for pos in ticks:
        ax.axvline(pos, lw=0.5, color="0.55")
    ax.axhline(0.0, lw=0.8, ls="--", color="0.25")
    if ticks and labels:
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    ax.set_xlim(float(x[0]), float(x[-1]))
    ax.set_xlabel("K-path")
    ax.set_ylabel(r"Energy - $E_F$ (eV)")
    ax.grid(False)


def find_related_dos_root(band_root: Path) -> Path | None:
    bases = []
    root = band_root.parent if band_root.is_file() else band_root
    bases.extend([root, root.parent, root.parent / "dos", root.parent / "pdos"])
    outdir = find_abacus_outdir(root)
    if outdir:
        bases.extend([outdir, outdir.parent, outdir.parent.parent / "dos"])
    seen: set[Path] = set()
    for base in bases:
        try:
            resolved = base.resolve()
        except OSError:
            resolved = base.absolute()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not base.exists():
            continue
        candidate = resolve_out_path(base)
        if find_first_file(candidate, ["PDOS", "PDOS.dat"], ["PDOS*"]):
            return base
    return None


def load_pdos_groups(
    root: Path,
    explicit_file: Path | None = None,
    selectors: list[str] | None = None,
) -> tuple[Path, np.ndarray, dict[tuple[str, str], np.ndarray]]:
    out_root = resolve_out_path(root)
    pdos_file = explicit_file or find_first_file(out_root, ["PDOS", "PDOS.dat"], ["PDOS*"])
    if not pdos_file:
        die(f"cannot find PDOS file under {out_root}")
    energies, groups = parse_pdos_file(pdos_file)
    selected = parse_selectors(selectors)
    if selected:
        groups = {key: val for key, val in groups.items() if key in selected}
    if not groups:
        die("selected PDOS channels are not present in the PDOS file")
    return pdos_file, energies, groups


def aggregate_pdos_by_orbital(groups: dict[tuple[str, str], np.ndarray]) -> dict[str, np.ndarray]:
    aggregated: dict[str, np.ndarray] = {}
    for (_, orbital), values in groups.items():
        if orbital not in ORBITAL_COLORS:
            continue
        if orbital not in aggregated:
            aggregated[orbital] = np.zeros_like(values)
        aggregated[orbital] += values
    return aggregated


def auto_pdos_xmax(pdos_groups: dict[str, np.ndarray], energies: np.ndarray, emin: float, emax: float) -> float:
    if not pdos_groups:
        return 1.0
    mask = (energies >= emin) & (energies <= emax)
    if not mask.any():
        mask = np.ones_like(energies, dtype=bool)
    values = np.concatenate([np.asarray(channel)[mask].reshape(-1) for channel in pdos_groups.values()])
    values = values[np.isfinite(values) & (values > 0.0)]
    if values.size == 0:
        return 1.0
    xmax = float(values.max())
    return max(xmax * 1.08, 1.0)


def set_energy_ticks(ax, emin: float, emax: float) -> None:
    ax.set_yticks(np.linspace(emin, emax, 5))


def cmd_plot_band(args) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = args.path.expanduser()
    band_data = load_band_plot_data(root, args.file, args.kpt, args.fermi)
    x = band_data["x"]
    bands = band_data["bands"]
    ticks = band_data["ticks"]
    labels = band_data["labels"]
    fermi = band_data["fermi"]
    fermi_log = band_data["fermi_log"]
    band_file = band_data["band_file"]
    block_note = band_data["block_note"]
    gap_info = analyze_band_gap(bands, x, ticks, labels)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=args.dpi)
    pband_file = args.pband_file or find_pband_file(root)
    pband_channels = load_projected_band_channels(pband_file, bands.shape)
    plot_band_axes(ax, x, bands, ticks, labels, pband_channels, linewidth=args.linewidth, color=args.color)
    if not pband_channels:
        input_info = find_input_info_file(root)
        if pband_file is None:
            note = "no projected band file found"
            if input_info and "out_proj_band                  1" not in input_info.read_text(errors="ignore"):
                note += f"; {input_info.name} shows out_proj_band=0"
            print(f"note: {note}; band plot is monochrome")
    ax.set_ylim(args.emin, args.emax)
    set_energy_ticks(ax, args.emin, args.emax)
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    fig.savefig(args.out)
    plt.close(fig)
    source = f" from {fermi_log}" if fermi_log else ""
    projection = f"; orbital colors from {pband_file}" if pband_channels and pband_file else ""
    if block_note:
        print(f"note: {block_note} in {band_file}")
    print(f"wrote BAND plot {args.out}; E_F = {fermi:.9f} eV{source}{projection}")
    print(format_band_gap_message(gap_info))


def cmd_plot_band_pdos(args) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    band_root = args.band_path.expanduser()
    band_data = load_band_plot_data(band_root, args.band_file, args.kpt, args.fermi)
    x = band_data["x"]
    bands = band_data["bands"]
    ticks = band_data["ticks"]
    labels = band_data["labels"]
    fermi = band_data["fermi"]
    fermi_log = band_data["fermi_log"]
    band_file = band_data["band_file"]
    block_note = band_data["block_note"]
    pband_file = args.pband_file or find_pband_file(band_root)
    pband_channels = load_projected_band_channels(pband_file, bands.shape)

    dos_root = args.dos_path.expanduser() if args.dos_path else find_related_dos_root(band_root)
    if dos_root is None:
        die(f"cannot find related PDOS directory for {band_root}; pass --dos-path")
    pdos_file, pdos_energies, pdos_groups = load_pdos_groups(dos_root, args.pdos_file, args.select)
    orbital_pdos = aggregate_pdos_by_orbital(pdos_groups)
    pdos_y = pdos_energies - float(fermi)
    pdos_xmax = float(args.pdos_max) if args.pdos_max is not None else auto_pdos_xmax(orbital_pdos, pdos_y, args.emin, args.emax)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, (band_ax, pdos_ax) = plt.subplots(
        1,
        2,
        figsize=(9.0, 5.0),
        dpi=args.dpi,
        sharey=True,
        constrained_layout=True,
        gridspec_kw={"width_ratios": [3.2, 1.05], "wspace": 0.06},
    )
    plot_band_axes(
        band_ax,
        x,
        bands,
        ticks,
        labels,
        pband_channels,
        linewidth=args.linewidth,
        color=args.color,
    )
    band_ax.set_ylim(args.emin, args.emax)
    set_energy_ticks(band_ax, args.emin, args.emax)
    if args.title:
        band_ax.set_title(args.title)

    for orbital in ("s", "p", "d", "f", "g"):
        values = orbital_pdos.get(orbital)
        if values is None:
            continue
        pdos_ax.plot(values, pdos_y, lw=1.2, color=ORBITAL_COLORS.get(orbital), label=orbital)
    pdos_ax.axhline(0.0, lw=0.8, ls="--", color="0.25")
    pdos_ax.set_xlabel("PDOS")
    pdos_ax.tick_params(axis="y", labelleft=False)
    pdos_ax.grid(False)
    pdos_ax.set_ylim(args.emin, args.emax)
    pdos_ax.set_xlim(0.0, pdos_xmax)
    if not args.no_legend and orbital_pdos:
        pdos_ax.legend(frameon=False, fontsize=8, loc="best")

    fig.savefig(args.out)
    plt.close(fig)
    source = f" from {fermi_log}" if fermi_log else ""
    projection = f"; orbital colors from {pband_file}" if pband_channels and pband_file else ""
    if block_note:
        print(f"note: {block_note} in {band_file}")
    print(f"wrote BAND+PDOS plot {args.out}; E_F = {float(fermi):.9f} eV{source}{projection}; PDOS from {pdos_file}")
    print(format_band_gap_message(analyze_band_gap(bands, x, ticks, labels)))


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
        elif find_first_file(root, ["TDOS.dat", "DOS1_smearing.dat", "DOS1"], ["TDOS*.dat", "DOS*_smearing.dat", "DOS*"]):
            kind = "dos"
        elif find_first_file(root, ["LDOS.txt"], ["LDOS*.txt"]):
            kind = "ldos"
        else:
            kind = "dos"

    if kind == "dos":
        dos_file = args.file or find_first_file(root, ["TDOS.dat", "DOS1_smearing.dat", "DOS1"], ["TDOS*.dat", "DOS*_smearing.dat", "DOS*"])
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
            ax.plot(x, values, lw=1.2, color=ORBITAL_COLORS.get(orbital), label=f"{species}-{orbital}")
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


def grid_contour_levels(
    plane: np.ndarray,
    count: int,
    vmin: float | None = None,
    vmax: float | None = None,
) -> np.ndarray:
    if count < 2:
        die("--levels must be at least 2")
    finite = plane[np.isfinite(plane)]
    if finite.size == 0:
        die("grid slice has no finite values")
    low = float(vmin) if vmin is not None else float(finite.min())
    high = float(vmax) if vmax is not None else float(finite.max())
    if low > high:
        die("--vmin cannot be larger than --vmax")
    if math.isclose(low, high, rel_tol=0.0, abs_tol=1.0e-14):
        pad = max(abs(low) * 1.0e-6, 1.0e-6)
        low -= pad
        high += pad
    return np.linspace(low, high, count)


def plot_grid_slice(
    values: np.ndarray,
    out: Path,
    label: str,
    axis: str,
    index: int | None,
    cmap: str,
    style: str = "image",
    levels: int = 16,
    vmin: float | None = None,
    vmax: float | None = None,
    title: str | None = None,
    contour_color: str = "black",
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/abacuskit-matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plane, used_index = cube_midplane(values, axis, index)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=180)
    data = plane.T
    if style == "image":
        mappable = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        x = np.arange(data.shape[1])
        y = np.arange(data.shape[0])
        xx, yy = np.meshgrid(x, y)
        contour_levels = grid_contour_levels(data, levels, vmin, vmax)
        if style in {"contourf", "both"}:
            mappable = ax.contourf(xx, yy, data, levels=contour_levels, cmap=cmap, extend="both")
        else:
            mappable = ax.contour(xx, yy, data, levels=contour_levels, colors=contour_color, linewidths=0.8)
            ax.clabel(mappable, inline=True, fontsize=7, fmt="%.2g")
        if style == "both":
            line_levels = contour_levels
            contours = ax.contour(xx, yy, data, levels=line_levels, colors=contour_color, linewidths=0.35)
            ax.clabel(contours, inline=True, fontsize=7, fmt="%.2g")
    fig.colorbar(mappable, ax=ax, label=label)
    ax.set_xlabel("grid")
    ax.set_ylabel("grid")
    ax.set_title(title or f"{label}, {axis} slice {used_index}")
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

    plot_grid_slice(
        values,
        args.out,
        label,
        args.axis,
        args.index,
        cmap,
        style=getattr(args, "style", "image"),
        levels=getattr(args, "levels", 16),
        vmin=getattr(args, "vmin", None),
        vmax=getattr(args, "vmax", None),
        title=getattr(args, "title", None),
        contour_color=getattr(args, "contour_color", "black"),
    )
    print(f"wrote {label} plot {args.out}")


def cmd_plot_elf(args) -> None:
    args.kind = "elf"
    args.minus = None
    args.minus_file = None
    args.cube_out = None
    if args.cmap is None:
        args.cmap = "viridis"
    cmd_plot_grid(args)


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


def parse_coordinate_range(text: str) -> tuple[float, float]:
    value = text.strip().replace(" ", "")
    match = re.fullmatch(r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\-([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)", value)
    if not match:
        die("coordinate range must be like a-b, for example 2-3")
    return float(match.group(1)), float(match.group(2))


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
    print("\n[16] Make candidate CIFs\n")
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
    print("\n[15] Prepare ABACUS jobs\n")
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
        "no_kspacing": False,
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
    kspacing = "off (use KPT file)" if state.get("no_kspacing") else state["kspacing"]
    print(
        f"""
---------- 20x: Generate ABACUS INPUT ----------
Current settings:
  Output          : {state["out"]}
  Calculation     : {state["kind"]}
  Suffix          : {state["suffix"]}
  Orbital basis   : {state["orbital_quality"]} LCAO
  Device / solver : {state["device"]} / {state["ks_solver"]}
  kspacing        : {kspacing}
  ecutwfc         : {state["ecutwfc"]}
  nspin           : {state["nspin"]}
  cal_stress      : {state["cal_stress"]}
  DOS/PDOS        : {state["dos"]}
  VDW             : {vdw}
  Dipole corr.    : {dipole}
  DFT+U           : {dftu}
  Extra INPUT     : {format_menu_value(state["set"])}

  201) Generate INPUT now using current settings
  202) Set calculation to scf
  203) Set calculation to relax
  204) Set orbital basis to precision
  205) Set orbital basis to efficiency
  206) Set nspin
  207) Set ecutwfc
  208) Set/enable kspacing
  209) Set device, gpu or cpu
  210) Set ks_solver
  211) Toggle cal_stress
  212) Toggle DOS/PDOS output
  213) Add extra INPUT key=value
  214) Change output INPUT path
  215) Change suffix
  216) Set relax parameters
  217) Clear extra INPUT settings
  218) Reset to defaults
  219) Set VDW correction, e.g. d3_bj
  220) Toggle dipole correction, default Z axis
  221) Apply DOS target template
  222) Apply PDOS target template
  223) Apply band structure target template
  224) Apply COHP output template
  225) Apply work-function/potential template
  226) Enable/edit DFT+U
  227) Apply DFT+U convergence-aid template
  228) Disable DFT+U
  229) Clear DFT+U convergence-aid settings
  230) Apply ELF cube-output template
  231) Apply charge-density cube-output template
  232) Apply hybrid HSE SCF template
  233) Apply hybrid HSE band/NSCF template
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
        state["no_kspacing"] = True
        state["orbital_quality"] = "precision"
        state["orbital_dir"] = None
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_dos", 1)
        set_extra_setting(state, "dos_sigma", 0.07)
        set_extra_setting(state, "dos_edelta_ev", 0.01)
        print("DOS target template applied. Precision LCAO basis is enabled and kspacing is disabled; prepare a dense KPT mesh and previous charge density.")
    elif target == "pdos":
        state["kind"] = "nscf"
        state["dos"] = True
        state["no_kspacing"] = True
        state["orbital_quality"] = "precision"
        state["orbital_dir"] = None
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_dos", 2)
        set_extra_setting(state, "dos_sigma", 0.07)
        set_extra_setting(state, "dos_edelta_ev", 0.01)
        print("PDOS target template applied. Precision LCAO basis is enabled and kspacing is disabled; prepare a dense KPT mesh.")
    elif target == "band":
        state["kind"] = "nscf"
        state["dos"] = False
        state["no_kspacing"] = True
        state["orbital_quality"] = "precision"
        state["orbital_dir"] = None
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_band", 1)
        set_extra_setting(state, "out_proj_band", 1)
        set_extra_setting(state, "smearing_method", "gaussian")
        set_extra_setting(state, "smearing_sigma", 0.02)
        print("Band target template applied. Precision LCAO basis is enabled and kspacing is disabled so ABACUS will use the line-mode KPT file.")
    elif target == "cohp":
        state["kind"] = "scf"
        state["basis_type"] = "lcao"
        state["no_kspacing"] = False
        set_extra_setting(state, "out_mat_hs", "1 8")
        set_extra_setting(state, "out_wfc_lcao", 1)
        set_extra_setting(state, "out_app_flag", 1)
        print("COHP output template applied. Run this LCAO SCF first, then use 132 for COHP post-processing.")
    elif target == "workfunc":
        state["kind"] = "scf"
        state["no_kspacing"] = False
        set_extra_setting(state, "out_pot", 2)
        set_extra_setting(state, "efield_flag", "true")
        set_extra_setting(state, "dip_cor_flag", "true")
        set_extra_setting(state, "efield_dir", 2)
        set_extra_setting(state, "efield_amp", 0)
        print("Work-function/potential template applied. Default dipole correction direction is Z.")
    elif target == "elf":
        state["kind"] = "scf"
        state["no_kspacing"] = False
        set_extra_setting(state, "out_elf", "1 3")
        print("ELF cube-output template applied. ABACUS will write elf.cube under OUT.<suffix>.")
    elif target == "charge":
        state["kind"] = "scf"
        state["no_kspacing"] = False
        set_extra_setting(state, "out_chg", "1 3")
        print("Charge-density cube-output template applied. ABACUS will write charge-density cube files.")
    elif target == "hybrid-scf":
        state["kind"] = "scf"
        state["dos"] = False
        state["no_kspacing"] = False
        state["basis_type"] = "lcao"
        state["orbital_quality"] = "precision"
        state["orbital_dir"] = None
        state["device"] = "cpu"
        state["ks_solver"] = "genelpa"
        set_extra_setting(state, "dft_functional", "hse")
        set_extra_setting(state, "exx_hybrid_step", 100)
        set_extra_setting(state, "out_chg", 1)
        print("Hybrid HSE SCF template applied. Precision LCAO, CPU device, and genelpa solver are enabled; ABACUS must be compiled with LibRI.")
    elif target == "hybrid-band":
        state["kind"] = "nscf"
        state["dos"] = False
        state["no_kspacing"] = True
        state["basis_type"] = "lcao"
        state["orbital_quality"] = "precision"
        state["orbital_dir"] = None
        state["device"] = "cpu"
        state["ks_solver"] = "genelpa"
        set_extra_setting(state, "dft_functional", "hse")
        set_extra_setting(state, "exx_hybrid_step", 100)
        set_extra_setting(state, "kpoint_file", "KPT")
        set_extra_setting(state, "init_chg", "file")
        set_extra_setting(state, "read_file_dir", "./")
        set_extra_setting(state, "out_band", 1)
        set_extra_setting(state, "cal_force", 0)
        set_extra_setting(state, "cal_stress", 0)
        set_extra_setting(state, "smearing_method", "gaussian")
        set_extra_setting(state, "smearing_sigma", 0.02)
        print("Hybrid HSE band/NSCF template applied. Precision LCAO, CPU device, and genelpa solver are enabled; kspacing is disabled for line-mode KPT.")
    else:
        die(f"unknown INPUT target template: {target}")


def interactive_input_template() -> None:
    state = default_input_state()
    while True:
        print_input_menu(state)
        choice = prompt_text("Enter 20x option", "201").lower()
        try:
            if choice in {"q", "quit", "exit"}:
                raise ProgramExit
            if choice in {"0", "200"}:
                return
            if choice == "201":
                run_input_from_state(state)
                print("INPUT generation finished. Exiting abacuskit.")
                raise ProgramExit
            if choice == "202":
                state["kind"] = "scf"
                state["no_kspacing"] = False
                print("Calculation set to scf.")
            elif choice == "203":
                state["kind"] = "relax"
                state["no_kspacing"] = False
                print("Calculation set to relax.")
            elif choice == "204":
                state["orbital_quality"] = "precision"
                print("Orbital basis set to precision.")
            elif choice == "205":
                state["orbital_quality"] = "efficiency"
                print("Orbital basis set to efficiency.")
            elif choice == "206":
                state["nspin"] = prompt_int("nspin", state["nspin"])
            elif choice == "207":
                state["ecutwfc"] = prompt_float("ecutwfc", state["ecutwfc"])
            elif choice == "208":
                state["kspacing"] = prompt_float("kspacing", state["kspacing"])
                state["no_kspacing"] = False
            elif choice == "209":
                state["device"] = prompt_choice("Device", ["gpu", "cpu"], state["device"])
            elif choice == "210":
                state["ks_solver"] = prompt_text("KS solver", state["ks_solver"])
            elif choice == "211":
                state["cal_stress"] = not state["cal_stress"]
                print(f"cal_stress set to {state['cal_stress']}.")
            elif choice == "212":
                state["dos"] = not state["dos"]
                print(f"DOS/PDOS output set to {state['dos']}.")
            elif choice == "213":
                item = prompt_text("Extra INPUT key=value, e.g. mixing_beta=0.05")
                if "=" not in item:
                    die("expected key=value, for example: mixing_beta=0.05")
                append_option(state, "set", item)
            elif choice == "214":
                state["out"] = prompt_path("Output INPUT file", "INPUT")
            elif choice == "215":
                state["suffix"] = prompt_text("ABACUS output suffix", state["suffix"])
            elif choice == "216":
                state["kind"] = "relax"
                state["relax_nmax"] = prompt_int("relax_nmax", state["relax_nmax"])
                state["force_thr_ev"] = prompt_float("force_thr_ev", state["force_thr_ev"])
                state["stress_thr"] = prompt_float("stress_thr", state["stress_thr"])
            elif choice == "217":
                state["set"] = None
                print("Extra INPUT settings cleared.")
            elif choice == "218":
                state.clear()
                state.update(default_input_state())
                print("INPUT settings reset to defaults.")
            elif choice == "219":
                method = prompt_choice("VDW method", ["d3_bj", "d3_0", "d2", "none"], "d3_bj")
                if method == "none":
                    remove_extra_setting(state, "vdw_method")
                    print("VDW correction disabled.")
                else:
                    set_extra_setting(state, "vdw_method", method)
                    print(f"VDW correction set to {method}.")
            elif choice == "220":
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
            elif choice == "221":
                apply_input_target_template(state, "dos")
            elif choice == "222":
                apply_input_target_template(state, "pdos")
            elif choice == "223":
                apply_input_target_template(state, "band")
            elif choice == "224":
                apply_input_target_template(state, "cohp")
            elif choice == "225":
                apply_input_target_template(state, "workfunc")
            elif choice == "226":
                apply_dftu_settings(state)
            elif choice == "227":
                apply_dftu_mixing_aid(state)
            elif choice == "228":
                clear_dftu_settings(state)
            elif choice == "229":
                clear_dftu_mixing_aid(state)
            elif choice == "230":
                apply_input_target_template(state, "elf")
            elif choice == "231":
                apply_input_target_template(state, "charge")
            elif choice == "232":
                apply_input_target_template(state, "hybrid-scf")
            elif choice == "233":
                apply_input_target_template(state, "hybrid-band")
            else:
                print("Unknown 20x option.")
        except SystemExit as exc:
            print(exc)


def interactive_check_abacus() -> None:
    print("\n[4] Check ABACUS job status in current directory\n")
    args = argparse.Namespace(jobs=[Path(".")], json=None, csv=None)
    cmd_check_abacus(args)
    print("ABACUS job check finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_launch_script() -> None:
    print("\n[19] Create ABACUS launch scripts\n")
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
    print(
        """
---------- 30x: Generate ABACUS KPT ----------
  301) Generate regular Gamma mesh KPT using kspacing=0.14
  302) Generate BAND high-symmetry line-mode KPT by SeeK-path
  0) Back to previous menu
  q) Quit abacuskit
"""
    )
    choice = prompt_text("Enter 30x option", "301").lower()
    if choice in {"q", "quit", "exit"}:
        raise ProgramExit
    if choice in {"0", "300"}:
        return
    if choice == "302":
        structure = find_default_structure_file()
        if not structure:
            die("cannot find STRU/POSCAR/CONTCAR/CIF in current directory")
        args = argparse.Namespace(
            structure=structure,
            out=Path("KPT"),
            high_symmetry_points=Path("HIGH_SYMMETRY_POINTS"),
            points_per_segment=20,
            symprec=1.0e-5,
            angle_tolerance=-1.0,
            threshold=1.0e-7,
            no_time_reversal=False,
        )
        cmd_kpt_path(args)
        print("High-symmetry KPT generation finished. Exiting abacuskit.")
        raise ProgramExit
    if choice != "301":
        print("Unknown 30x option.")
        return

    structure = find_default_structure_file()
    if not structure:
        die("cannot find STRU/POSCAR/CONTCAR/CIF in current directory")
    args = argparse.Namespace(
        structure=structure,
        out=Path("KPT"),
        kspacing=0.14,
        mode="mesh",
        high_symmetry_points=Path("HIGH_SYMMETRY_POINTS"),
        points_per_segment=20,
        symprec=1.0e-5,
        angle_tolerance=-1.0,
        threshold=1.0e-7,
        no_time_reversal=False,
    )
    cmd_kpt_auto(args)
    print("KPT generation finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_conv_test() -> None:
    print("\n[17] Prepare convergence-test jobs\n")
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
    print("\n[18] Collect ABACUS metrics / report\n")
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
    print("\n[11] Plot DOS / PDOS / LDOS\n")
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
    print("DOS/PDOS/LDOS plot finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_plot_band() -> None:
    print("\n[12] Plot BAND in current directory\n")
    args = argparse.Namespace(
        path=Path("."),
        file=None,
        pband_file=None,
        kpt=None,
        fermi=None,
        out=Path("band.png"),
        emin=-10.0,
        emax=10.0,
        title=None,
        linewidth=0.8,
        color="C0",
        dpi=300,
    )
    cmd_plot_band(args)
    print("BAND plot finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_plot_band_pdos() -> None:
    print("\n[24] Plot BAND + PDOS in current directory\n")
    args = argparse.Namespace(
        band_path=Path("."),
        band_file=None,
        dos_path=None,
        pdos_file=None,
        pband_file=None,
        kpt=None,
        select=None,
        fermi=None,
        out=Path("band_pdos.png"),
        emin=-10.0,
        emax=10.0,
        title=None,
        linewidth=0.8,
        color="C0",
        dpi=300,
        no_legend=False,
        pdos_max=None,
    )
    cmd_plot_band_pdos(args)
    print("BAND+PDOS plot finished. Exiting abacuskit.")
    raise ProgramExit


def run_interactive_plot_grid(kind: str, default_out: str) -> None:
    args = argparse.Namespace(
        path=Path("."),
        kind=kind,
        file=None,
        minus=None,
        minus_file=None,
        cube_out=None,
        axis="z",
        index=None,
        cmap=None,
        style="contourf" if kind == "elf" else "image",
        levels=16,
        vmin=0.0 if kind == "elf" else None,
        vmax=1.0 if kind == "elf" else None,
        title=None,
        contour_color="black",
        out=Path(default_out),
    )
    if kind == "diff":
        args.minus = prompt_path("Subtracted ABACUS job, OUT.* directory, or cube file")
        minus_file_text = prompt_text("Explicit subtracted cube file, empty for auto", "")
        args.minus_file = Path(minus_file_text).expanduser() if minus_file_text else None
        cube_out_text = prompt_text("Output difference cube", "charge_diff.cube")
        args.cube_out = Path(cube_out_text).expanduser() if cube_out_text else None
    cmd_plot_grid(args)
    print("Grid plot finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_plot_charge() -> None:
    print("\n[8] Plot charge density in current directory\n")
    run_interactive_plot_grid("charge", "charge.png")


def interactive_plot_charge_diff() -> None:
    print("\n[9] Plot charge-density difference in current directory\n")
    run_interactive_plot_grid("diff", "charge_diff.png")


def interactive_plot_elf() -> None:
    print("\n[10] Plot ELF in current directory\n")
    run_interactive_plot_grid("elf", "elf.png")


def interactive_cohp() -> None:
    print(
        """
---------- 13x: ABACUS LCAO COHP ----------
  131) Generate COHP-ready SCF INPUT
  132) List atom orbital channels / global NAO ranges
  133) Calculate COHP/COOP from OUT.ABACUS
  0) Back to previous menu
  q) Quit abacuskit
"""
    )
    choice = prompt_text("Enter 13x option", "131").lower()
    if choice in {"q", "quit", "exit"}:
        raise ProgramExit
    if choice in {"0", "130"}:
        return
    if choice == "131":
        state = default_input_state()
        apply_input_target_template(state, "cohp")
        state["out"] = prompt_path("Output COHP INPUT file", "INPUT.cohp")
        run_input_from_state(state)
        print("COHP INPUT generation finished. Run ABACUS with this INPUT before post-processing.")
        raise ProgramExit
    if choice == "132":
        args = argparse.Namespace(
            out_dir=prompt_path("ABACUS OUT.* directory", "OUT.ABACUS"),
            stru=None,
            input=None,
            orbital_dir=None,
        )
        stru_text = prompt_text("Explicit STRU path, empty for auto", "")
        input_text = prompt_text("Explicit INPUT path, empty for auto", "")
        orbital_text = prompt_text("Explicit orbital_dir, empty for INPUT/default", "")
        args.stru = Path(stru_text).expanduser() if stru_text else None
        args.input = Path(input_text).expanduser() if input_text else None
        args.orbital_dir = Path(orbital_text).expanduser() if orbital_text else None
        cmd_cohp_orbitals(args)
        raise ProgramExit
    if choice == "133":
        use_atom = prompt_yes_no("Use atom index + shell selector?", True)
        args = argparse.Namespace(
            out_dir=prompt_path("ABACUS OUT.* directory", "OUT.ABACUS"),
            atom_i_index=None,
            atom_j_index=None,
            atom_i_orbs=None,
            atom_j_orbs=None,
            stru=None,
            input=None,
            orbital_dir=None,
            method=prompt_choice("Method", ["COHP", "COOP"], "COHP"),
            spin=prompt_choice("Spin channel", ["sum", "up", "down"], "sum"),
            de=prompt_float("Energy grid step de/eV", 0.05),
            no_smooth=not prompt_yes_no("Apply Gaussian smoothing?", True),
            smooth_nstddev=prompt_float("Smoothing sigma in de units", 4.0),
            emin=prompt_float("Plot Emin relative to E_Fermi/eV", -10.0),
            emax=prompt_float("Plot Emax relative to E_Fermi/eV", 10.0),
            width=None,
            invert=prompt_yes_no("Plot -COHP/-COOP convention?", True),
            output_prefix=prompt_path("Output prefix", "COHP"),
        )
        if use_atom:
            args.atom_i_index = prompt_int("Atom I index, 1-based", 1)
            args.atom_j_index = prompt_int("Atom J index, 1-based", 2)
            args.atom_i_orbs = prompt_text("Atom I shells, e.g. all or 3d", "all")
            args.atom_j_orbs = prompt_text("Atom J shells, e.g. all or 2p", "all")
            stru_text = prompt_text("Explicit STRU path, empty for auto", "")
            input_text = prompt_text("Explicit INPUT path, empty for auto", "")
            orbital_text = prompt_text("Explicit orbital_dir, empty for INPUT/default", "")
            args.stru = Path(stru_text).expanduser() if stru_text else None
            args.input = Path(input_text).expanduser() if input_text else None
            args.orbital_dir = Path(orbital_text).expanduser() if orbital_text else None
        else:
            args.atom_i_orbs = prompt_text("Atom/group I global NAO indices, comma-separated")
            args.atom_j_orbs = prompt_text("Atom/group J global NAO indices, comma-separated")
        cmd_cohp(args)
        raise ProgramExit
    print("Unknown 13x option.")


def interactive_bader() -> None:
    print("\n[14] Bader charge analysis\n")
    bader_text = prompt_text("bader executable, empty for PATH/ABACUSKIT_BADER", "")
    cube_text = prompt_text("Explicit charge cube, empty for auto", "")
    ref_text = prompt_text("Reference cube for bader -ref, empty for none", "")
    args = argparse.Namespace(
        path=prompt_path("ABACUS job, OUT.* directory, or cube file", "."),
        cube=Path(cube_text).expanduser() if cube_text else None,
        reference_cube=Path(ref_text).expanduser() if ref_text else None,
        bader=Path(bader_text).expanduser() if bader_text else None,
        work_dir=prompt_path("Bader work directory", "bader_work"),
        total_cube=None,
        bader_arg=[],
        out=prompt_path("CSV output", "bader.csv"),
        json=prompt_path("JSON output", "bader.json"),
        no_json=False,
    )
    cmd_bader(args)
    print("Bader analysis finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_collect_deepmd() -> None:
    print("\n[20] Collect ABACUS outputs to DeepMD data\n")
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
    print("\n[21] Make DeepMD training input\n")
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
    print("\n[22] Init workflow skeleton\n")
    args = argparse.Namespace(out=prompt_path("Workflow root directory", "abacus_deepmd_project"))
    cmd_init_workflow(args)


def apns_resource_status() -> dict[str, Path | None]:
    return {
        "pseudo_dir": find_apns_dir(APNS_RESOURCE_NAMES["pseudo"]),
        "orbital_efficiency_dir": find_apns_dir(APNS_RESOURCE_NAMES["orbital_efficiency"]),
        "orbital_precision_dir": find_apns_dir(APNS_RESOURCE_NAMES["orbital_precision"]),
    }


def print_apns_resource_status(resources: dict[str, Path | None]) -> None:
    labels = {
        "pseudo_dir": "Pseudopotentials",
        "orbital_efficiency_dir": "Orbitals efficiency",
        "orbital_precision_dir": "Orbitals precision",
    }
    print("\nDetected APNS resource paths:")
    for key, label in labels.items():
        value = resources.get(key)
        print(f"  {label:<20}: {value if value else 'not found'}")


def save_user_config(updates: dict[str, str]) -> None:
    config = dict(USER_CONFIG)
    config.update(updates)
    USER_CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    USER_CONFIG_PATH.write_text(json.dumps(config, indent=2) + "\n")


def interactive_config_apns_paths() -> None:
    print("\n[23] Search and save APNS pseudopotential/orbital paths\n")
    resources = apns_resource_status()
    print_apns_resource_status(resources)

    missing = [key for key, path in resources.items() if not path or not Path(path).expanduser().is_dir()]
    if missing:
        die("cannot save APNS config; missing valid directories: " + ", ".join(missing))

    save_user_config({key: str(Path(path).expanduser()) for key, path in resources.items() if path})
    print(f"Saved APNS paths to {USER_CONFIG_PATH}")
    print("Future abacuskit runs will use this config unless command-line arguments or environment variables override it.")
    raise ProgramExit


def cmd_fix_stru_range(args) -> None:
    fixed = fix_stru_range(
        stru=args.stru,
        out=args.out or args.stru,
        axis=args.axis,
        lower=args.min,
        upper=args.max,
        backup=not args.no_backup,
    )
    target = args.out or args.stru
    print(f"fixed {fixed} atoms in xyz directions; wrote {target}")


def cmd_rotate_vacuum_z_to_y(args) -> None:
    rotate_stru_vacuum_z_to_y(
        stru=args.stru,
        out=args.out or args.stru,
        backup=not args.no_backup,
    )
    target = args.out or args.stru
    print(f"rotated vacuum direction from Z to Y; wrote {target}")


def cmd_stru2cif(args) -> None:
    atoms = canonicalize_axis_aligned_atoms(parse_stru_atoms(args.stru))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write(args.out, atoms, format="cif")
    print(f"converted {args.stru} to CIF {args.out}")


def cmd_cohp_orbitals(args) -> None:
    out_dir = args.out_dir
    stru = args.stru or (out_dir.parent / "STRU")
    inp = args.input or (out_dir.parent / "INPUT")
    if not stru.is_file():
        die(f"cannot find STRU for orbital mapping: {stru}")
    orbital_map = build_orbital_map(
        stru_path=stru,
        input_path=inp if inp.is_file() else None,
        orbital_dir=args.orbital_dir,
    )
    print(format_orbital_map(orbital_map))


def cmd_cohp(args) -> None:
    try:
        atom_i_orbs, atom_j_orbs, resolution = resolve_orbital_arguments(
            out_dir=args.out_dir,
            atom_i_orbs=args.atom_i_orbs,
            atom_j_orbs=args.atom_j_orbs,
            atom_i_index=args.atom_i_index,
            atom_j_index=args.atom_j_index,
            stru_path=args.stru,
            input_path=args.input,
            orbital_dir=args.orbital_dir,
        )
        metadata = run_cohp(
            out_dir=args.out_dir,
            atom_i_orbs=atom_i_orbs,
            atom_j_orbs=atom_j_orbs,
            method=args.method,
            spin=args.spin,
            de=args.de,
            smooth=not args.no_smooth,
            smooth_nstddev=args.smooth_nstddev,
            invert=args.invert,
            output_prefix=args.output_prefix,
            emin=args.emin,
            emax=args.emax,
            width=args.width,
        )
    except (OSError, ValueError) as exc:
        die(str(exc))

    if resolution is not None:
        sel_i, sel_j, orbital_map = resolution
        print(
            f"atom I: {sel_i.symbol} #{sel_i.atom_index}, {sel_i.selector} -> "
            f"{len(sel_i.indices)} NAOs"
        )
        print(
            f"atom J: {sel_j.symbol} #{sel_j.atom_index}, {sel_j.selector} -> "
            f"{len(sel_j.indices)} NAOs"
        )
        print(f"orbital map: {len(orbital_map.atoms)} atoms, {orbital_map.total_orbitals} NAOs")
    print(f"E_Fermi = {metadata['efermi_ev']:.6f} eV")
    print(f"ICOHP = {metadata['icohp_raw_ev']:.6f} eV")
    print(f"-ICOHP = {metadata['minus_icohp_ev']:.6f} eV")
    print(f"raw COHP: {metadata['files']['raw']}")
    print(f"E-E_Fermi COHP: {metadata['files']['shifted']}")
    print(f"ICOHP summary: {metadata['files']['icohp']}")
    if metadata["files"]["plot"]:
        print(f"plot: {metadata['files']['plot']}")
    else:
        print("plot: skipped because matplotlib is not installed")
    print(f"metadata: {metadata['files']['metadata']}")


def cmd_bader(args) -> None:
    path = resolve_out_path(args.path)
    bader_program = args.bader if args.bader else None
    json_out = None if args.no_json else args.json
    try:
        metadata = run_bader_analysis(
            path=path,
            cube=args.cube,
            work_dir=args.work_dir,
            total_cube=args.total_cube,
            bader=bader_program,
            reference_cube=args.reference_cube,
            extra_args=args.bader_arg,
        )
        rows = metadata["rows"]
        write_bader_csv(args.out, rows)
        if json_out:
            write_bader_json(json_out, metadata)
    except (OSError, RuntimeError, ValueError) as exc:
        die(str(exc))

    total_charge = sum(float(row["charge"]) for row in rows)
    total_bader = sum(float(row["bader_electrons"]) for row in rows)
    print(f"Bader atoms: {len(rows)}")
    print(f"Bader electrons: {total_bader:.8f}")
    print(f"Net valence charge sum: {total_charge:.8f}")
    print(f"charge cube: {metadata['files']['charge_cube']}")
    if metadata["generated_total_cube"]:
        print("spin channels: summed SPIN1_CHG.cube + SPIN2_CHG.cube")
    print(f"ACF.dat: {metadata['files']['acf']}")
    print(f"CSV: {args.out}")
    if json_out:
        print(f"JSON: {json_out}")


def interactive_fix_stru_range() -> None:
    print("\n[5] Fix STRU atoms by coordinate range\n")
    stru = Path("STRU")
    if not stru.is_file():
        die("cannot find STRU in current directory")
    axis = prompt_choice("Coordinate axis", ["z", "x", "y"], "z")
    range_text = prompt_text(f"{axis} coordinate range a-b", "0-0")
    lower, upper = parse_coordinate_range(range_text)
    args = argparse.Namespace(stru=stru, out=None, axis=axis, min=lower, max=upper, no_backup=False)
    cmd_fix_stru_range(args)
    print("STRU coordinate-range fixing finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_rotate_vacuum_z_to_y() -> None:
    print("\n[6] Rotate STRU vacuum direction from Z to Y\n")
    stru = Path("STRU")
    if not stru.is_file():
        die("cannot find STRU in current directory")
    args = argparse.Namespace(stru=stru, out=None, no_backup=False)
    cmd_rotate_vacuum_z_to_y(args)
    print("STRU vacuum-direction rotation finished. Exiting abacuskit.")
    raise ProgramExit


def interactive_stru2cif() -> None:
    print("\n[7] Convert STRU to CIF\n")
    stru = Path("STRU")
    if not stru.is_file():
        die("cannot find STRU in current directory")
    args = argparse.Namespace(stru=stru, out=Path("STRU.cif"))
    cmd_stru2cif(args)
    print("STRU -> CIF conversion finished. Exiting abacuskit.")
    raise ProgramExit


def print_interactive_menu() -> None:
    print_terminal_logo()
    print(
        f"""
============== abacuskit {__version__} ==============
Author: {__author__}
Affiliation: {__affiliation__}

  1) CIF -> ABACUS STRU
  2) Generate ABACUS INPUT
  3) Generate ABACUS KPT
  4) Check ABACUS job status
  5) Fix STRU atoms by coordinate range
  6) Rotate STRU vacuum direction Z -> Y
  7) Convert STRU to CIF

  8) Plot charge density
  9) Plot charge-density difference
  10) Plot ELF
  11) Plot DOS / PDOS / LDOS
  12) Auto plot BAND for current directory
  13) ABACUS LCAO COHP
  14) Bader charge analysis

  15) Prepare ABACUS jobs
  16) Make candidate CIFs
  17) Prepare convergence-test jobs
  18) Collect ABACUS metrics / report
  19) Create ABACUS launch scripts
  20) Collect ABACUS outputs to DeepMD data
  21) Make DeepMD training input
  22) Init workflow skeleton

  23) Search/save APNS pseudopotential and orbital paths
  24) Plot BAND + PDOS in current directory
  h) Show command-line help
  q) Quit abacuskit
  0) Exit
"""
    )


def interactive_menu() -> None:
    actions = {
        "1": interactive_cif2stru,
        "2": interactive_input_template,
        "3": interactive_kpt,
        "4": interactive_check_abacus,
        "5": interactive_fix_stru_range,
        "6": interactive_rotate_vacuum_z_to_y,
        "7": interactive_stru2cif,
        "8": interactive_plot_charge,
        "9": interactive_plot_charge_diff,
        "10": interactive_plot_elf,
        "11": interactive_plot_dos,
        "12": interactive_plot_band,
        "13": interactive_cohp,
        "14": interactive_bader,
        "15": interactive_prepare_abacus,
        "16": interactive_make_candidates,
        "17": interactive_conv_test,
        "18": interactive_collect_report,
        "19": interactive_launch_script,
        "20": interactive_collect_deepmd,
        "21": interactive_make_train,
        "22": interactive_init_workflow,
        "23": interactive_config_apns_paths,
        "24": interactive_plot_band_pdos,
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

    p = sub.add_parser("fix-stru-range", help="fix STRU atoms whose coordinate falls in an axis range")
    p.add_argument("stru", type=Path, nargs="?", default=Path("STRU"))
    p.add_argument("--axis", choices=["x", "y", "z"], required=True)
    p.add_argument("--min", type=float, required=True)
    p.add_argument("--max", type=float, required=True)
    p.add_argument("--out", type=Path, help="output STRU; default overwrites input and writes a .bak backup")
    p.add_argument("--no-backup", action="store_true", help="do not write backup when overwriting input")
    p.set_defaults(func=cmd_fix_stru_range)

    p = sub.add_parser("rotate-vacuum-z-to-y", help="rotate STRU so the Z vacuum direction becomes Y")
    p.add_argument("stru", type=Path, nargs="?", default=Path("STRU"))
    p.add_argument("--out", type=Path, help="output STRU; default overwrites input and writes a .bak backup")
    p.add_argument("--no-backup", action="store_true", help="do not write backup when overwriting input")
    p.set_defaults(func=cmd_rotate_vacuum_z_to_y)

    p = sub.add_parser("stru2cif", help="convert an ABACUS STRU file to CIF")
    p.add_argument("stru", type=Path, nargs="?", default=Path("STRU"))
    p.add_argument("-o", "--out", type=Path, default=Path("STRU.cif"))
    p.set_defaults(func=cmd_stru2cif)

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

    p = sub.add_parser("input-template", help="write an ABACUS INPUT template")
    p.add_argument("--kind", choices=["scf", "relax", "nscf"], required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--suffix", default="ABACUS")
    p.add_argument("--pseudo-dir", type=Path, default=DEFAULT_PSEUDO_DIR)
    p.add_argument("--orbital-dir", type=Path)
    p.add_argument("--orbital-quality", choices=sorted(DEFAULT_ORBITAL_DIRS), default="efficiency")
    p.add_argument("--basis-type", choices=["lcao", "pw"], default="lcao")
    p.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    p.add_argument("--ks-solver", default="cusolver")
    p.add_argument("--kspacing", type=float, default=0.14)
    p.add_argument("--no-kspacing", action="store_true", help="do not write kspacing, so ABACUS uses the KPT file")
    p.add_argument("--ecutwfc", type=float, default=100)
    p.add_argument("--nspin", type=int, default=1)
    p.add_argument("--cal-stress", action="store_true")
    p.add_argument("--relax-nmax", type=int, default=100)
    p.add_argument("--force-thr-ev", type=float, default=0.04)
    p.add_argument("--stress-thr", type=float, default=1.0)
    p.add_argument("--dos", action="store_true", help="include DOS/PDOS output parameters")
    p.add_argument("--hybrid-hse-scf", action="store_true", help="apply the LCAO HSE SCF template used by menu option 232")
    p.add_argument("--hybrid-hse-band", action="store_true", help="apply the LCAO HSE band/NSCF template used by menu option 233")
    p.add_argument("--set", action="append", help="extra INPUT key=value; can be repeated")
    p.add_argument("--no-comments", action="store_true")
    p.set_defaults(func=cmd_input_template)

    p = sub.add_parser("kpt", help="write an ABACUS KPT file")
    p.add_argument("--mesh", type=int, nargs=3, required=True, metavar=("NX", "NY", "NZ"))
    p.add_argument("--shift", type=int, nargs=3, default=(0, 0, 0), metavar=("SX", "SY", "SZ"))
    p.add_argument("--model", choices=["gamma", "mp"], default="gamma")
    p.add_argument("--out", type=Path, default=Path("KPT"))
    p.set_defaults(func=cmd_kpt)

    p = sub.add_parser("kpt-path", help="write a VASPKIT-like high-symmetry line-mode KPT using SeeK-path")
    p.add_argument("structure", type=Path, help="structure file readable by ASE, e.g. STRU, CIF, or POSCAR")
    p.add_argument("--out", type=Path, default=Path("KPT"))
    p.add_argument("--high-symmetry-points", type=Path, default=Path("HIGH_SYMMETRY_POINTS"))
    p.add_argument("--points-per-segment", type=int, default=20)
    p.add_argument("--symprec", type=float, default=1.0e-5)
    p.add_argument("--angle-tolerance", type=float, default=-1.0)
    p.add_argument("--threshold", type=float, default=1.0e-7)
    p.add_argument("--no-time-reversal", action="store_true")
    p.set_defaults(func=cmd_kpt_path)

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

    p = sub.add_parser("plot-band", help="plot ABACUS BANDS_*.dat and shift Fermi level to 0 eV")
    p.add_argument("path", type=Path, help="ABACUS job directory, OUT.* directory, or BANDS_*.dat file")
    p.add_argument("--file", type=Path, help="explicit BANDS_*.dat file")
    p.add_argument("--pband-file", type=Path, help="explicit projected band XML file, e.g. pbands1.xml")
    p.add_argument("--kpt", type=Path, help="explicit line-mode KPT file for high-symmetry labels")
    p.add_argument("--fermi", type=float, help="explicit Fermi energy in eV; default reads running_*.log")
    p.add_argument("--out", type=Path, default=Path("band.png"))
    p.add_argument("--emin", type=float, default=-10.0, help="minimum plotted energy after Fermi shift")
    p.add_argument("--emax", type=float, default=10.0, help="maximum plotted energy after Fermi shift")
    p.add_argument("--title", help="plot title")
    p.add_argument("--linewidth", type=float, default=0.8)
    p.add_argument("--color", default="C0")
    p.add_argument("--dpi", type=int, default=300)
    p.set_defaults(func=cmd_plot_band)

    p = sub.add_parser("plot-band-pdos", help="plot band structure with vertical PDOS")
    p.add_argument("band_path", type=Path, help="ABACUS band job directory, OUT.* directory, or BANDS_*.dat file")
    p.add_argument("--band-file", type=Path, help="explicit BANDS_*.dat or band.txt file")
    p.add_argument("--dos-path", type=Path, help="ABACUS DOS/PDOS job directory or OUT.* directory")
    p.add_argument("--pdos-file", type=Path, help="explicit PDOS file")
    p.add_argument("--pband-file", type=Path, help="explicit projected band XML file, e.g. pbands1.xml")
    p.add_argument("--kpt", type=Path, help="explicit line-mode KPT file for high-symmetry labels")
    p.add_argument("--select", action="append", help="PDOS selector, e.g. C=p --select H=s --select O=p --select Ni=d")
    p.add_argument("--fermi", type=float, help="explicit Fermi energy in eV; default reads running_*.log")
    p.add_argument("--out", type=Path, default=Path("band_pdos.png"))
    p.add_argument("--emin", type=float, default=-10.0, help="minimum plotted energy after Fermi shift")
    p.add_argument("--emax", type=float, default=10.0, help="maximum plotted energy after Fermi shift")
    p.add_argument("--title", help="plot title")
    p.add_argument("--linewidth", type=float, default=0.8)
    p.add_argument("--color", default="C0")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--no-legend", action="store_true", help="hide PDOS legend")
    p.add_argument("--pdos-max", type=float, help="manual PDOS x-axis maximum; default auto-scales from plotted window")
    p.set_defaults(func=cmd_plot_band_pdos)

    p = sub.add_parser("cohp-orbitals", help="list atom shell channels and global NAO ranges for COHP")
    p.add_argument("out_dir", type=Path, nargs="?", default=Path("OUT.ABACUS"))
    p.add_argument("--stru", type=Path, help="ABACUS STRU path; default searches beside OUT.*")
    p.add_argument("--input", type=Path, help="ABACUS INPUT path; default searches beside OUT.*")
    p.add_argument("--orbital-dir", type=Path, help="directory containing ABACUS numerical orbital files")
    p.set_defaults(func=cmd_cohp_orbitals)

    p = sub.add_parser("cohp", help="calculate built-in ABACUS LCAO COHP/COOP from OUT.* outputs")
    p.add_argument("out_dir", type=Path, nargs="?", default=Path("OUT.ABACUS"))
    p.add_argument("--atom-i-orbs", required=True, help="global NAO indices, or shell selector with --atom-i-index")
    p.add_argument("--atom-j-orbs", required=True, help="global NAO indices, or shell selector with --atom-j-index")
    p.add_argument("--atom-i-index", type=int, help="1-based atom index for atom/group I")
    p.add_argument("--atom-j-index", type=int, help="1-based atom index for atom/group J")
    p.add_argument("--stru", type=Path, help="ABACUS STRU path for atom-index shell selection")
    p.add_argument("--input", type=Path, help="ABACUS INPUT path for orbital_dir discovery")
    p.add_argument("--orbital-dir", type=Path, help="directory containing ABACUS numerical orbital files")
    p.add_argument("--method", choices=["COHP", "COOP"], default="COHP")
    p.add_argument("--spin", choices=["sum", "up", "down"], default="sum")
    p.add_argument("--de", type=float, default=0.05)
    p.add_argument("--no-smooth", action="store_true")
    p.add_argument("--smooth-nstddev", type=float, default=4.0)
    p.add_argument("--emin", type=float, default=-10.0)
    p.add_argument("--emax", type=float, default=10.0)
    p.add_argument("--width", type=float)
    p.add_argument("--invert", action="store_true", default=True, help="plot -COHP/-COOP convention")
    p.add_argument("--no-invert", action="store_false", dest="invert", help="plot raw COHP/COOP sign")
    p.add_argument("--output-prefix", type=Path, default=Path("COHP"))
    p.set_defaults(func=cmd_cohp)

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
    p.add_argument("--style", choices=["image", "contour", "contourf", "both"], default="image", help="2D slice style")
    p.add_argument("--levels", type=int, default=16, help="number of contour levels for contour/contourf/both")
    p.add_argument("--vmin", type=float, help="minimum plotted value for color scale/contour levels")
    p.add_argument("--vmax", type=float, help="maximum plotted value for color scale/contour levels")
    p.add_argument("--title", help="plot title")
    p.add_argument("--contour-color", default="black", help="contour line color for contour/both styles")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_plot_grid)

    p = sub.add_parser("plot-elf", help="plot 2D ELF cube slices as contour or filled contour maps")
    p.add_argument("path", type=Path, help="ABACUS job directory, OUT.* directory, or ELF cube file")
    p.add_argument("--file", type=Path, help="explicit ELF cube file")
    p.add_argument("--axis", choices=["x", "y", "z"], default="z")
    p.add_argument("--index", type=int, help="slice index; default is the middle slice")
    p.add_argument("--style", choices=["contourf", "contour", "both", "image"], default="contourf", help="ELF 2D plot style")
    p.add_argument("--levels", type=int, default=16, help="number of contour levels")
    p.add_argument("--vmin", type=float, default=0.0, help="minimum plotted ELF value")
    p.add_argument("--vmax", type=float, default=1.0, help="maximum plotted ELF value")
    p.add_argument("--cmap", default="viridis", help="matplotlib colormap name")
    p.add_argument("--title", help="plot title")
    p.add_argument("--contour-color", default="black", help="contour line color for contour/both styles")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=cmd_plot_elf)

    p = sub.add_parser("bader", help="calculate Bader charges from ABACUS charge-density cube output")
    p.add_argument("path", type=Path, nargs="?", default=Path("."), help="ABACUS job directory, OUT.* directory, or cube file")
    p.add_argument("--cube", type=Path, help="explicit charge-density cube; otherwise auto-detect CHG/SPIN*_CHG")
    p.add_argument("--reference-cube", "--ref", type=Path, help="optional reference cube passed to bader as -ref")
    p.add_argument("--bader", type=Path, help="bader executable path; default uses ABACUSKIT_BADER or PATH")
    p.add_argument("--work-dir", type=Path, default=Path("bader_work"), help="directory where bader writes ACF.dat")
    p.add_argument("--total-cube", type=Path, help="output path for summed spin cube; default is work-dir/TOTAL_CHG.cube")
    p.add_argument("--bader-arg", action="append", default=[], help="extra argument passed to bader; repeat as needed")
    p.add_argument("--out", type=Path, default=Path("bader.csv"), help="CSV summary output")
    p.add_argument("--json", type=Path, default=Path("bader.json"), help="JSON metadata output")
    p.add_argument("--no-json", action="store_true", help="do not write JSON metadata")
    p.set_defaults(func=cmd_bader)

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
