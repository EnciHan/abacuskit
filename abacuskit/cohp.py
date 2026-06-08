"""ABACUS LCAO COHP post-processing helpers."""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from pathlib import Path

import numpy as np

RY_TO_EV = 13.605693009
SHELL_ORDER = ["s", "p", "d", "f", "g"]
SHELL_MULTIPLICITY = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}


@dataclass
class OrbitalAtom:
    atom_index: int
    symbol: str
    orbital_start: int
    orbital_stop: int
    shells: dict[str, list[int]]


@dataclass
class OrbitalMap:
    atoms: list[OrbitalAtom]
    total_orbitals: int
    stru_path: Path
    input_path: Path | None
    orbital_dir: Path | None


@dataclass
class OrbitalSelection:
    atom_index: int
    symbol: str
    selector: str
    normalized_selectors: list[str]
    indices: list[int]


def _first_existing(paths):
    for path in paths:
        if path is not None and Path(path).exists():
            return Path(path)
    return None


def _strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def _section_lines(lines: list[str], section_name: str, keep_blank: bool = False) -> list[str]:
    headers = {
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "ATOMIC_POSITIONS",
    }
    start = None
    for idx, line in enumerate(lines):
        if _strip_comment(line).upper() == section_name:
            start = idx + 1
            break
    if start is None:
        return []

    out: list[str] = []
    for line in lines[start:]:
        cleaned = _strip_comment(line)
        if not cleaned:
            if keep_blank:
                continue
            if out:
                break
            continue
        if cleaned.upper() in headers:
            break
        out.append(cleaned)
    return out


def cxx_complex(text: str) -> complex:
    token = text.strip()
    if token.startswith("(") and token.endswith(")") and "," in token:
        real, imag = token[1:-1].split(",", 1)
        return complex(float(real), float(imag))
    return complex(float(token), 0.0)


def read_hs_matrix(path: Path) -> np.ndarray:
    tokens = Path(path).read_text(errors="ignore").replace("\n", " ").split()
    if not tokens:
        raise ValueError(f"empty matrix file: {path}")
    size = int(tokens[0])
    expected = size * (size + 1) // 2
    values = tokens[1:]
    if len(values) < expected:
        raise ValueError(f"matrix file {path} has {len(values)} values, expected at least {expected}")

    matrix = np.zeros((size, size), dtype=np.complex128)
    cursor = 0
    for i in range(size):
        for j in range(i, size):
            matrix[i, j] = cxx_complex(values[cursor])
            cursor += 1
    matrix = matrix + matrix.conj().T
    matrix[np.diag_indices(size)] *= 0.5
    return matrix


def read_wfc_nao(path: Path) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
    lines = [line.strip() for line in Path(path).read_text(errors="ignore").splitlines()]
    nband = 0
    nlocal = 0
    for line in lines:
        if line.endswith("(number of bands)"):
            nband = int(line.split()[0])
        elif line.endswith("(number of orbitals)"):
            nlocal = int(line.split()[0])
    if nband <= 0 or nlocal <= 0:
        raise ValueError(f"cannot read band/orbital count from {path}")

    coeffs = np.zeros((nlocal, nband), dtype=np.complex128)
    energies = np.zeros(nband, dtype=float)
    occupations = np.zeros(nband, dtype=float)
    kvec = None
    band = -1
    raw: list[float] = []

    def flush() -> None:
        nonlocal raw
        if band < 0 or not raw:
            raw = []
            return
        if len(raw) == nlocal:
            coeffs[:, band] = np.asarray(raw, dtype=float)
        elif len(raw) == 2 * nlocal:
            coeffs[:, band] = np.asarray(
                [complex(raw[i], raw[i + 1]) for i in range(0, len(raw), 2)],
                dtype=np.complex128,
            )
        else:
            raise ValueError(
                f"unexpected coefficient count in {path}: band {band + 1} has {len(raw)}, "
                f"expected {nlocal} or {2 * nlocal}"
            )
        raw = []

    for line in lines:
        if not line:
            continue
        if line.endswith("(band)"):
            flush()
            band = int(line.split()[0]) - 1
            continue
        if line.endswith("(Ry)") and band >= 0:
            energies[band] = float(line.split()[0])
            continue
        if line.endswith("(Occupations)") and band >= 0:
            occupations[band] = float(line.split()[0])
            continue
        if line.endswith(")"):
            continue
        fields = line.split()
        if band < 0 and len(fields) == 3:
            try:
                kvec = np.asarray([float(x) for x in fields], dtype=float)
            except ValueError:
                pass
        elif band >= 0:
            raw.extend(float(x) for x in fields)
    flush()
    return coeffs, kvec, energies, occupations


def read_kpoint_weights(path: Path, nmatrix: int) -> np.ndarray:
    path = Path(path)
    if not path.exists():
        return np.full(nmatrix, 1.0 / nmatrix)
    rows: list[list[float]] = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if not fields or not fields[0][0].isdigit():
            continue
        try:
            rows.append([float(x) for x in fields])
        except ValueError:
            continue
    if not rows:
        return np.full(nmatrix, 1.0 / nmatrix)
    first_width = len(rows[0])
    table = [row for row in rows if len(row) == first_width]
    weights = np.asarray([row[-1] for row in table], dtype=float)
    if len(weights) == 0:
        return np.full(nmatrix, 1.0 / nmatrix)
    return weights


def read_fermi_energy(path: Path) -> float:
    path = Path(path)
    if not path.exists():
        return 0.0
    values: list[float] = []
    patterns = [
        re.compile(r"^\s*E_Fermi\s*[:=]?\s*([-+0-9.eE]+)"),
        re.compile(r"\bEFERMI\b\s*[:=]?\s*([-+0-9.eE]+)", re.IGNORECASE),
        re.compile(r"\bFermi\b.*?([-+0-9.eE]+)\s*eV", re.IGNORECASE),
    ]
    for line in path.read_text(errors="ignore").splitlines():
        for pattern in patterns:
            match = pattern.search(line)
            if match:
                values.append(float(match.group(1)))
                break
    return values[-1] if values else 0.0


def parse_input_orbital_dir(input_path: Path | None) -> Path | None:
    if input_path is None or not Path(input_path).exists():
        return None
    input_path = Path(input_path)
    for line in input_path.read_text(errors="ignore").splitlines():
        fields = _strip_comment(line).split()
        if len(fields) >= 2 and fields[0].lower() == "orbital_dir":
            path = Path(fields[1]).expanduser()
            return path if path.is_absolute() else (input_path.parent / path).resolve()
    return None


def parse_stru_metadata(stru_path: Path) -> tuple[list[str], dict[str, str], list[str]]:
    stru_path = Path(stru_path)
    lines = stru_path.read_text(errors="ignore").splitlines()
    species = [line.split()[0] for line in _section_lines(lines, "ATOMIC_SPECIES") if line.split()]
    if not species:
        raise ValueError(f"cannot parse ATOMIC_SPECIES from {stru_path}")
    orbital_files = _section_lines(lines, "NUMERICAL_ORBITAL")
    if len(orbital_files) != len(species):
        raise ValueError(
            f"NUMERICAL_ORBITAL entries ({len(orbital_files)}) do not match species ({len(species)}) in {stru_path}"
        )

    position_lines = _section_lines(lines, "ATOMIC_POSITIONS", keep_blank=True)
    if len(position_lines) < 4:
        raise ValueError(f"cannot parse ATOMIC_POSITIONS from {stru_path}")
    atoms: list[str] = []
    cursor = 1
    while cursor < len(position_lines):
        symbol = position_lines[cursor].split()[0]
        if cursor + 2 >= len(position_lines):
            raise ValueError(f"incomplete ATOMIC_POSITIONS block for {symbol} in {stru_path}")
        count = int(float(position_lines[cursor + 2].split()[0]))
        atoms.extend([symbol] * count)
        cursor += 3 + count
    return species, dict(zip(species, orbital_files)), atoms


def shell_counts_from_orbital_file(path: Path) -> dict[str, int]:
    text = Path(path).read_text(errors="ignore")
    shells: dict[str, int] = {}
    for shell in SHELL_ORDER:
        match = re.search(rf"Number\s+of\s+{shell.upper()}orbital\s*-->\s*(\d+)", text, re.IGNORECASE)
        if match:
            shells[shell] = int(match.group(1))
    return shells


def shell_counts_from_filename(path: Path) -> dict[str, int]:
    name = Path(path).name.lower()
    return {shell: int(count) for count, shell in re.findall(r"(\d+)([spdfg])", name)}


def parse_orbital_shells(path: Path) -> dict[str, int]:
    path = Path(path)
    if path.exists():
        shells = shell_counts_from_orbital_file(path)
        if shells:
            return shells
    shells = shell_counts_from_filename(path)
    if shells:
        return shells
    raise ValueError(f"cannot determine orbital shell counts from {path}")


def resolve_orbital_path(entry: str, stru_path: Path, orbital_dir: Path | None) -> Path:
    path = Path(entry).expanduser()
    candidates = [path] if path.is_absolute() else [(Path(stru_path).parent / path).resolve()]
    if not path.is_absolute() and orbital_dir is not None:
        candidates.extend([(Path(orbital_dir) / path.name).resolve(), (Path(orbital_dir) / path).resolve()])
    return _first_existing(candidates) or candidates[0]


def shell_ranges(start: int, shell_counts: dict[str, int]) -> dict[str, list[int]]:
    ranges: dict[str, list[int]] = {}
    cursor = start
    for shell in SHELL_ORDER:
        if shell not in shell_counts:
            continue
        width = shell_counts[shell] * SHELL_MULTIPLICITY[shell]
        ranges[shell] = list(range(cursor, cursor + width))
        cursor += width
    ranges["all"] = list(range(start, cursor))
    return ranges


def build_orbital_map(stru_path: Path, input_path: Path | None = None, orbital_dir: Path | None = None) -> OrbitalMap:
    stru_path = Path(stru_path)
    input_path = Path(input_path) if input_path is not None and Path(input_path).exists() else None
    orbital_dir = Path(orbital_dir).expanduser() if orbital_dir else parse_input_orbital_dir(input_path)
    if orbital_dir is not None and not orbital_dir.is_absolute():
        orbital_dir = (stru_path.parent / orbital_dir).resolve()

    _, orbital_by_symbol, atom_symbols = parse_stru_metadata(stru_path)
    shells_by_symbol = {}
    for symbol, entry in orbital_by_symbol.items():
        shells_by_symbol[symbol] = parse_orbital_shells(resolve_orbital_path(entry, stru_path, orbital_dir))

    atoms: list[OrbitalAtom] = []
    cursor = 0
    for atom_index, symbol in enumerate(atom_symbols, start=1):
        ranges = shell_ranges(cursor, shells_by_symbol[symbol])
        atoms.append(
            OrbitalAtom(
                atom_index=atom_index,
                symbol=symbol,
                orbital_start=cursor,
                orbital_stop=cursor + len(ranges["all"]),
                shells=ranges,
            )
        )
        cursor += len(ranges["all"])
    return OrbitalMap(atoms=atoms, total_orbitals=cursor, stru_path=stru_path, input_path=input_path, orbital_dir=orbital_dir)


def format_orbital_map(orbital_map: OrbitalMap) -> str:
    lines = ["atom_index symbol orbital_start orbital_stop shells"]
    for atom in orbital_map.atoms:
        shells = ",".join(shell for shell in SHELL_ORDER + ["all"] if shell in atom.shells)
        lines.append(f"{atom.atom_index} {atom.symbol} {atom.orbital_start} {atom.orbital_stop} {shells}")
    return "\n".join(lines)


def parse_global_indices(text: str) -> list[int]:
    indices = [int(token.strip()) for token in text.split(",") if token.strip()]
    if not indices:
        raise ValueError("empty orbital index list")
    if any(index < 0 for index in indices):
        raise ValueError("global NAO indices must be non-negative")
    return indices


def selector_shell(token: str) -> str | None:
    if token == "all":
        return "all"
    match = re.fullmatch(r"(?:\d+)?([spdfg])", token)
    return match.group(1) if match else None


def resolve_atom_orbitals(orbital_map: OrbitalMap, atom_index: int, selector: str = "all") -> OrbitalSelection:
    if atom_index < 1 or atom_index > len(orbital_map.atoms):
        raise ValueError(f"atom index {atom_index} is outside 1..{len(orbital_map.atoms)}")
    tokens = [token.strip().lower() for token in (selector or "all").split(",") if token.strip()]
    atom = orbital_map.atoms[atom_index - 1]
    selected: set[int] = set()
    normalized: list[str] = []
    for token in tokens:
        shell = selector_shell(token)
        if shell is None:
            raise ValueError(f"invalid orbital selector: {token}")
        if shell not in atom.shells:
            available = ",".join(shell for shell in SHELL_ORDER + ["all"] if shell in atom.shells)
            raise ValueError(f"atom {atom_index} ({atom.symbol}) has no {shell} shell; available: {available}")
        if shell not in normalized:
            normalized.append(shell)
        selected.update(atom.shells[shell])
    return OrbitalSelection(atom_index, atom.symbol, selector, normalized, sorted(selected))


def infer_companion_file(out_dir: Path, explicit_path: Path | None, filename: str) -> Path | None:
    if explicit_path:
        return Path(explicit_path)
    out_dir = Path(out_dir)
    return _first_existing([out_dir / filename, out_dir.parent / filename])


def resolve_orbital_arguments(
    out_dir: Path,
    atom_i_orbs: str,
    atom_j_orbs: str,
    atom_i_index: int | None = None,
    atom_j_index: int | None = None,
    stru_path: Path | None = None,
    input_path: Path | None = None,
    orbital_dir: Path | None = None,
) -> tuple[list[int], list[int], tuple[OrbitalSelection, OrbitalSelection, OrbitalMap] | None]:
    if atom_i_index is None and atom_j_index is None:
        return parse_global_indices(atom_i_orbs), parse_global_indices(atom_j_orbs), None
    if atom_i_index is None or atom_j_index is None:
        raise ValueError("atom-i-index and atom-j-index must be used together")
    stru_path = infer_companion_file(out_dir, stru_path, "STRU")
    input_path = infer_companion_file(out_dir, input_path, "INPUT")
    if stru_path is None:
        raise FileNotFoundError("cannot find STRU; pass --stru explicitly")
    orbital_map = build_orbital_map(stru_path, input_path=input_path, orbital_dir=orbital_dir)
    sel_i = resolve_atom_orbitals(orbital_map, atom_i_index, atom_i_orbs or "all")
    sel_j = resolve_atom_orbitals(orbital_map, atom_j_index, atom_j_orbs or "all")
    return sel_i.indices, sel_j.indices, (sel_i, sel_j, orbital_map)


def cohp_values_for_kpoint(
    hmat: np.ndarray,
    smat: np.ndarray,
    coeffs: np.ndarray,
    energies_ry: np.ndarray,
    atom_i_orbs: list[int],
    atom_j_orbs: list[int],
    method: str,
) -> tuple[np.ndarray, np.ndarray]:
    if hmat.shape != smat.shape or hmat.shape[0] != hmat.shape[1]:
        raise ValueError("H/S matrix dimensions do not match")
    if coeffs.shape[0] != hmat.shape[0]:
        raise ValueError(f"wavefunction has {coeffs.shape[0]} orbitals but matrix has {hmat.shape[0]}")
    matrix = smat if method == "COOP" else hmat
    block = matrix[np.ix_(atom_i_orbs, atom_j_orbs)]
    values = np.einsum("ib,ij,jb->b", coeffs[atom_i_orbs, :].conj(), block, coeffs[atom_j_orbs, :]).real
    return energies_ry * RY_TO_EV, values


def accumulate_cohp(
    h_mats: list[np.ndarray],
    s_mats: list[np.ndarray],
    coeffs: list[np.ndarray],
    energies: list[np.ndarray],
    weights: np.ndarray,
    atom_i_orbs: list[int],
    atom_j_orbs: list[int],
    method: str,
) -> tuple[np.ndarray, np.ndarray]:
    all_energy: list[float] = []
    all_values: list[float] = []
    for hmat, smat, cmat, evals, weight in zip(h_mats, s_mats, coeffs, energies, weights):
        energy_ev, values = cohp_values_for_kpoint(hmat, smat, cmat, evals, atom_i_orbs, atom_j_orbs, method)
        all_energy.extend(float(x) for x in energy_ev)
        all_values.extend(float(weight) * float(x) for x in values)
    energy = np.asarray(all_energy)
    values = np.asarray(all_values)
    order = np.argsort(energy)
    energy = energy[order]
    values = values[order]
    unique = np.unique(energy)
    summed = np.asarray([np.sum(values[energy == item]) for item in unique], dtype=float)
    return unique, summed


def zero_pad_spectrum(energy: np.ndarray, values: np.ndarray, de: float) -> tuple[np.ndarray, np.ndarray]:
    if de <= 0:
        raise ValueError("de must be positive")
    emin = float(np.min(energy))
    emax = float(np.max(energy))
    span = emax - emin
    padding = max(de, 0.05 * span)
    emin -= padding
    emax += padding
    grid = np.arange(emin, emax + de * 0.5, de)
    out = np.zeros(len(grid), dtype=float)
    for x, y in zip(energy, values):
        if x < grid[0] or x > grid[-1]:
            continue
        idx = int(np.abs(grid - x).argmin())
        if abs(grid[idx] - x) <= de / 2:
            out[idx] += y
    return grid, out


def gaussian_smooth(energy: np.ndarray, values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return values
    diff = energy[:, None] - energy[None, :]
    kernel = np.exp(-(diff * diff) / (2.0 * sigma * sigma))
    return kernel @ values


def integrate_to_fermi(energy: np.ndarray, values: np.ndarray, efermi: float) -> float:
    energy = np.asarray(energy, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(energy)
    energy = energy[order]
    values = values[order]
    mask = energy <= efermi
    if not np.any(mask):
        return 0.0
    integ_energy = energy[mask]
    integ_values = values[mask]
    if integ_energy[-1] < efermi and np.any(energy > efermi):
        right = int(np.argmax(energy > efermi))
        left = right - 1
        if left >= 0:
            value_at_fermi = np.interp(efermi, [energy[left], energy[right]], [values[left], values[right]])
            integ_energy = np.append(integ_energy, efermi)
            integ_values = np.append(integ_values, value_at_fermi)
    if len(integ_energy) < 2:
        return 0.0
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(integ_values, integ_energy))
    widths = np.diff(integ_energy)
    heights = 0.5 * (integ_values[:-1] + integ_values[1:])
    return float(np.sum(widths * heights))


def sorted_matrix_files(out_dir: Path, suffix: str) -> list[Path]:
    def key(path: Path) -> int:
        match = re.search(r"data-(\d+)-[HS]$", path.name)
        return int(match.group(1)) if match else 0

    return sorted(Path(out_dir).glob(f"data-*-{suffix}"), key=key)


def find_wfc_file(out_dir: Path, index: int) -> Path | None:
    return _first_existing(
        [
            out_dir / f"WFC_NAO_K{index}.txt",
            out_dir / f"WFC_NAO_K{index}_ION1.txt",
            out_dir / f"LOWF_K_{index}.txt",
            out_dir / f"LOWF_K_{index}.dat",
            out_dir / f"WFC_NAO_GAMMA{index}.txt",
            out_dir / f"WFC_NAO_GAMMA{index}_ION1.txt",
            out_dir / f"LOWF_GAMMA_S{index}.dat",
        ]
    )


def load_cohp_inputs(out_dir: Path, spin: str = "sum"):
    out_dir = Path(out_dir)
    h_files = sorted_matrix_files(out_dir, "H")
    s_files = sorted_matrix_files(out_dir, "S")
    if not h_files or len(h_files) != len(s_files):
        raise FileNotFoundError(f"cannot find matching data-*-H/data-*-S in {out_dir}")
    h_mats = [read_hs_matrix(path) for path in h_files]
    s_mats = [read_hs_matrix(path) for path in s_files]

    wfc_files = []
    for idx in range(1, len(h_files) + 1):
        wfc = find_wfc_file(out_dir, idx)
        if wfc is None:
            raise FileNotFoundError(f"cannot find WFC_NAO/LOWF file for k index {idx} in {out_dir}")
        wfc_files.append(wfc)
    wfc_data = [read_wfc_nao(path) for path in wfc_files]
    coeffs = [item[0] for item in wfc_data]
    energies = [item[2] for item in wfc_data]

    base_weights = read_kpoint_weights(out_dir / "kpoints", len(h_files))
    if len(h_files) == len(base_weights):
        weights = base_weights
        nspin = 1
    elif len(h_files) == 2 * len(base_weights):
        weights = np.concatenate([base_weights, base_weights])
        nspin = 2
    else:
        raise ValueError(f"cannot map {len(h_files)} matrix files to {len(base_weights)} k-point weights")

    spin = spin.lower()
    if spin not in {"sum", "up", "down"}:
        raise ValueError("spin must be one of: sum, up, down")
    if nspin == 1 and spin in {"up", "down"}:
        raise ValueError(f"requested spin={spin}, but output appears to be nspin=1")
    if nspin == 2 and spin in {"up", "down"}:
        nk = len(base_weights)
        selected = slice(0, nk) if spin == "up" else slice(nk, 2 * nk)
        h_mats = h_mats[selected]
        s_mats = s_mats[selected]
        coeffs = coeffs[selected]
        energies = energies[selected]
        weights = weights[selected]

    return h_mats, s_mats, coeffs, energies, weights, read_fermi_energy(out_dir / "running_scf.log")


def plot_cohp(path: Path, energy: np.ndarray, values: np.ndarray, efermi: float, method: str, invert: bool, emin, emax, width):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    x = -values if invert else values
    y = energy - efermi
    mask = np.ones(len(y), dtype=bool)
    if emin is not None:
        mask &= y >= emin
    if emax is not None:
        mask &= y <= emax
    if not np.any(mask):
        raise ValueError("no COHP data points in requested energy window")
    x = x[mask]
    y = y[mask]
    if width is None:
        width = max(abs(float(np.min(x))), abs(float(np.max(x))), 1.0e-12)

    plt.figure(figsize=(4.5, 7.0))
    plt.plot(x, y, lw=1.0)
    plt.axvline(0.0, color="black", lw=0.6)
    if (emin is None or emin <= 0.0) and (emax is None or emax >= 0.0):
        plt.axhline(0.0, color="black", lw=0.6, ls="--")
    plt.fill_betweenx(y, x, 0.0, where=(y <= 0.0), alpha=0.25)
    plt.xlim(-width, width)
    if emin is not None or emax is not None:
        plt.ylim(emin, emax)
    plt.xlabel(("-" if invert else "") + method)
    plt.ylabel("Energy - E_Fermi (eV)")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()
    return path


def run_cohp(
    out_dir: Path,
    atom_i_orbs: list[int],
    atom_j_orbs: list[int],
    method: str = "COHP",
    spin: str = "sum",
    de: float = 0.05,
    smooth: bool = True,
    smooth_nstddev: float = 4.0,
    invert: bool = True,
    output_prefix: Path | str = "COHP",
    emin: float | None = -10.0,
    emax: float | None = 10.0,
    width: float | None = None,
) -> dict:
    method = method.upper()
    if method not in {"COHP", "COOP"}:
        raise ValueError("built-in COHP currently supports COHP and COOP")
    h_mats, s_mats, coeffs, energies, weights, efermi = load_cohp_inputs(out_dir, spin=spin)
    nlocal = h_mats[0].shape[0]
    requested = atom_i_orbs + atom_j_orbs
    if max(requested) >= nlocal:
        raise ValueError(f"requested NAO index {max(requested)} but output has {nlocal} orbitals")

    energy, values = accumulate_cohp(h_mats, s_mats, coeffs, energies, weights, atom_i_orbs, atom_j_orbs, method)
    energy, values = zero_pad_spectrum(energy, values, de)
    if smooth:
        values = gaussian_smooth(energy, values, smooth_nstddev * de)

    prefix = Path(output_prefix)
    raw_path = prefix.with_suffix(".dat")
    shifted_path = prefix.with_name(prefix.stem + "_EminusEf").with_suffix(".dat")
    icohp_path = prefix.with_suffix(".icohp.txt")
    json_path = prefix.with_suffix(".json")
    png_path = prefix.with_suffix(".png")
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(raw_path, np.column_stack([energy, values]), header=f"Energy(eV) {method}")
    np.savetxt(shifted_path, np.column_stack([energy - efermi, values]), header=f"Energy-E_Fermi(eV) {method}")
    icohp = integrate_to_fermi(energy, values, efermi)
    icohp_path.write_text(
        "\n".join(
            [
                f"method {method}",
                f"efermi_ev {efermi:.12g}",
                f"icohp_raw_eV {icohp:.12g}",
                f"minus_icohp_eV {-icohp:.12g}",
            ]
        )
        + "\n"
    )
    plotted = plot_cohp(png_path, energy, values, efermi, method, invert, emin, emax, width)
    metadata = {
        "method": method,
        "spin": spin,
        "efermi_ev": float(efermi),
        "icohp_raw_ev": float(icohp),
        "minus_icohp_ev": float(-icohp),
        "atom_i_orbs": atom_i_orbs,
        "atom_j_orbs": atom_j_orbs,
        "de_ev": float(de),
        "smooth": bool(smooth),
        "smooth_nstddev": float(smooth_nstddev),
        "invert_plot": bool(invert),
        "files": {
            "raw": str(raw_path),
            "shifted": str(shifted_path),
            "icohp": str(icohp_path),
            "plot": str(plotted) if plotted else None,
            "metadata": str(json_path),
        },
    }
    json_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
