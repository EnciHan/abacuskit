#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/home/hec/apps/abacuskit")

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from abacuskit.cli import (
    find_abacus_input,
    find_abacus_outdir,
    find_running_log,
    iter_job_dirs,
    natural_key,
    parse_abacus_metrics,
)


ABACUSKIT = Path("/home/hec/bin/abacuskit")
DEFAULT_ROOT = Path("/home/hec/data/abacus/cu2o_surface")


st.set_page_config(page_title="abacuskit visual", page_icon="A", layout="wide")


def is_job_dir(path: Path) -> bool:
    return bool(find_abacus_input(path) or find_abacus_outdir(path) or find_running_log(path))


def scan_jobs(root: Path, recursive: bool, max_depth: int) -> list[Path]:
    root = root.expanduser()
    if not root.exists():
        return []
    if not recursive:
        try:
            return iter_job_dirs([root])
        except SystemExit:
            return [root] if is_job_dir(root) else []

    jobs: list[Path] = []
    base_depth = len(root.resolve().parts)
    for current, dirs, _ in os.walk(root):
        current_path = Path(current)
        depth = len(current_path.resolve().parts) - base_depth
        if depth > max_depth:
            dirs[:] = []
            continue
        skip = {".git", "__pycache__", "bader_work"}
        dirs[:] = [item for item in dirs if item not in skip and not item.startswith(".")]
        if is_job_dir(current_path):
            jobs.append(current_path)
            dirs[:] = []
    return sorted(jobs, key=lambda p: natural_key(str(p.relative_to(root) if p != root else p.name)))


def safe_metrics(job: Path) -> dict:
    try:
        return parse_abacus_metrics(job)
    except Exception as exc:
        return {
            "job": str(job),
            "calculation": "",
            "outdir": "",
            "log": "",
            "finished": False,
            "converged": False,
            "failed": True,
            "energy_ev": None,
            "energy": None,
            "message": f"parse failed: {exc}",
        }


@st.cache_data(show_spinner=False)
def load_metrics(root_text: str, recursive: bool, max_depth: int) -> tuple[list[str], list[dict]]:
    root = Path(root_text).expanduser()
    jobs = scan_jobs(root, recursive, max_depth)
    return [str(job) for job in jobs], [safe_metrics(job) for job in jobs]


def display_df(rows: list[dict], root: Path) -> pd.DataFrame:
    records = []
    for row in rows:
        record = {
            "job": short_path(Path(row.get("job", "")), root),
            "type": row.get("calculation"),
            "finished": row.get("finished"),
            "converged": row.get("converged") or row.get("converge"),
            "failed": row.get("failed"),
            "energy_eV": row.get("energy_ev") if row.get("energy_ev") is not None else row.get("energy"),
            "ecutwfc": row.get("ecutwfc"),
            "kpt": row.get("kpt"),
            "nspin": row.get("nspin"),
            "efermi": row.get("efermi"),
            "time_s": row.get("total_time"),
            "scf_steps": row.get("scf_steps"),
            "message": row.get("message"),
        }
        records.append(record)
    return pd.DataFrame(records)


def short_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def read_tail(path: Path, lines: int = 120) -> str:
    if not path.is_file():
        return ""
    data = path.read_text(errors="ignore").splitlines()
    return "\n".join(data[-lines:])


def parse_force_series(log: Path) -> pd.DataFrame:
    if not log.is_file():
        return pd.DataFrame(columns=["step", "force_eV_A"])
    text = log.read_text(errors="ignore")
    values = [float(x) for x in re.findall(r"Largest force is\s+([-+0-9.eE]+)", text)]
    return pd.DataFrame({"step": range(1, len(values) + 1), "force_eV_A": values})


def parse_scf_series(log: Path) -> pd.DataFrame:
    if not log.is_file():
        return pd.DataFrame(columns=["iter", "energy_eV", "ediff_eV", "drho"])
    rows = []
    pattern = re.compile(
        r"^\s*(?:[A-Za-z]+)?\s*([0-9]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)",
        re.MULTILINE,
    )
    for match in pattern.finditer(log.read_text(errors="ignore")):
        rows.append(
            {
                "iter": int(match.group(1)),
                "energy_eV": float(match.group(2)),
                "ediff_eV": float(match.group(3)),
                "drho": float(match.group(4)),
            }
        )
    return pd.DataFrame(rows)


def find_files(root: Path, patterns: list[str], max_items: int = 250) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        files.extend(root.rglob(pattern))
    result = sorted({p for p in files if p.is_file()}, key=lambda p: natural_key(str(p)))
    return result[:max_items]


def bader_summary(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = list(csv.DictReader(path.open()))
    df = pd.DataFrame(rows)
    for col in ["valence_electrons", "bader_electrons", "charge", "x", "y", "z", "min_dist", "atomic_volume"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    summary = (
        df.groupby("symbol", as_index=False)
        .agg(
            atoms=("atom_index", "count"),
            valence_electrons=("valence_electrons", "sum"),
            bader_electrons=("bader_electrons", "sum"),
            net_charge=("charge", "sum"),
        )
        .sort_values("symbol")
    )
    return df, summary


def run_plot_grid(
    path: Path,
    kind: str,
    axis: str,
    index: int | None,
    out: Path,
    file_path: Path | None = None,
    minus_file: Path | None = None,
) -> subprocess.CompletedProcess:
    cmd = [str(ABACUSKIT), "plot-grid", str(path), "--kind", kind, "--axis", axis, "--out", str(out)]
    if index is not None:
        cmd.extend(["--index", str(index)])
    if file_path:
        cmd.extend(["--file", str(file_path)])
    if minus_file:
        cmd.extend(["--minus-file", str(minus_file)])
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def metric_cards(rows: list[dict]) -> None:
    total = len(rows)
    converged = sum(1 for row in rows if row.get("converged") or row.get("converge"))
    failed = sum(1 for row in rows if row.get("failed"))
    running = sum(1 for row in rows if not row.get("finished"))
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Jobs", total)
    c2.metric("Converged", converged)
    c3.metric("Failed", failed)
    c4.metric("Not finished", running)


def download_json_button(rows: list[dict]) -> None:
    data = json.dumps(rows, indent=2, ensure_ascii=False)
    st.download_button("Download metrics JSON", data, file_name="abacuskit_metrics.json", mime="application/json")


def main() -> None:
    st.title("abacuskit visual")
    st.caption("Local ABACUS result browser powered by abacuskit.")

    with st.sidebar:
        st.header("Workspace")
        root_text = st.text_input("Root directory", str(DEFAULT_ROOT))
        recursive = st.checkbox("Recursive scan", value=False)
        max_depth = st.slider("Max depth", 1, 6, 3, disabled=not recursive)
        if st.button("Refresh", use_container_width=True):
            load_metrics.clear()

    root = Path(root_text).expanduser()
    job_paths, rows = load_metrics(str(root), recursive, max_depth)

    if not root.exists():
        st.error(f"Path does not exist: {root}")
        return
    if not rows:
        st.warning("No ABACUS jobs were found under this path.")
        return

    metric_cards(rows)
    overview, details, bader, images, grid, files = st.tabs(
        ["Overview", "Job Details", "Bader", "Grid Images", "Plot Cube", "Files"]
    )

    with overview:
        df = display_df(rows, root)
        st.dataframe(df, use_container_width=True, hide_index=True)
        download_json_button(rows)

        energy_df = df.dropna(subset=["energy_eV"]).copy()
        if not energy_df.empty:
            st.subheader("Energy")
            st.line_chart(energy_df.set_index("job")["energy_eV"])

    with details:
        selected = st.selectbox("Job", job_paths, format_func=lambda p: short_path(Path(p), root))
        job = Path(selected)
        row = rows[job_paths.index(selected)]
        st.json({k: v for k, v in row.items() if k != "INPUT"}, expanded=False)

        log = Path(row.get("log") or "")
        c1, c2 = st.columns(2)
        with c1:
            forces = parse_force_series(log)
            st.subheader("Force")
            if forces.empty:
                st.info("No largest-force series found.")
            else:
                st.line_chart(forces.set_index("step")["force_eV_A"])
                st.dataframe(forces.tail(12), use_container_width=True, hide_index=True)
        with c2:
            scf = parse_scf_series(log)
            st.subheader("SCF residual")
            if scf.empty:
                st.info("No SCF table found.")
            else:
                st.line_chart(scf.set_index("iter")["drho"])
                st.dataframe(scf.tail(12), use_container_width=True, hide_index=True)

        st.subheader("Log tail")
        st.code(read_tail(log), language="text")

    with bader:
        bader_files = find_files(root, ["bader.csv", "*bader*.csv"])
        if not bader_files:
            st.info("No Bader CSV files found.")
        else:
            chosen = st.selectbox("Bader CSV", bader_files, format_func=lambda p: short_path(p, root))
            atoms, summary = bader_summary(chosen)
            st.subheader("Element Summary")
            st.dataframe(summary, use_container_width=True, hide_index=True)
            st.subheader("Atoms")
            st.dataframe(atoms, use_container_width=True, hide_index=True)

    with images:
        image_files = find_files(root, ["*.png", "*.jpg", "*.jpeg"])
        if not image_files:
            st.info("No image files found.")
        else:
            chosen = st.selectbox("Image", image_files, format_func=lambda p: short_path(p, root))
            st.image(str(chosen), caption=str(chosen), use_container_width=True)

    with grid:
        cube_files = find_files(root, ["*.cube"], max_items=120)
        if not cube_files:
            st.info("No cube files found.")
        else:
            selected_cube = st.selectbox("Input cube/job", cube_files, format_func=lambda p: short_path(p, root))
            kind = st.selectbox("Kind", ["auto", "charge", "elf", "cube"])
            axis = st.selectbox("Axis", ["z", "x", "y"])
            index_text = st.text_input("Slice index, empty for middle", "")
            out = st.text_input("Output PNG", str(root / "abacuskit_visual_slice.png"))
            if st.button("Generate plot", type="primary"):
                try:
                    index = int(index_text) if index_text.strip() else None
                    result = run_plot_grid(selected_cube, kind, axis, index, Path(out))
                    if result.returncode != 0:
                        st.error(result.stderr or result.stdout)
                    else:
                        st.success(result.stdout.strip() or f"Wrote {out}")
                        if Path(out).is_file():
                            st.image(str(out), caption=out, use_container_width=True)
                except Exception as exc:
                    st.error(str(exc))

    with files:
        st.subheader("Common outputs")
        common = find_files(
            root,
            ["INPUT", "STRU", "KPT", "abacus.log", "running_*.log", "*.json", "*.csv", "*.cube", "*.png"],
            max_items=500,
        )
        file_df = pd.DataFrame(
            [
                {
                    "path": short_path(p, root),
                    "size_MB": round(p.stat().st_size / 1024 / 1024, 3),
                    "mtime": pd.to_datetime(p.stat().st_mtime, unit="s"),
                }
                for p in common
            ]
        )
        st.dataframe(file_df, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
