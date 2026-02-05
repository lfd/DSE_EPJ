#!/usr/bin/env python3
"""
Reproducible benchmark transpilation runner (CLI).

"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from multiprocessing import Pool, cpu_count

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import (
    OptimizeCliffords,
    Optimize1qGates,
    Optimize1qGatesSimpleCommutation,
    Optimize1qGatesDecomposition,
    RemoveDiagonalGatesBeforeMeasure,
    CommutativeCancellation,
    CommutativeInverseCancellation,
    CXCancellation,
    HoareOptimizer,
)

# Your topology util (must provide increase_coupling_density)
import src.topology_functions as tf

# Optional: for building a simple "starting backend" coupling map in the same shape as your old code.
try:
    from qiskit.providers.fake_provider import ConfigurableFakeBackend
except Exception:
    ConfigurableFakeBackend = None


# -------------------------
# PassManager setup presets (ported from your code)
# -------------------------
def setup1(pm: PassManager) -> None:
    pm.append(Optimize1qGates())
    pm.append(OptimizeCliffords())
    pm.append(Optimize1qGatesSimpleCommutation())


def setup2(pm: PassManager) -> None:
    pm.append(Optimize1qGatesDecomposition())
    pm.append(CXCancellation())


def setup3(pm: PassManager) -> None:
    setup1(pm)
    pm.append(RemoveDiagonalGatesBeforeMeasure())


def setup4(pm: PassManager) -> None:
    pm.append(Optimize1qGatesDecomposition())
    pm.append(CommutativeCancellation())


def setup5(pm: PassManager, benchmark_name: str) -> None:
    """
    Your original setup5 had a special-case list for which HoareOptimizer is skipped.
    I preserved that behavior.
    """
    setup3(pm)
    skip_hoare = {
        "20QBT_100CYC_QSE_8.qasm",
        "cycle10_2_110.qasm",
        "plus63mod4096_163.qasm",
        "q=6_s=19994_2qbf=05_1.qasm",
        "q=8_s=39992_2qbf=08_1.qasm",
        "shor_15.qasm",
        "shor_35.qasm",
        "square_root_7.qasm",
    }
    if benchmark_name not in skip_hoare:
        try:
            pm.append(HoareOptimizer())
        except Exception:
            # Keep going even if HoareOptimizer errors out
            pass
    pm.append(CommutativeCancellation())
    pm.append(CommutativeInverseCancellation())


SETUP_FUNCS: Dict[int, Callable[..., None]] = {
    1: setup1,
    2: setup2,
    3: setup3,
    4: setup4,
    # setup5 needs benchmark_name, so it's handled specially
}


# -------------------------
# Config
# -------------------------
@dataclass
class Config:
    circuits_dir: Path
    out_path: Path

    processes: int = max(1, cpu_count() // 2)

    # Device set
    include_bristlecone: bool = True
    include_custom_density_maps: bool = True
    custom_density_values: List[float] = field(default_factory=lambda: [0.013895, 0.03, 0.05, 0.1, 0.3, 0.5, 0.8])
    custom_density_n_qubits: int = 128

    # Transpiler knobs
    optimization_levels: List[int] = field(default_factory=lambda: [0, 1, 2])
    routing_methods: List[str] = field(default_factory=lambda: ["stochastic", "sabre"])
    layout_methods: List[str] = field(default_factory=lambda: ["trivial", "dense", "sabre"])
    setups: List[int] = field(default_factory=lambda: [0, 1, 2, 3, 4, 5])  # 0 = no extra setup PM

    # basis gates you used
    basis_gates: List[str] = field(default_factory=lambda: ["x", "y", "z", "rx", "ry", "rz", "cx", "cy"])


# -------------------------
# Devices / coupling maps
# Each entry: (device_name, coupling_map, num_phys_qubits)
# -------------------------
def bristlecone_coupling_map_72() -> List[List[int]]:
    # This is exactly your list (kept as-is).
    return [
        [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11],
        [12, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20], [20, 21],
        [21, 22], [22, 23], [24, 25], [25, 26], [26, 27], [27, 28], [28, 29], [29, 30], [31, 32],
        [32, 33], [33, 34], [34, 35], [36, 37], [37, 38], [38, 39], [39, 40], [40, 41], [42, 43],
        [43, 44], [44, 45], [45, 46], [46, 47], [48, 49], [49, 50], [50, 51], [51, 52], [52, 53],
        [53, 54], [54, 55], [55, 56], [56, 57], [57, 58], [58, 59], [60, 61], [61, 62], [62, 63],
        [63, 64], [64, 65], [65, 66], [66, 67], [67, 68], [68, 69], [69, 70], [70, 71],
        [1, 12], [1, 14], [14, 3], [3, 16], [16, 5], [5, 18], [18, 7], [7, 20], [9, 20], [9, 22],
        [22, 11], [13, 24], [13, 26], [15, 26], [15, 28], [28, 17], [30, 17], [30, 19], [19, 32],
        [32, 21], [21, 34], [34, 23], [25, 36], [25, 38], [38, 27], [40, 27], [40, 29], [42, 29],
        [42, 31], [31, 44], [44, 33], [33, 46], [46, 35], [37, 48], [37, 50], [50, 39], [39, 52],
        [52, 41], [41, 54], [54, 43], [43, 56], [45, 56], [45, 58], [58, 47], [49, 60], [49, 62],
        [62, 51], [51, 64], [64, 53], [53, 66], [66, 55], [55, 68], [68, 57], [57, 70], [70, 59],

        # reverse edges
        [1, 0], [2, 1], [3, 2], [4, 3], [5, 4], [6, 5], [7, 6], [8, 7], [9, 8], [10, 9], [11, 10],
        [13, 12], [14, 13], [15, 14], [16, 15], [17, 16], [18, 17], [19, 18], [20, 19], [21, 20],
        [22, 21], [23, 22], [25, 24], [26, 25], [27, 26], [28, 27], [29, 28], [30, 29], [32, 31],
        [33, 32], [34, 33], [35, 34], [37, 36], [38, 37], [39, 38], [40, 39], [41, 40], [43, 42],
        [44, 43], [45, 44], [46, 45], [47, 46], [49, 48], [50, 49], [51, 50], [52, 51], [53, 52],
        [54, 53], [55, 54], [56, 55], [57, 56], [58, 57], [59, 58], [61, 60], [62, 61], [63, 62],
        [64, 63], [65, 64], [66, 65], [67, 66], [68, 67], [69, 68], [70, 69], [71, 70],
        [12, 1], [14, 1], [3, 14], [16, 3], [5, 16], [18, 5], [7, 18], [20, 7], [20, 9], [22, 9],
        [11, 22], [24, 13], [26, 13], [26, 15], [28, 15], [17, 28], [17, 30], [19, 30], [32, 19],
        [21, 32], [34, 21], [23, 34], [36, 25], [38, 25], [27, 38], [27, 40], [29, 40], [29, 42],
        [31, 42], [44, 31], [33, 44], [46, 33], [35, 46], [48, 37], [50, 37], [39, 50], [52, 39],
        [41, 52], [54, 41], [43, 54], [56, 43], [56, 45], [58, 45], [47, 58], [60, 49], [62, 49],
        [51, 62], [64, 51], [53, 64], [66, 53], [55, 66], [68, 55], [57, 68], [70, 57], [59, 70],
    ]


def line_coupling_map(n: int) -> List[List[int]]:
    """Bidirectional line coupling map of length n (like your ConfigurableFakeBackend example)."""
    edges = []
    for i in range(n - 1):
        edges.append([i, i + 1])
    for i in range(n - 1):
        edges.append([i + 1, i])
    return edges


def build_devices(cfg: Config) -> List[Tuple[str, List[List[int]], int]]:
    devices: List[Tuple[str, List[List[int]], int]] = []

    if cfg.include_bristlecone:
        devices.append(("Google Bristlecone", bristlecone_coupling_map_72(), 72))

    if cfg.include_custom_density_maps:
        base = line_coupling_map(cfg.custom_density_n_qubits)

        # Use your tf.increase_coupling_density which expects a coupling map + density
        for d in cfg.custom_density_values:
            devices.append((f"Custom_{d}", tf.increase_coupling_density(base, d), cfg.custom_density_n_qubits))

    return devices


# -------------------------
# Core worker
# -------------------------
def process_one_benchmark(args) -> List[List[object]]:
    """
    Worker function for multiprocessing.
    Returns rows for that benchmark.
    """
    qasm_path, bench_name, cfg, devices = args

    # Load
    circuit = QuantumCircuit.from_qasm_file(str(qasm_path))
    depth_before = circuit.depth()
    gates_before = circuit.size()
    gate_types_before = dict(circuit.count_ops())
    twoq_before = gate_types_before.get("cz", 0) + gate_types_before.get("cx", 0)

    rows: List[List[object]] = []

    for (dev_name, coupling_map, num_phys_qubits) in devices:
        if circuit.num_qubits > num_phys_qubits:
            continue

        for ol in cfg.optimization_levels:
            for rt in cfg.routing_methods:
                for lm in cfg.layout_methods:
                    for setup_id in cfg.setups:
                        pm = PassManager()

                        # apply setup variant
                        if setup_id == 0:
                            pass
                        elif setup_id == 5:
                            setup5(pm, bench_name)
                        else:
                            SETUP_FUNCS[setup_id](pm)

                        # In your original code you did pm.run(circuit) but didn't use the returned circuit.
                        # Keep it for reproducibility (it may raise on some circuits).
                        try:
                            _ = pm.run(circuit)
                        except Exception:
                            # If the pass manager fails, skip this setup for this benchmark/device combo
                            continue

                        try:
                            tr = transpile(
                                circuit,
                                basis_gates=cfg.basis_gates,
                                coupling_map=coupling_map,
                                optimization_level=ol,
                                layout_method=lm,
                                routing_method=rt,
                            )
                        except Exception:
                            continue

                        swaps = tr.count_ops().get("swap", 0)
                        depth_after = tr.depth()
                        gates_after = tr.size()
                        gate_types_after = dict(tr.count_ops())

                        twoq_after = (
                            gate_types_after.get("cz", 0)
                            + gate_types_after.get("cx", 0)
                            + 3 * swaps
                        )

                        rows.append([
                            bench_name.replace(".qasm", ""),
                            dev_name,
                            ol,
                            rt,
                            lm,
                            setup_id,
                            gates_before,
                            gates_after,
                            gate_types_before,
                            gate_types_after,
                            twoq_before,
                            twoq_after,
                            swaps,
                            depth_before,
                            depth_after,
                        ])

    return rows


# -------------------------
# CLI
# -------------------------
def parse_list_of_ints(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_list_of_floats(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_list_of_strs(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def parse_args() -> Config:
    p = argparse.ArgumentParser(description="Reproducible benchmark transpilation DSE runner")

    p.add_argument("--circuits-dir", type=str, required=True, help="Directory containing .qasm benchmarks")
    p.add_argument("--out", type=str, required=True, help="Output path (.xlsx or .csv)")

    p.add_argument("--processes", type=int, default=max(1, cpu_count() // 2))

    # Devices
    p.add_argument("--include-bristlecone", action="store_true", default=True)
    p.add_argument("--no-bristlecone", action="store_false", dest="include_bristlecone")

    p.add_argument("--include-custom-density-maps", action="store_true", default=True)
    p.add_argument("--no-custom-density-maps", action="store_false", dest="include_custom_density_maps")

    p.add_argument("--custom-density-values", type=str, default="0.013895,0.03,0.05,0.1,0.3,0.5,0.8")
    p.add_argument("--custom-density-n-qubits", type=int, default=128)

    # Sweep knobs
    p.add_argument("--optimization-levels", type=str, default="0,1,2")
    p.add_argument("--routing-methods", type=str, default="stochastic,sabre")
    p.add_argument("--layout-methods", type=str, default="trivial,dense,sabre")
    p.add_argument("--setups", type=str, default="0,1,2,3,4,5", help="Comma list; 0 means no extra passmanager setup")

    # Basis gates
    p.add_argument("--basis-gates", type=str, default="x,y,z,rx,ry,rz,cx,cy")

    a = p.parse_args()

    cfg = Config(
        circuits_dir=Path(a.circuits_dir).resolve(),
        out_path=Path(a.out).resolve(),
        processes=a.processes,
        include_bristlecone=a.include_bristlecone,
        include_custom_density_maps=a.include_custom_density_maps,
        custom_density_values=parse_list_of_floats(a.custom_density_values),
        custom_density_n_qubits=a.custom_density_n_qubits,
        optimization_levels=parse_list_of_ints(a.optimization_levels),
        routing_methods=parse_list_of_strs(a.routing_methods),
        layout_methods=parse_list_of_strs(a.layout_methods),
        setups=parse_list_of_ints(a.setups),
        basis_gates=parse_list_of_strs(a.basis_gates),
    )

    if not cfg.circuits_dir.exists():
        raise FileNotFoundError(f"circuits-dir not found: {cfg.circuits_dir}")

    return cfg


def main() -> None:
    cfg = parse_args()

    qasm_files = sorted([p for p in cfg.circuits_dir.iterdir() if p.is_file() and p.suffix == ".qasm"])
    if not qasm_files:
        raise RuntimeError(f"No .qasm files in {cfg.circuits_dir}")

    devices = build_devices(cfg)

    work = [(p, p.name, cfg, devices) for p in qasm_files]

    all_rows: List[List[object]] = []
    with Pool(processes=cfg.processes) as pool:
        for rows in pool.imap_unordered(process_one_benchmark, work):
            all_rows.extend(rows)

    df = pd.DataFrame(
        all_rows,
        columns=[
            "benchmark_name",
            "device_name",
            "optimization_level",
            "routing_method",
            "layout_method",
            "setup_number",
            "gates_before",
            "gates_after",
            "gate_types_before",
            "gate_types_after",
            "twoqgates_before",
            "twoqgates_after",
            "swaps",
            "depth_before",
            "depth_after",
        ],
    )

    cfg.out_path.parent.mkdir(parents=True, exist_ok=True)

    if cfg.out_path.suffix.lower() == ".xlsx":
        df.to_excel(cfg.out_path, index=False)
    elif cfg.out_path.suffix.lower() == ".csv":
        df.to_csv(cfg.out_path, index=False)
    else:
        raise ValueError("--out must end with .xlsx or .csv")

    print(f"Wrote {len(df)} rows to {cfg.out_path}")


if __name__ == "__main__":
    main()
