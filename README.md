# Design Space Exploration for Quantum Circuit Compilation

This repository contains the code and experimental setup for the paper:

> **Hardware–Software Co-Design in Quantum Circuit Compilation:  
> A Multi-Layer Design Space Exploration**  
> _[Hila Safi, Medina Bandic, Christoph Niedermeier, Carmen G. Almudever, Sebastian Feld, Wolfgang Mauerer 2025]_

We investigate how compilation strategies and hardware characteristics
interact to affect quantum circuit performance, with a focus on **layout methods**,
**qubit routing techniques**, **optimisation levels**, **topological connectivity**, and
**device-specific noise variants** (including crosstalk).  
Our experiments span both **noisy simulations** and **quantum error correction (QEC)** scenarios.


---

## Requirements

We maintain **two separate requirement files** to keep Python dependencies reproducible for each type of experiment:

- `requirements.txt` — for **device experiments**
- `requirements_mapper.txt` — for **mapper/compiler experiments**

Both files are **fully pinned** to ensure long-term reproducibility.

---

## Running with Docker (Recommended)

Using Docker ensures the environment is identical across machines.

### 1. Build the image: Device Experiments

From the repository root:

```bash
docker build -t dse-device -f src/dockerfile . --no-cache

docker run -it dse-device

```

# Build the image: mapper experiments
```bash

docker build -t dse-mapper -f src/dockerfile_mapper . --no-cache

docker run -it dse-mapper

```
---

## Running locally (Without Docker)

# Activate your environment - Device Experiments
```bash
pip install -r requirements.txt
python device_experiments/run_device_experiment.py
```
# Activate your environment - Device Experiments

```bash
pip install -r requirements_mapper.txt - Compilation Experiments
python mapper_experiments/run_mapper_experiment.py
```

## Repository Structure

```text
repo-root/
├── circuits/                  # Input quantum circuits (.qasm)
│   ├── circuit1.qasm
│   ├── circuit2.qasm
│   └── ...
├── src/
│   └── topology_functions.py  # Topology + coupling helpers
├── results/                   # Auto-generated outputs (created if missing)
│   ├── *.csv                  # Experiment results
│   └── *.log                  # Logs
├── run_device_experiments.py  # Main experiment runner
├── README.md
└── requirements.txt
```



## Further instruction on Device Experiments

# What run_device_experiments.py does
1. Loads all .qasm files from circuits/
2. Sweeps over parameters
3. Transpiles each circuit for each configuration
4. Computes Fidelities
5. Writes one CSV per circuit per configuration  to results/

# Running with Defaults (Quickstart) - from repository root
```bash
python -m src.run_device_experiments

```
This uses defaults defined in the code:
- **Backend:** Fake127QPulseV1  
- **Backend sizes:** 6x6, 12x12  
- **Optimisation Level:** 0  
- **Crosstalk version:** topology  
- **Connectivity density:** 0.013895, 0.015, 0.018, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
- **Coupling Map:** Sycamore  
- **Circuits directory:** circuits/  
- **Results directory:** results/


# Command-Line Configurations
Example: Change Backend Size and Optimisation Level

```bash
python -m src.run_device_experiments \
  --backend-sizes 6x6 8x8 \
  --opt-levels 0 3
```

Example: Custom Connectivity Density Sweep

```bash
python -m src.run_device_experiments \
  --connectivity-density 0.01 0.02 0.05 0.1
```

Example: Overwrite Gate Errors (JSON)
```bash
python -m src.run_device_experiments \
  --gate-errors-json '{"cx": 0.005, "x": 0.0003}'
```

# Output Files
We generate per experiment configuration a .csv file which includes all metrics and fidelity results
and a dedicated .log file.

# Parallel Execution
Experiments run in parallel using Python multiprocessing. The default behavior uses
CPU_COUNT // 2 processes.

Override manually with:
```bash
python -m src.run_device_experiments --processes 4
```


## Further instruction on Mapper Experiments

# What DSE_mapper.py does
1. Loads all .qasm files from circuits/
2. Sweeps over parameters
3. Transpiles each circuit for each configuration
4. Writes all results into a single CSV or Excel file

# Compiler Design Parameters

The design space explored corresponds to:

- **Optimisation Level:** `0, 1, 2`
- **Layout Method:** `sabre`, `dense`, `trivial`
- **Routing Technique:** `sabre`, `stochastic`
- **Additional Optimisation Setups:** `0–5`  
  (`0` = no extra setup, `1–5` = custom PassManager presets)
- **Scheduling Method:** `ALAP` (Qiskit default)


# Running with Defaults (Quickstart) - from repository root
```bash
python run_mapper_experiments.py \
  --circuits-dir circuits/ \
  --out results/mapper_results.csv
```

# Device Topology Options
1. include-bristlecone / --no-bristlecone
2. include-custom-density-maps / --no-custom-density-maps
3. custom-density-values 0.013895,0.03,0.05,0.1,0.3,0.5,0.8
4. custom-density-n-qubits 128

# Compiler Design Space Options
1. optimization-levels 0,1,2
2. routing-methods stochastic,sabre
3. layout-methods trivial,dense,sabre
4. setups 0,1,2,3,4,5
5. basis-gates x,y,z,rx,ry,rz,cx,cy





