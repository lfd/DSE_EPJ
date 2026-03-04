import argparse
import csv
import json
import logging
import multiprocessing as mp
import time
import numpy as np
import networkx as nx
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union
from qiskit.converters import circuit_to_dag
from functools import reduce
from qiskit import QuantumCircuit, transpile
from qiskit.providers.fake_provider import Fake127QPulseV1
from src.libs.topology_functions import create_sycamore_topology, create_heavy_hex_IBMQ, increase_coupling_density
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# -------------------------
# Repo-relative paths
# -------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CIRCUITS_DIR = PROJECT_ROOT / "circuits"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "results"


# -------------------------
# Config object: all knobs with defaults
# -------------------------
@dataclass
class ExperimentConfig:
    # Paths
    circuits_dir: Path = DEFAULT_CIRCUITS_DIR
    results_dir: Path = DEFAULT_RESULTS_DIR

    # Experiment sweep knobs
    backend_sizes: List[Tuple[int, int]] = field(default_factory=lambda: [(6,6)]) #[(6, 6), (12, 12)])
    opt_levels: List[int] = field(default_factory=lambda: [0])
    connectivity_density: List[float] = field(default_factory=lambda: [ 0.013895, 0.015])
  #      0.013895, 0.015, 0.018, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
  #      0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
  #  ])

    # Crosstalk
    crosstalk_version: str = "topology"  # "cxneighbors", "topology",
    crosstalk_fidelity: List[float] = field(default_factory=lambda: [0.98847])
    neighbor_fidelity: Optional[float] = 0.9997285  # can be None -> auto rule below
    coupling_map: str = "sycamore"  # "sycamore" or "heavy_hex"

    # Gates and errors
    gate_set: List[str] = field(default_factory=lambda: ['id', 'rz', 'sx', 'x', 'cx', 'swap', 'cz'])

    gate_errors: Dict[str, float] = field(default_factory=lambda: {
        'cx': 0.00479,
        'cz': 0.00220,
        'swap': 0.00479,
        'x': 0.0002715,
        'rz': 0.0002715,
        'id': 0.0002715,
        'sx': 0.0002715,
        'measure': 0.0002715,
    })

    gate_fidelity: Dict[str, float] = field(default_factory=lambda: {
        'cx': 0.99522,
        'cz': 0.9980,
        'swap': 0.99522,
        'x': 0.9997285,
        'rz': 0.9997285,
        'id': 0.9997285,
        'sx': 0.9997285,
        'measure': 0.9997285,
    })

    depolarization_error: float = 0.99661

    # Backend selection
    backend_name: str = "Fake127QPulseV1"  # you can add other backends

    # Parallelism
    processes: Optional[int] = None  # None => cpu_count//2


# -------------------------
# Your experiment class
# -------------------------
class ExperimentsNew:
    def __init__(
        self,
        *,
        config: ExperimentConfig,
        backend,
        coupling_map: Sequence[Tuple[int, int]],
        crosstalk_index: int,
        optimization_level: int,
        rows: int,
        cols: int,
    ):
        self.config = config
        self.backend = backend
        self.coupling_map = coupling_map

        self.i = crosstalk_index
        self.optimization_level = optimization_level
        self.rows = rows
        self.cols = cols

        # Defaults from config
        self.connectivity_density = config.connectivity_density
        self.gate_set = config.gate_set
        self.gate_errors = config.gate_errors
        self.gate_fidelity = config.gate_fidelity
        self.depolarization_error = config.depolarization_error
        self.crosstalk_version = config.crosstalk_version
        self.crosstalk_fidelity = config.crosstalk_fidelity

        # Neighbor fidelity rule: use user-provided or auto rule
        if config.neighbor_fidelity is not None:
            self.neighbor_fidelity = config.neighbor_fidelity
        else:
            # Your old rule: j==0 => 0.9997285 else 1.1*crosstalk_fidelity capped
            if self.i == 0:
                self.neighbor_fidelity = 0.9997285
            else:
                nf = self.crosstalk_fidelity[self.i] * 1.1
                self.neighbor_fidelity = 0.9997285 if nf >= 1.0 else nf

        self.properties = backend.properties()
        self.gate_set_backend = backend.configuration().basis_gates

        # Ensure output folder exists
        self.config.results_dir.mkdir(parents=True, exist_ok=True)

    def estimate_thermal_fidelity(self, circuit, t1s, t2s, gate_times):
        """
        Estimate fidelity loss from thermal relaxation.

        Parameters:
        - circuit: Qiskit QuantumCircuit object
        - t1s: dict mapping qubit index to T1 time (in ns)
        - t2s: dict mapping qubit index to T2 time (in ns)
        - gate_times: dict mapping gate name to duration in ns (e.g., {'id': 50, 'cx': 200, 'rz': 0, ...})

        Returns:
        - estimated fidelity (float between 0 and 1)
        """
        from collections import defaultdict
        qubit_times = defaultdict(float)

        # Track time each qubit is active
        for instruction, qargs, _ in circuit.data:
            gate_name = instruction.name
            duration = gate_times.get(gate_name, 0)
            for q in qargs:
                qubit_times[q.index] += duration

        total_fidelity = 1.0

        for q, total_time in qubit_times.items():
            t1 = t1s[q]
            t2 = t2s[q]
            t_phi = 1 / (1 / t2 - 1 / (2 * t1)) if t2 != 0 else float('inf')

            f_t1 = np.exp(-total_time / t1)
            f_tphi = np.exp(-total_time / t_phi)
            f_q = f_t1 * f_tphi

            total_fidelity *= f_q  # multiply fidelities across qubits

        return total_fidelity

    def calculate_thermal_relaxation_fidelity(self, qubits, properties):
        fidelity = 1.0
        for qubit in qubits:
            if qubit < len(properties.qubits):
                t1 = properties.qubits[qubit][1].value  # T1 relaxation time
                t2 = properties.qubits[qubit][2].value  # T2 dephasing time
                if t1 and t2:
                    thermal_noise = 1 - (1 / t1 + 1 / (2 * t2))
                    fidelity *= thermal_noise
        return fidelity

    def load_and_prepare_circuit_measure(self, quantum_circuit):

        # Check if the circuit has any measurement gates
        if not any(inst[0].name == "measure" for inst in quantum_circuit.data):
            # Add measurements to all qubits if none exist
            quantum_circuit.measure_all()
        return quantum_circuit

    def remove_measure_gates(self, circuit):
        new_circuit = QuantumCircuit(circuit.num_qubits, circuit.num_clbits)

        for instr, qargs, cargs in circuit.data:
            if instr.name != 'measure':
                new_circuit.append(instr, qargs, cargs)

        return new_circuit

    def calculate_fidelity_depolarization(self, transpiled_circuit, single_qubit_error, two_qubit_error):
        fidelity = 1.0
        # circuit = self.remove_measure_gates(transpiled_circuit)
        for instr, qargs, _ in transpiled_circuit.data:
            if instr.name in ['u1', 'u2', 'u3', 'h', 'x', 'y', 'z', 'rx', 'ry', 'rz', 'measure']:
                fidelity *= (1 - single_qubit_error)
            elif instr.name in ["cx", "swap", "cz"]:
                fidelity *= (1 - two_qubit_error)

        return fidelity

    def calculate_gate_fidelity(self, gate_error):
        overall_depolarization_fidelity = gate_error

        return overall_depolarization_fidelity

    def calculate_depolarization_fidelity(self, qubits, properties):
        fidelity = 1.0
        for qubit in qubits:
            if qubit < len(properties.qubits):
                # Iterate through the gates to find depolarizing errors for single-qubit and two-qubit gates
                for gate in properties.gates:
                    if gate.qubits == [qubit]:  # Single-qubit gate
                        for param in gate.parameters:
                            if param.name == 'gate_error':
                                gate_error = param.value
                                depolarizing_noise = 1 - gate_error
                                fidelity *= depolarizing_noise

                    elif qubit in gate.qubits and len(gate.qubits) == 2:  # Two-qubit gate
                        for param in gate.parameters:
                            if param.name == 'gate_error':
                                gate_error = param.value
                                depolarizing_noise = 1 - gate_error
                                fidelity *= depolarizing_noise
        return fidelity

    def calculate_depolarization_fidelity_(self, qubits, properties):
        import math

        log_fidelity = 0.0  # Using log space for stability

        for gate in properties.gates:
            gate_qubits = gate.qubits
            num_qubits_in_gate = len(gate_qubits)

            # Check if the gate operates on the selected qubits
            if any(qubit in qubits for qubit in gate_qubits):
                for param in gate.parameters:
                    if param.name == 'gate_error':
                        gate_error = param.value

                        # Skip invalid gate errors
                        if gate_error is None or gate_error >= 1.0 or gate_error < 0.0:
                            continue

                        # Distinguish between single-qubit and two-qubit gates
                        if num_qubits_in_gate == 1:
                            # Single-qubit depolarizing model
                            depolarizing_noise = 1 - gate_error
                        elif num_qubits_in_gate == 2:
                            # Two-qubit depolarizing noise affects fidelity more
                            depolarizing_noise = 1 - gate_error
                        else:
                            # Ignore gates with more than 2 qubits (rare)
                            continue

                        # Prevent log(0)
                        depolarizing_noise = max(depolarizing_noise, 1e-10)

                        # Accumulate in log space
                        log_fidelity += math.log(depolarizing_noise)

        # Convert back from log space
        fidelity = math.exp(log_fidelity)

        # Prevent fidelity from being exactly zero
        return max(fidelity, 1e-10)

        # include_gate, include_thermal, include_depolarization : Bool
        # fidelity_all, fidelity_gate, fidelity_thermal, fidelity_depo : real, 0:1.0

    def calculate_other_errors(self, circuit, include_gate, include_thermal, include_depolarization_error, fidelity_all,
                               fidelity_gate, fidelity_thermal, fidelity_depo):
        include_gate_errors = include_gate
        include_thermal_relaxation = include_thermal
        include_depolarization = include_depolarization_error

        fidelity_all_errors = fidelity_all
        fidelity_gate_errors_only = fidelity_gate
        fidelity_thermal_relaxation_only = fidelity_thermal
        fidelity_depolarization_only = fidelity_depo

        for instruction in circuit.data:
            try:
                operation = instruction[0]  # Access the Instruction object

                # Extract the qubit indices correctly
                qubits = [q._index for q in instruction[1]]  # Use `_index` for Qubit indices

            except Exception as e:
                print(f"Error processing instruction {instruction}: {e}")
                continue

            # Thermal Relaxation Noise
            if include_thermal_relaxation:
                thermal_relaxation_fidelity = self.calculate_thermal_relaxation_fidelity(qubits, self.properties)
                fidelity_all_errors *= thermal_relaxation_fidelity
                fidelity_thermal_relaxation_only *= thermal_relaxation_fidelity

                # Depolarization Noise
            if include_depolarization:
                #         result = self.calculate_fidelity_depolarization(circuit, 0.0005, 0.005)
                result = self.calculate_depolarization_fidelity(qubits, self.properties)
                fidelity_depolarization_only = result
                fidelity_all_errors *= result

            if include_gate_errors:
                fidelity_gate_errors_only = self.calculate_gate_fidelity(self.depolarization_error)

        results = {
            "fidelity_all_errors": fidelity_all_errors,
            "fidelity_gate_errors_only": fidelity_gate_errors_only,
            "fidelity_thermal_relaxation_only": fidelity_thermal_relaxation_only,
            "fidelity_depolarization_only": fidelity_depolarization_only
        }

        return results

    def save_error(self, file_csv, error_message, key, backend):
        # Create the filename
        file_base = f'{file_csv}{key.strip(".json")}_{backend}_{self.optimization_level}'
        error_filename = f'{file_base}_error.txt'

        # Save error message in a text file
        with open(error_filename, 'w') as errorfile:
            errorfile.write(error_message)

    def harmonic_mean(self, f1, f2):
        return 2 / ((1 / f1) + (1 / f2))

    def multiply_list_elements(self, elements):
        # Return the product of all elements in the list
        return reduce(lambda x, y: x * y, elements)

    def get_fidelity(self, crosstalk_fidelities, i):
        return crosstalk_fidelities[i]

    def calculate_total_harmonic_mean(self, harmonic_values):
        n = len(harmonic_values)
        if n == 0:
            return None  # Or return 0 if you prefer
        inverse_sum = sum(1 / hv for hv in harmonic_values)
        total_harmonic_mean = n / inverse_sum
        return total_harmonic_mean

    def find_simultaneous_two_qubit_gates_2_(self, circuit):
        dag = circuit_to_dag(circuit)
        simultaneous_two_qubit_gates = []
        neighbors = []
        all_simultaneous_two_qubit_gates = []  # Track only simultaneous two-qubit gates
        all_neighbors = []

        for layer in dag.layers():
            two_qubit_gates_in_layer = []
            layer_qubits = set()

            for op_node in layer['graph'].op_nodes():
                if len(op_node.qargs) == 2:  # Check if the gate is a two-qubit gate
                    gate_str = f"{op_node.name} on qubits {op_node.qargs[0].index} and {op_node.qargs[1].index}"
                    two_qubit_gates_in_layer.append(gate_str)
                    layer_qubits.update(op_node.qargs)

            if len(two_qubit_gates_in_layer) > 1:  # Only consider layers with simultaneous two-qubit gates
                simultaneous_two_qubit_gates.append(two_qubit_gates_in_layer)
                all_simultaneous_two_qubit_gates.extend(
                    [op_node.name for op_node in layer['graph'].op_nodes() if len(op_node.qargs) == 2])
                neighbors_in_layer = self.find_neighbors_2(layer_qubits, circuit.qregs)
                neighbors.append(neighbors_in_layer)
                all_neighbors.extend(q.index for q in layer_qubits)

        unique_neighbors = list(set(all_neighbors))
        fidelities = [self.get_fidelity(circuit, qubit) for qubit in unique_neighbors]

        return simultaneous_two_qubit_gates, neighbors, all_simultaneous_two_qubit_gates, all_neighbors, unique_neighbors, fidelities

    def find_simultaneous_two_qubit_gates_2(self, circuit):
        dag = circuit_to_dag(circuit)
        simultaneous_two_qubit_gates = []
        neighbors = []
        all_simultaneous_two_qubit_gates = []
        all_neighbors = []

        for layer in dag.layers():
            two_qubit_gates_in_layer = []
            layer_qubits = set()

            for op_node in layer["graph"].op_nodes():
                if len(op_node.qargs) == 2:
                    qubit_indices = [circuit.find_bit(q).index for q in op_node.qargs]
                    gate_str = f"{op_node.name} on qubits {qubit_indices[0]} and {qubit_indices[1]}"
                    two_qubit_gates_in_layer.append(gate_str)
                    layer_qubits.update(qubit_indices)

            if len(two_qubit_gates_in_layer) > 1:
                simultaneous_two_qubit_gates.append(two_qubit_gates_in_layer)
                all_simultaneous_two_qubit_gates.extend(
                    [op_node.name for op_node in layer["graph"].op_nodes() if len(op_node.qargs) == 2])

                # Fetch neighbors from the coupling map
                for q_index in layer_qubits:
                    all_neighbors.extend(self.find_neighbors_2(q_index, self.coupling_map))

        unique_neighbors = list(set(all_neighbors))
        fidelities = [self.get_fidelity(self.crosstalk_fidelity, self.i) for q in unique_neighbors]

        return simultaneous_two_qubit_gates, neighbors, all_simultaneous_two_qubit_gates, all_neighbors, unique_neighbors, fidelities

    def find_neighbors_2_(self, layer_qubits, qregs):
        neighbors = []
        all_qubits = {q for qreg in qregs for q in qreg}

        for qubit in layer_qubits:
            neighboring_qubits = [q.index for q in all_qubits if q != qubit and (q in layer_qubits)]
            neighbors.append(f"Qubit {qubit.index} neighbors: {neighboring_qubits}")

        return neighbors

    def find_neighbors_2(self, qubit_index, coupling_map):
        neighbors = []
        for pair in coupling_map:
            if qubit_index in pair:
                neighbors.append(pair[1] if pair[0] == qubit_index else pair[0])
        return neighbors

    def calculate_crosstalk_fidelity_cx_neighbors_(self, transpiled_circ, apply_neighbor_noise, neighbor_fidelity):
        qubit_graph = nx.Graph()

        # Build the qubit graph from the circuit
        for gate in transpiled_circ.data:
            if len(gate.qubits) == 2:
                # Convert qubits to a consistent string representation
                q1 = str(gate.qubits[0])
                q2 = str(gate.qubits[1])
                qubit_graph.add_edge(q1, q2, gate=gate)

        circuit_fidelity = 1.0
        processed_edges = set()  # Keep track of processed edge pairs

        for edge1 in qubit_graph.edges:
            for edge2 in qubit_graph.edges:
                # Sort edges for consistent representation
                edge1_sorted = tuple(sorted(edge1))
                edge2_sorted = tuple(sorted(edge2))

                # Skip the same edge or processed edge pairs
                if edge1_sorted == edge2_sorted or (edge1_sorted, edge2_sorted) in processed_edges:
                    continue

                # Check if edges share a common qubit
                shared_qubit = set(edge1).intersection(edge2)
                if shared_qubit:
                    # Retrieve fidelities for the qubits in the edges
                    node1, node2 = edge1
                    node3, node4 = edge2

                    f1 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                    f2 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                    f3 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                    f4 = self.get_fidelity(self.crosstalk_fidelity, self.i)

                    # Calculate harmonic means for both edges
                    harmonic_value1 = self.harmonic_mean(f1, f2)
                    harmonic_value2 = self.harmonic_mean(f3, f4)

                    # Update fidelity for this pair of connected gates
                    # circuit_fidelity *= (harmonic_value1 * harmonic_value2)
                    circuit_fidelity = circuit_fidelity * (harmonic_value1 * harmonic_value2)

                    # Adjust for neighbors of the current edge
                    if apply_neighbor_noise:
                        neighbors = (set(qubit_graph.neighbors(node1)) | set(qubit_graph.neighbors(node2)) |
                                     set(qubit_graph.neighbors(node3)) | set(qubit_graph.neighbors(node4)))
                        for neighbor in neighbors:
                            # Exclude nodes in the current edge pair
                            if neighbor not in edge1 and neighbor not in edge2:
                                circuit_fidelity *= neighbor_fidelity

                    # Mark edge pair as processed
                    processed_edges.add((edge1_sorted, edge2_sorted))

        return circuit_fidelity

    def calculate_crosstalk_fidelity_cx_neighbors__(self, transpiled_circ, apply_neighbor_noise, neighbor_fidelity):
        qubit_graph = nx.Graph()
        two_qubit_gates = set()  # Track valid two-qubit operations

        # Build the full qubit graph from the circuit (including single- and two-qubit gates)
        for gate in transpiled_circ.data:
            qubits = [str(q) for q in gate.qubits]
            if len(qubits) == 2:
                qubit_graph.add_edge(qubits[0], qubits[1], gate=gate)
                two_qubit_gates.add((qubits[0], qubits[1]))  # Store as two-qubit gate pair
            elif len(qubits) == 1:
                qubit_graph.add_node(qubits[0])  # Include single-qubit operations

        circuit_fidelity = 1.0
        processed_edges = set()  # Track processed edge pairs

        for edge1 in qubit_graph.edges:
            for edge2 in qubit_graph.edges:
                # Sort edges for consistent representation
                edge1_sorted = tuple(sorted(edge1))
                edge2_sorted = tuple(sorted(edge2))

                # Skip the same edge or already processed edge pairs
                if edge1_sorted == edge2_sorted or (edge1_sorted, edge2_sorted) in processed_edges:
                    continue

                # Ensure both edges represent two-qubit gates
                if edge1_sorted in two_qubit_gates and edge2_sorted in two_qubit_gates:
                    # Check if the two two-qubit gates share a common qubit (i.e., they are neighbors)
                    shared_qubit = set(edge1).intersection(edge2)
                    if shared_qubit:
                        node1, node2 = edge1
                        node3, node4 = edge2

                        # Get fidelity values for each gate
                        f1 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f2 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f3 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f4 = self.get_fidelity(self.crosstalk_fidelity, self.i)

                        # Calculate harmonic means for both edges
                        harmonic_value1 = self.harmonic_mean(f1, f2)
                        harmonic_value2 = self.harmonic_mean(f3, f4)

                        # Apply crosstalk noise only if the two two-qubit operations are neighbors
                        circuit_fidelity *= (harmonic_value1 * harmonic_value2)

                        # Apply neighbor noise only if they are neighbors
                        if apply_neighbor_noise:
                            # Collect all involved qubits
                            involved_qubits = set(edge1) | set(edge2)

                            for qubit in involved_qubits:
                                for neighbor in qubit_graph.neighbors(qubit):
                                    # Apply noise only to direct neighboring qubits not part of the two-qubit operations
                                    if neighbor not in involved_qubits:
                                        circuit_fidelity *= neighbor_fidelity

                        # Mark edge pair as processed
                        processed_edges.add((edge1_sorted, edge2_sorted))

        return circuit_fidelity

    def calculate_crosstalk_fidelity_cx_neighbors(self, transpiled_circ, apply_neighbor_noise, neighbor_fidelity):
        qubit_graph = nx.Graph()
        two_qubit_gates = set()  # Track valid two-qubit operations

        # Build the full qubit graph from the circuit (including single- and two-qubit gates)
        for gate in transpiled_circ.data:
            qubits = [str(q) for q in gate.qubits]
            if len(qubits) == 2:
                qubit_graph.add_edge(qubits[0], qubits[1], gate=gate)
                two_qubit_gates.add(tuple(sorted([qubits[0], qubits[1]])))  # Store as sorted pairs
            elif len(qubits) == 1:
                qubit_graph.add_node(qubits[0])  # Include single-qubit operations

        circuit_fidelity = 1.0
        processed_edges = set()  # Track processed edge pairs

        for edge1 in qubit_graph.edges:
            for edge2 in qubit_graph.edges:
                # Sort edges for consistent representation
                edge1_sorted = tuple(sorted(edge1))
                edge2_sorted = tuple(sorted(edge2))

                # Skip the same edge or already processed edge pairs
                if edge1_sorted == edge2_sorted or (edge1_sorted, edge2_sorted) in processed_edges:
                    continue

                # Ensure both edges represent two-qubit gates
                if edge1_sorted in two_qubit_gates and edge2_sorted in two_qubit_gates:
                    # Check if the two two-qubit gates share a common qubit (i.e., they are neighbors)
                    shared_qubit = set(edge1).intersection(edge2)
                    if shared_qubit:
                        # Get fidelity values for each gate
                        f1 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f2 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f3 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                        f4 = self.get_fidelity(self.crosstalk_fidelity, self.i)

                        # Calculate harmonic means for both edges
                        harmonic_value1 = self.harmonic_mean(f1, f2)
                        harmonic_value2 = self.harmonic_mean(f3, f4)

                        # Apply crosstalk noise only if the two two-qubit operations are neighbors
                        circuit_fidelity *= (harmonic_value1 * harmonic_value2)

                        # Apply neighbor noise only if they are neighbors
                        if apply_neighbor_noise:
                            involved_qubits = set(edge1) | set(edge2)

                            for qubit in involved_qubits:
                                for neighbor in qubit_graph.neighbors(qubit):
                                    # Apply noise only to direct neighboring qubits not part of the two-qubit operations
                                    if neighbor not in involved_qubits:
                                        circuit_fidelity *= neighbor_fidelity

                        # Mark edge pair as processed
                        processed_edges.add((edge1_sorted, edge2_sorted))

        return circuit_fidelity

    def calculate_crosstalk_fidelity_proximity(self, transpiled_circ, apply_neighbor_noise, neighbor_fidelity,
                                               max_distance=2):
        """
        Calculates the fidelity of a circuit with crosstalk noise occurring based on physical qubit proximity.

        Args:
            transpiled_circ (QuantumCircuit): Transpiled quantum circuit.
            apply_neighbor_noise (bool): Whether to apply noise to neighboring single-qubit gates.
            neighbor_fidelity (float): Fidelity penalty for single-qubit neighboring operations.
            max_distance (int): Maximum physical distance between two CX gates to consider crosstalk.

        Returns:
            float: Circuit fidelity after applying crosstalk and neighbor noise.
        """
        # fake_backend = Fake127QPulseV1()
        qubit_coordinates = {i: (q[0].value, q[1].value) for i, q in enumerate(self.backend.properties().qubits)}
        cx_gates = []
        for gate in transpiled_circ.data:
            qubits = [transpiled_circ.find_bit(q).index for q in gate.qubits]
            if len(qubits) == 2 and gate.operation.name == "cx":
                cx_gates.append(tuple(sorted(qubits)))

        circuit_fidelity = 1.0
        processed_pairs = set()

        for i, edge1 in enumerate(cx_gates):
            for j, edge2 in enumerate(cx_gates):
                if i >= j:
                    continue

                distances = [
                    np.linalg.norm(np.array(qubit_coordinates[edge1[k]]) - np.array(qubit_coordinates[edge2[m]]))
                    for k in range(2) for m in range(2)
                ]

                if min(distances) <= max_distance:
                    f1 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                    f2 = self.get_fidelity(self.crosstalk_fidelity, self.i)
                    harmonic_value = self.harmonic_mean(f1, f2)
                    circuit_fidelity *= harmonic_value

                    if apply_neighbor_noise:
                        involved_qubits = set(edge1) | set(edge2)
                        for qubit in involved_qubits:
                            for neighbor, coord in qubit_coordinates.items():
                                if neighbor not in involved_qubits:
                                    neighbor_dist = np.linalg.norm(np.array(qubit_coordinates[qubit]) - np.array(coord))
                                    if neighbor_dist <= max_distance:
                                        circuit_fidelity *= neighbor_fidelity

                    processed_pairs.add((edge1, edge2))

        return circuit_fidelity

    def run_experiment_for_file(self, filename: str, run_tag: str):
        circuit_path = self.config.circuits_dir / filename
        if not circuit_path.exists():
            raise FileNotFoundError(f"Missing circuit file: {circuit_path}")

        # Stable filenames (portable)
        log_path = self.config.results_dir / (
            f"logs_{run_tag}_{filename}_{self.crosstalk_version}_{self.crosstalk_fidelity[self.i]}_opt{self.optimization_level}.log"
        )
        csv_path = self.config.results_dir / (
            f"csv_{run_tag}_{filename[:-5]}_{self.crosstalk_version}_{self.crosstalk_fidelity[self.i]}_opt{self.optimization_level}.csv"
        )

        logging.basicConfig(filename=str(log_path), level=logging.INFO)

        # Write header once
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                'connectivity', 'crosstalk_fidelity', 'crosstalk_version', 'optimisation_level',
                'num_qubits_qc', 'num_qubits_transpiled', 'depth', 'depth_transpiled_circ',
                'time', 'gate_counts', 'gate_names', 'gate_counts_transpiled_circ', 'gate_names_transpiled_circ',
                'fidelity_crosstalk', 'fidelity_crosstalk_nn', 'total_harmonic', 'harmonic_values',
                'fidelity_all_errors', 'fidelity_gate_errors', 'fidelity_thermal_relaxation',
                'fidelity_depolarization', 'fidelity', 'neighbor_fidelity'
            ])

        quantum_circuit = QuantumCircuit.from_qasm_file(str(circuit_path))
        quantum_circuit = self.load_and_prepare_circuit_measure(quantum_circuit)

        # Sweep over connectivity
        for conn in self.connectivity_density:
            cmap_ext = increase_coupling_density(self.coupling_map, conn)

            start = time.time()
            transpiled_circ = transpile(
                quantum_circuit,
                coupling_map=cmap_ext,
                optimization_level=self.optimization_level,
                basis_gates=self.gate_set
            )
            elapsed = time.time() - start

            # IMPORTANT: I’m leaving your metric computation placeholders here.
            # Plug in your existing computations exactly as you do now:
            depth_qc = quantum_circuit.depth()
            depth_tr = transpiled_circ.depth()

            gate_names = list(quantum_circuit.count_ops().keys())
            gate_counts = list(quantum_circuit.count_ops().values())
            gate_names_tr = list(transpiled_circ.count_ops().keys())
            gate_counts_tr = list(transpiled_circ.count_ops().values())

            # --- compute crosstalk + other errors (ported from your old code) ---
            if self.crosstalk_version == "cxneighbors":
                crosstalk_fidelity_calculated = self.calculate_crosstalk_fidelity_cx_neighbors(
                    transpiled_circ, True, self.neighbor_fidelity
                )
                fidelity_crosstalk = crosstalk_fidelity_calculated

                crosstalk_fidelity_nn = self.calculate_crosstalk_fidelity_cx_neighbors(
                    transpiled_circ, False, self.neighbor_fidelity
                )
                fidelity_crosstalk_nn = crosstalk_fidelity_nn

                other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                fidelity_all_errors = other_errors["fidelity_all_errors"]
                fidelity_gate_errors = other_errors["fidelity_gate_errors_only"]
                fidelity_thermal = other_errors["fidelity_thermal_relaxation_only"]
                fidelity_depo = other_errors["fidelity_depolarization_only"]

                fidelity_total = fidelity_all_errors * fidelity_crosstalk

                harmonic_values = ["n.a"]
                total_harmonic = "n.a"

            elif self.crosstalk_version == "topology":
                crosstalk_fidelity_calculated = self.calculate_crosstalk_fidelity_proximity(
                    transpiled_circ, True, self.neighbor_fidelity, max_distance=2
                )
                fidelity_crosstalk = crosstalk_fidelity_calculated
                fidelity_crosstalk_nn = "n.a"
                harmonic_values = ["n.a"]
                total_harmonic = "n.a"

                # Guard: if your method ever returns non-float
                final_metric = fidelity_crosstalk if isinstance(fidelity_crosstalk, (int, float)) else 1.0

                other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                fidelity_all_errors = other_errors["fidelity_all_errors"]
                fidelity_gate_errors = other_errors["fidelity_gate_errors_only"]
                fidelity_thermal = other_errors["fidelity_thermal_relaxation_only"]
                fidelity_depo = other_errors["fidelity_depolarization_only"]

                fidelity_total = fidelity_all_errors * final_metric

            else:
                # Your "simultaneous two-qubit" harmonic metric branch
                try:
                    simultaneous_gates, neighbors, all_two_qubit_gates, all_neighbors, unique_neighbors, fidelities = \
                        self.find_simultaneous_two_qubit_gates_2(transpiled_circ)

                    f = self.crosstalk_fidelity[self.i]
                    harmonic_value = self.harmonic_mean(f, f)

                    harmonic_values = [harmonic_value for _ in all_two_qubit_gates]

                    n = len(fidelities)
                    inverse_sum = sum(1 / v for v in fidelities if v != 0) if n > 0 else 0

                    if inverse_sum != 0:
                        harmonic_values.append(n / inverse_sum)

                    total_harmonic = self.calculate_total_harmonic_mean(harmonic_values)

                    m = len(harmonic_values)
                    final_metric = (total_harmonic ** m) if (m > 1 and total_harmonic is not None) else "n.a"

                    fidelity_crosstalk = final_metric
                    fidelity_crosstalk_nn = "n.a"

                    other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                    fidelity_all_errors = other_errors["fidelity_all_errors"]
                    fidelity_gate_errors = other_errors["fidelity_gate_errors_only"]
                    fidelity_thermal = other_errors["fidelity_thermal_relaxation_only"]
                    fidelity_depo = other_errors["fidelity_depolarization_only"]

                    fidelity_total = fidelity_all_errors * (
                        final_metric if isinstance(final_metric, (int, float)) else 1.0)

                except Exception as e:
                    # If anything goes wrong, fail gracefully but log it
                    logging.exception(f"Error in crosstalk else-branch: {e}")
                    fidelity_crosstalk = "n.a"
                    fidelity_crosstalk_nn = "n.a"
                    harmonic_values = ["n.a"]
                    total_harmonic = "n.a"

                    fidelity_all_errors = 1.0
                    fidelity_gate_errors = 1.0
                    fidelity_thermal = 1.0
                    fidelity_depo = 1.0
                    fidelity_total = 1.0

            with open(csv_path, "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow([
                    conn, self.crosstalk_fidelity[self.i], self.crosstalk_version, self.optimization_level,
                    quantum_circuit.num_qubits, transpiled_circ.num_qubits,
                    depth_qc, depth_tr, elapsed,
                    gate_counts, gate_names,
                    gate_counts_tr, gate_names_tr,
                    fidelity_crosstalk, fidelity_crosstalk_nn, total_harmonic, harmonic_values,
                    fidelity_all_errors, fidelity_gate_errors, fidelity_thermal, fidelity_depo,
                    fidelity_total, self.neighbor_fidelity
                ])

    def run_experiment(self, run_tag: str):
        circuit_dir = self.config.circuits_dir
        qasm_files = sorted([p.name for p in circuit_dir.iterdir() if p.is_file() and p.suffix == ".qasm"])

        if not qasm_files:
            raise RuntimeError(f"No .qasm files found in {circuit_dir}")

        procs = self.config.processes
        if procs is None:
            procs = max(1, mp.cpu_count() // 2)

        with mp.Pool(processes=procs) as pool:
            for fn in qasm_files:
                pool.apply_async(self.run_experiment_for_file, args=(fn, run_tag))
            pool.close()
            pool.join()


# -------------------------
# Backend factory (extend later for real backends)
# -------------------------
def make_backend(name: str):
    if name == "Fake127QPulseV1":
        return Fake127QPulseV1()
    raise ValueError(f"Unknown backend: {name}")

# -------------------------
# Create coupling map
# -------------------------
def make_coupling_map(coupling_map: str, rows: int, cols: int):
    if coupling_map == "sycamore":
        return create_sycamore_topology(rows, cols)
    if coupling_map == "heavy_hex":
        return create_heavy_hex_IBMQ(rows, cols)
    raise ValueError(f"Unknown topology: {coupling_map}")



# -------------------------
# CLI parsing
# -------------------------
def parse_backend_sizes(values: List[str]) -> List[Tuple[int, int]]:
    out = []
    for v in values:
        if "x" not in v:
            raise ValueError(f"Invalid backend size '{v}'. Use like 6x6 or 12x12.")
        a, b = v.split("x", 1)
        out.append((int(a), int(b)))
    return out


def parse_args() -> ExperimentConfig:
    p = argparse.ArgumentParser(description="Reproducible Qiskit crosstalk experiments runner")

    p.add_argument("--circuits-dir", type=str, default=str(DEFAULT_CIRCUITS_DIR))
    p.add_argument("--results-dir", type=str, default=str(DEFAULT_RESULTS_DIR))

    p.add_argument("--backend", type=str, default="Fake127QPulseV1")
    p.add_argument("--backend-sizes", nargs="*", default=["6x6", "12x12"])
    p.add_argument("--opt-levels", nargs="*", type=int, default=[0])

    p.add_argument("--crosstalk-version", type=str, default="topology")
    p.add_argument("--crosstalk-fidelity", nargs="*", type=float, default=[0.98847])
    p.add_argument("--neighbor-fidelity", type=float, default=0.9997285)

    p.add_argument("--connectivity-density", nargs="*", type=float, default=None)
    p.add_argument("--gate-set", nargs="*", default=None)

    # Optional JSON dict overrides (easy for users)
    p.add_argument("--gate-errors-json", type=str, default=None,
                   help="JSON string or path to JSON file with gate errors dict")
    p.add_argument("--gate-fidelity-json", type=str, default=None,
                   help="JSON string or path to JSON file with gate fidelity dict")

    p.add_argument("--depolarization-error", type=float, default=0.99661)
    p.add_argument("--processes", type=int, default=None)
    p.add_argument("--coupling-map", type=str, choices=["sycamore", "heavy_hex"],
        default="sycamore", help="Coupling map used to build the coupling map."
    )


    args = p.parse_args()

    cfg = ExperimentConfig()
    cfg.circuits_dir = Path(args.circuits_dir).resolve()
    cfg.results_dir = Path(args.results_dir).resolve()

    cfg.backend_name = args.backend
    cfg.backend_sizes = parse_backend_sizes(args.backend_sizes)
    cfg.opt_levels = args.opt_levels

    cfg.crosstalk_version = args.crosstalk_version
    cfg.crosstalk_fidelity = args.crosstalk_fidelity
    cfg.neighbor_fidelity = args.neighbor_fidelity
    cfg.coupling_map = args.coupling_map


    if args.connectivity_density is not None:
        cfg.connectivity_density = args.connectivity_density

    if args.gate_set is not None:
        cfg.gate_set = args.gate_set

    cfg.depolarization_error = args.depolarization_error
    cfg.processes = args.processes

    def _load_json_override(val: Optional[str]) -> Optional[dict]:
        if not val:
            return None
        candidate = Path(val)
        if candidate.exists():
            return json.loads(candidate.read_text())
        return json.loads(val)

    ge = _load_json_override(args.gate_errors_json)
    if ge is not None:
        cfg.gate_errors = ge

    gf = _load_json_override(args.gate_fidelity_json)
    if gf is not None:
        cfg.gate_fidelity = gf

    return cfg


def main():
    cfg = parse_args()

    backend = make_backend(cfg.backend_name)

    # Loop over size / opt / crosstalk index combos
    for (rows, cols) in cfg.backend_sizes:
        coupling_map = make_coupling_map(cfg.coupling_map, rows, cols)


        for opt in cfg.opt_levels:
            for j in range(len(cfg.crosstalk_fidelity)):
                exp = ExperimentsNew(
                    config=cfg,
                    backend=backend,
                    coupling_map=coupling_map,
                    crosstalk_index=j,
                    optimization_level=opt,
                    rows=rows,
                    cols=cols,
                )
                run_tag = f"{cfg.coupling_map}_{rows}x{cols}_opt{opt}_R1"
                exp.run_experiment(run_tag)


if __name__ == "__main__":
    main()
