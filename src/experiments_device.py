from qiskit.providers.fake_provider import Fake127QPulseV1
from qiskit import transpile, QuantumCircuit
from topology_functions import *
from functools import reduce
from qiskit.converters import circuit_to_dag
import logging
import csv
import multiprocessing
import os
import sys
import time


class ExperimentsNew:
    def __init__(self,
                 crosstalk_version,
                 crosstalk_fidelity,
                 i,
                 neighbor_fidelity,
                 coupling_map,
                 optimization_level,
                 directory,
                 gate_errors,
                 connectivity_density,
                 depolarization_error,
                 backend,
                 gate_fidelity=None,
                 gate_set=None,
                 ):

        if connectivity_density is None:
            self.connectivity_density = [0.013895, 0.015, 0.018, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.1, 0.15,
                                         0.2, 0.25,
                                         0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        #   self.connectivity_density =  [0.01]
        if gate_set is None:
            self.gate_set = ['id', 'rz', 'sx', 'x', 'cx', 'swap', 'cz']
        if gate_fidelity is None:
            self.gate_fidelity = {
                'cx': 0.99522,
                'cz': 0.9980,
                'swap': 0.99522,
                'x': 0.9997285,
                'rz': 0.9997285,
                'id': 0.9997285,
                'sx': 0.9997285,
                'measure': 0.9997285,
            }
        self.connectivity_density = connectivity_density
        self.depolarization_error = depolarization_error
        self.crosstalk_version = crosstalk_version
        self.directory = directory
        self.gate_errors = gate_errors
        self.crosstalk_fidelity = crosstalk_fidelity
        self.i = i
        self.neighbor_fidelity = neighbor_fidelity
        self.coupling_map = coupling_map
        self.optimization_level = optimization_level
        self.gate_set = gate_set
        self.backend = backend

        self.properties = backend.properties()
        self.gate_set_backend = backend.configuration().basis_gates

    # Function to calculate depolarization and thermal relaxation fidelities
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
            print(qubit)
            if qubit < len(properties.qubits):
                t1 = properties.qubits[qubit][1].value  # T1 relaxation time
                t2 = properties.qubits[qubit][2].value  # T2 dephasing time
                if t1 and t2:
                    print("i have a t1 and t2")
                    print(properties.qubits[qubit][1])
                    print(qubit)
                    print(properties.qubits[qubit][1].value)
                    print(properties.qubits[qubit][2].value)
                    thermal_noise = 1 - (1 / t1 + 1 / (2 * t2))  # Simplified formula
                    fidelity *= thermal_noise
        return fidelity

    def load_and_prepare_circuit_measure(self, quantum_circuit):

        # Check if the circuit has any measurement gates
        if not any(inst[0].name == "measure" for inst in quantum_circuit.data):
            # Add measurements to all qubits if none exist
            quantum_circuit.measure_all()
            print("Measurement gates added to the circuit.")
        else:
            print("Circuit already contains measurement gates.")

        return quantum_circuit

    def remove_measure_gates(self, circuit):
        # Erstelle einen neuen Circuit ohne die Messgates
        new_circuit = QuantumCircuit(circuit.num_qubits, circuit.num_clbits)

        for instr, qargs, cargs in circuit.data:
            if instr.name != 'measure':
                new_circuit.append(instr, qargs, cargs)

        return new_circuit

    def calculate_fidelity_depolarization(self, transpiled_circuit, single_qubit_error, two_qubit_error):
        # Initiale Fidelity ist 1 (perfekter Zustand)
        fidelity = 1.0
        # circuit = self.remove_measure_gates(transpiled_circuit)

        # Iteriere über alle Gates im Circuit
        for instr, qargs, _ in transpiled_circuit.data:
            #  print(instr.name)
            if instr.name in ['u1', 'u2', 'u3', 'h', 'x', 'y', 'z', 'rx', 'ry', 'rz', 'measure']:
                # Ein-Qubit-Gate
                fidelity *= (1 - single_qubit_error)
            elif instr.name in ["cx", "swap", "cz"]:
                # Zwei-Qubit-Gate
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
                    #circuit_fidelity *= (harmonic_value1 * harmonic_value2)
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

    def run_experiment_for_file(self, filename, save_file):
        cur_dir = os.getcwd()
        try:
            if cur_dir not in sys.path:
                sys.path.append(cur_dir)
        except Exception as e:
            print(f"Error occurred: {e}")
        directory = self.directory
        circuit_directory = os.path.join(cur_dir, directory)
        simulator_logs = "logs_" + save_file + "_" + filename + "_" + self.crosstalk_version + "_" + str(
            self.crosstalk_fidelity[self.i])
        file_csv = "csv_" + save_file + "_" + filename[:-5] + "_" + self.crosstalk_version + "_" + str(
            self.crosstalk_fidelity[self.i]) + '.csv'
        logging.basicConfig(filename=simulator_logs, level=logging.INFO)

        with open(file_csv, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                ['connectivity', 'crosstalk_fidelity', 'crosstalk_version', 'optimisation_level',
                 'num_qubits_qc', 'num_qubits_transpiled', 'depth', 'depth_transpiled_circ',
                 'time', 'gate_counts', 'gate_names ', 'gate_counts_transpiled_circ', 'gate_names_transpiled_circ',
                 'fidelity_crosstalk', "fidelity_crosstalk_nn", 'total_harmonic', 'harmonic_values',
                 "fidelity_all_errors",
                 "fidelity_gate_errors", "fidelity_thermal_relaxation", "fidelity_depolarization", "fidelity",
                 "neighbor_fidelity"])

        file_path = os.path.join(circuit_directory, filename)
        quantum_circuit = QuantumCircuit.from_qasm_file(file_path)
        quantum_circuit = self.load_and_prepare_circuit_measure(quantum_circuit)
        # backend = Aer.get_backend('qasm_simulator')

        for connectivity in range(len(self.connectivity_density)):
            cmap_ext = increase_coupling_density(self.coupling_map, self.connectivity_density[connectivity])
            start = time.time()

            transpiled_circ = transpile(quantum_circuit, coupling_map=cmap_ext,
                                        optimization_level=self.optimization_level, basis_gates=self.gate_set
                                        )

            end = time.time()
            times = []
            depth = []
            depth_transpiled = []
            gate_counts = []
            gate_counts_transpiled = []
            gate_names = []
            gate_names_transpiled = []
            crosstalk_fidelities = []
            #crosstalk fidelity without neighbors
            crosstalk_fidelities_nn = []
            # all error without cx noise
            fidelity_all_errors = []
            fidelity_gate_errors = []
            fidelity_thermal_relaxation = []
            fidelity_depolarization = []
            # overall fidelity
            fidelity = []

            if self.crosstalk_version == 'cxneighbors':

                crosstalk_fidelity_calculated = self.calculate_crosstalk_fidelity_cx_neighbors(transpiled_circ,
                                                                                               True,
                                                                                               self.neighbor_fidelity)
                crosstalk_fidelities.append(crosstalk_fidelity_calculated)

                crosstalk_fidelity_nn = self.calculate_crosstalk_fidelity_cx_neighbors(transpiled_circ, False,
                                                                                       self.neighbor_fidelity)
                crosstalk_fidelities_nn.append(crosstalk_fidelity_nn)
                other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                fidelity_all_errors.append(other_errors["fidelity_all_errors"])
                fidelity_gate_errors.append(other_errors["fidelity_gate_errors_only"])
                fidelity_thermal_relaxation.append(other_errors["fidelity_thermal_relaxation_only"])
                fidelity_depolarization.append(other_errors["fidelity_depolarization_only"])
                fidelity.append(other_errors["fidelity_all_errors"] * crosstalk_fidelity_calculated)

                c_values = ["n.a"]
                total_harmonic = "n.a"

            elif self.crosstalk_version == 'topology':
                crosstalk_fidelity_calculated = self.calculate_crosstalk_fidelity_proximity(transpiled_circ,
                                                                                            True,
                                                                                            self.neighbor_fidelity,
                                                                                            max_distance=2)

                crosstalk_fidelities.append(crosstalk_fidelity_calculated)
                crosstalk_fidelity_nn = ["n.a"]
                crosstalk_fidelities_nn.append(crosstalk_fidelity_nn)
                c_values = ["n.a"]
                total_harmonic = "n.a"

                if isinstance(crosstalk_fidelity_calculated, (int, float)):
                    final_metric = crosstalk_fidelity_calculated
                else:
                    final_metric = 1.0  # Default to 1.0 if calculation fails

                other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                fidelity_all_errors.append(other_errors["fidelity_all_errors"])
                fidelity_gate_errors.append(other_errors["fidelity_gate_errors_only"])
                fidelity_thermal_relaxation.append(other_errors["fidelity_thermal_relaxation_only"])
                fidelity_depolarization.append(other_errors["fidelity_depolarization_only"])
                fidelity.append(other_errors["fidelity_all_errors"] * final_metric)



            else:

                try:
                    simultaneous_gates, neighbors, all_two_qubit_gates, all_neighbors, unique_neighbors, fidelities = \
                        self.find_simultaneous_two_qubit_gates_2(transpiled_circ)
                    c_values = []

                    f = self.crosstalk_fidelity[self.i]
                    harmonic_value = self.harmonic_mean(f, f)

                    # Use a list comprehension to create the list of harmonic values
                    c_values = [harmonic_value for _ in all_two_qubit_gates]

                    n = len(fidelities)
                    if n > 0:
                        inverse_sum = sum(1 / v for v in fidelities if v != 0)
                    else:
                        inverse_sum = 0

                    if inverse_sum != 0:
                        c_values.append(n / inverse_sum)

                    total_harmonic = self.calculate_total_harmonic_mean(c_values)
                    n = len(c_values)
                    final_metric = total_harmonic ** n if n > 1 else 'n.a'

                    crosstalk_fidelities.append(final_metric)
                    crosstalk_fidelity_nn = ["n.a"]
                    crosstalk_fidelities_nn.append(crosstalk_fidelity_nn)

                    other_errors = self.calculate_other_errors(transpiled_circ, True, True, True, 1.0, 1.0, 1.0, 1.0)
                    fidelity_all_errors.append(other_errors["fidelity_all_errors"])
                    fidelity_gate_errors.append(other_errors["fidelity_gate_errors_only"])
                    fidelity_thermal_relaxation.append(other_errors["fidelity_thermal_relaxation_only"])
                    fidelity_depolarization.append(other_errors["fidelity_depolarization_only"])
                    fidelity.append(other_errors["fidelity_all_errors"] * final_metric)

                except Exception as e:

                    print(e)

            times.append(end - start)
            depth.append(quantum_circuit.depth())
            depth_transpiled.append(transpiled_circ.depth())
            for gate, count in quantum_circuit.count_ops().items():
                gate_names.append(gate)
                gate_counts.append(count)

            for gate, count in transpiled_circ.count_ops().items():
                gate_names_transpiled.append(gate)
                gate_counts_transpiled.append(count)

            with open(file_csv, 'a', newline='') as file:
                print("we are writing our file")
                writer = csv.writer(file)
                writer.writerow([self.connectivity_density[connectivity], self.crosstalk_fidelity[self.i],
                                 self.crosstalk_version, self.optimization_level,
                                 quantum_circuit.num_qubits, transpiled_circ.num_qubits,
                                 depth, depth_transpiled, times, gate_counts, gate_names,
                                 gate_counts_transpiled, gate_names_transpiled,
                                 crosstalk_fidelities, crosstalk_fidelities_nn, total_harmonic, c_values,
                                 fidelity_all_errors, fidelity_gate_errors, fidelity_thermal_relaxation,
                                 fidelity_depolarization, fidelity, self.neighbor_fidelity
                                 ])

    def run_experiment(self, save_file):
        cur_dir = os.getcwd()
        if cur_dir not in sys.path:
            sys.path.append(cur_dir)
        circuit_directory = os.path.join(cur_dir, self.directory)
        files = os.listdir(circuit_directory)
        qasm_files = [f for f in files if f.endswith('.qasm')]

        # Parallel execution using multiprocessing
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() // 2)
        for filename in qasm_files:
            pool.apply_async(self.run_experiment_for_file, args=(filename, save_file))

        pool.close()
        pool.join()


if __name__ == "__main__":

    backend = Fake127QPulseV1()
    #  backend_sizes = [[6, 4], [6, 5], [8, 4], [6, 6], [8, 5], [8, 6]]
    #  backend_sizes = [[6, 4], [6, 5]]
    #  backend_sizes_syc = [[6, 6], [7, 7], [8, 8], [8, 10], [9, 10], [10, 10], [11, 11], [12, 12]]
    backend_sizes = [[6, 6], [12, 12]]

    # coupling_map_sycamore = create_sycamore_topology(6, 6)
    #coupling_map = create_heavy_hex_IBMQ(6, 5)

    connectivity_density = [0.013895, 0.015, 0.018, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.1, 0.15, 0.2, 0.25,
                            0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    #  connectivity_density = [0.01, 0.02]
    #cmap_ext = increase_coupling_density(coupling_map, connectivity_density[0])
    gate_set = ['id', 'rz', 'sx', 'x', 'cx', 'swap', 'cz']
    # Initialize the base gate fidelities
    gate_errors = {
        'cx': 0.00479,
        'cz': 0.00220,
        'swap': 0.00479,
        'x': 0.0002715,
        'rz': 0.0002715,
        'id': 0.0002715,
        'sx': 0.0002715,
        'measure': 0.0002715,
    }

    gate_fidelity = {
        'cx': 0.99522,
        'cz': 0.9980,
        'swap': 0.99522,
        'x': 0.9997285,
        'rz': 0.9997285,
        'id': 0.9997285,
        'sx': 0.9997285,
        'measure': 0.9997285,
    }

    directory = (r'circuits')

    #  crosstalk_fidelity = [1.0, 0.98847, 0.97694, 0.96541, 0.94361, 0.9218, 0.9, 0.82, 0.74, 0.66]
    crosstalk_fidelity = [0.98847]
    depolarization_error = 0.99661
    gate_error = 0.00253075
    #cxneighbors
    #9997285
    # for backend_size in backend_sizes_syc
    for backend_size in backend_sizes:
        rows, cols = backend_size
        #      coupling_map = create_heavy_hex_IBMQ(rows, cols)
        coupling_map = create_sycamore_topology(rows, cols)

        # for opt_level in [0, 3]:
        for opt_level in [0]:
            for j in range(len(crosstalk_fidelity)):
                if j == 0:
                    neighbor_fidelity = 0.9997285
                else:
                    neighbor_fidelity = crosstalk_fidelity[j] + 0.1 * crosstalk_fidelity[j]
                    if neighbor_fidelity >= 1.0:
                        neighbor_fidelity = 0.9997285

                exp = ExperimentsNew(
                    "topology", crosstalk_fidelity, j, neighbor_fidelity, coupling_map,
                    optimization_level=opt_level, directory=directory, gate_errors=gate_errors,
                    connectivity_density=connectivity_density, depolarization_error=depolarization_error,
                    backend=backend
                )
                filename = f"syc_syc_syc_top_bz_{rows}_{cols}_opt{opt_level}_R1"
                exp.run_experiment(filename)

#    for j in range(len(crosstalk_fidelity)):
#       if j == 0:
#            neighbor_fidelity = 0.9997285
#       else:
#           neighbor_fidelity = crosstalk_fidelity[j] + 0.1 * crosstalk_fidelity[j]
#           if neighbor_fidelity >= 1.0:
#               neighbor_fidelity = 0.9997285

#      exp = ExperimentsNew("cxneighbors", crosstalk_fidelity, j, 0.9997285, coupling_map,
#                        optimization_level=3, directory=directory, gate_errors=gate_errors,
#                        connectivity_density=connectivity_density, depolarization_error=depolarization_error)

#       exp.run_experiment("hhex_6_5_R1_opt3")

# Create a random circuit for testing (replace this with your circuit) --> replace with own circuits
# num_qubits = 100
# depth = 20
# circuit = random_circuit(num_qubits=num_qubits, depth=depth, max_operands=3)
# transpiled_circuit = transpile(circuit, coupling_map=cmap_ext, basis_gates=gate_set, optimization_level=3)
'TODO problem im never in the gate errors depo error and thermal error insdie'
#   fidelity = exp.calculate_crosstalk_fidelity_cx_neighbors(transpiled_circuit, False, 0.999999)
#   results = exp.calculate_other_errors(circuit, True, True, True, 1.0, 1.0, 1.0, 1.0)

# print(fidelity)
# print(results)

#  simultaneous_gates, neighbors, all_two_qubit_gates, all_neighbors, unique_neighbors, fidelities = \
#      exp.find_simultaneous_two_qubit_gates_2(transpiled_circuit)
#  c_values = []
#  for i in all_two_qubit_gates:
#      f = crosstalk_fidelity[0]
#      harmonic_value = exp.harmonic_mean(f, f)
#      c_values.append(harmonic_value)

#   n = len(fidelities)
#   if n > 0:
#       inverse_sum = sum(1 / v for v in fidelities if v != 0)
#   else:
#       inverse_sum = 0#

#   if inverse_sum != 0:
#       c_values.append(n / inverse_sum)


#  total_harmonic = exp.calculate_total_harmonic_mean(c_values)
#  n = len(c_values)
#  final_metric = total_harmonic ** n if n > 1 else 'n.a'

#  crosstalk_fidelity_result.append(final_metric)

#   print(crosstalk_fidelity_result)


#cmap_ext = increase_coupling_density(coupling_map, connectivity_density[20])
# properties = backend.properties()
# gate_set = backend.configuration().basis_gates

# Transpile the circuit for the fake backend
# transpiled_circuit = transpile(circuit, coupling_map=c, basis_gates=gate_set,
#                                optimization_level=3)

#results = exp.calculate_other_errors(circuit, True, True, True, 1.0, 1.0, 1.0, 1.0)

# Print results
#print("Fidelity Results (Original Circuit):", results)

#check if error or fidelity
# complete .csv file
