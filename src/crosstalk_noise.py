from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Gate
from qiskit.transpiler import CouplingMap
from topology_functions import *
from qiskit_aer.noise import NoiseModel, depolarizing_error, thermal_relaxation_error


def neighboring_noise_model_for_active_qubits_id(base_error_prob, neighbor_error_prob, coupling_map, active_qubits):
    noise_model = NoiseModel()
    added_errors = set()

    for connection in coupling_map:
        q1, q2 = connection

        # Depolarizing error for qubits directly involved applied to 'id' gate
        direct_error = depolarizing_error(base_error_prob, 1)
        if ('id', q1) not in added_errors:
            noise_model.add_quantum_error(direct_error, ['id'], [q1])
            added_errors.add(('id', q1))
        if ('id', q2) not in added_errors:
            noise_model.add_quantum_error(direct_error, ['id'], [q2])
            added_errors.add(('id', q2))

        # Handling neighbors
        for qubit in [q1, q2]:
            neighbors = [target for source, target in coupling_map if source == qubit] + \
                        [source for source, target in coupling_map if target == qubit]
            # Filter to include only active qubits
            neighbors = [neighbor for neighbor in set(neighbors) - {q1, q2} if neighbor in active_qubits]

            # Apply depolarizing error to each active neighboring qubit
            for neighbor in neighbors:
                neighbor_error = depolarizing_error(neighbor_error_prob, 1)
                for gate in ['x', 'h', 'rz', 'sx']:
                    if (gate, neighbor) not in added_errors:
                        noise_model.add_quantum_error(neighbor_error, [gate], [neighbor])
                        added_errors.add((gate, neighbor))

    return noise_model


def create_circ_with_add_ids_when2qubit(transpiled_circ, n_qubits):
    qc = QuantumCircuit(n_qubits)
    id_inserted = set()  # To keep track of where id gates are inserted due to two-qubit operations

    for gate, qubits, _ in transpiled_circ.data:
        indices = [qubit._index for qubit in qubits]

        # Apply the gate
        if gate.name in ['cx', 'cz']:  # Consider two-qubit gates that are relevant
            qc.append(gate, qubits)
            for q in indices:
                # Add id gates to each qubit involved in the two-qubit gate
                if q not in id_inserted:
                    qc.id(q)
                    id_inserted.add(q)

        elif gate.name in ['x', 'h', 'rz']:  # Consider single-qubit gates
            qc.append(gate, qubits)

    qc.measure_all()
    return qc


def create_circ_with_add_ids_when22qubit(transpiled_circ, n_qubits):
    qc = QuantumCircuit(n_qubits)
    id_inserted = set()  # To keep track of where id gates are inserted due to two-qubit operations
    last_gate_was_two_qubit = False

    for gate, qubits, _ in transpiled_circ.data:
        indices = [qubit._index for qubit in qubits]

        if gate.name in ['cx', 'cz']:
            if last_gate_was_two_qubit:  # If the last gate was also a two-qubit gate, add id gates
                for q in indices:
                    if q not in id_inserted:
                        qc.id(q)
                        id_inserted.add(q)
            # Apply the current two-qubit gate
            qc.append(gate, qubits)
            last_gate_was_two_qubit = True
        else:
            last_gate_was_two_qubit = False
            # Apply single-qubit gates directly
            if gate.name in ['x', 'h', 'rz']:
                qc.append(gate, qubits)

    qc.measure_all()
    return qc


if __name__ == "__main__":# Example usage
    # Add CX gates
    # Initialize the circuit with sufficient qubits
    qc = QuantumCircuit(20)  # At least as many qubits as the highest index used

    # Now add your gates
    qc.cx(0, 1)
    #  qc.cx(2, 3)
    qc.h(1)
    qc.x(3)
    qc.y(0)
    qc.z(5)
    qc.h(7)
    qc.s(10)
    qc.t(11)

    # Draw the circuit
    print(qc.draw())


    cmap = create_heavy_hex_IBMQ(6, 3)

  #  new_qc = add_identity_to_gates(qc)
  #  print(new_qc)