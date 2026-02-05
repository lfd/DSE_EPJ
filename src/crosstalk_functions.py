from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from qiskit_aer.noise import NoiseModel, depolarizing_error, phase_damping_error
from topology_functions import *


def find_active_qubits(circuit):
    dag = circuit_to_dag(circuit)
    active_qubits = [qubit.index for qubit in circuit.qubits
                     if qubit not in dag.idle_wires()]
    return active_qubits


def get_active_qubits(circuit):
    active_qubits = set()
    for instr, qargs, _ in circuit.data:
        for qubit in qargs:
            active_qubits.add(qubit._index)  # Access the protected member _index
    return active_qubits

def get_neighbors(q1_index, q2_index, coupling_map):
    # Find neighbors of qubit
    neighbors = set()
    for qubit in [q1_index, q2_index]:

        for source, target in coupling_map:
            if source == qubit:
                neighbors.add(target)
            elif target == qubit:
                neighbors.add(source)

    if q2_index != q1_index:
        neighbors.remove(q1_index)
        neighbors.remove(q2_index)

    return neighbors


def neighboring_noise_model_for_active_qubits_cx(base_error_prob, neighbor_error_prob, coupling_map, active_qubits):
    noise_model = NoiseModel()

    for connection in coupling_map:
        q1, q2 = connection

        # Depolarizing error for qubits directly involved
        direct_error = depolarizing_error(base_error_prob, 2)
        noise_model.add_quantum_error(direct_error, ['cx', 'cz', 'swap'], [q1, q2])

        # Handling neighbors
        for qubit in [q1, q2]:
            neighbors = [target for source, target in coupling_map if source == qubit] + \
                        [source for source, target in coupling_map if target == qubit]
            # Filter to include only active qubits
            neighbors = [neighbor for neighbor in set(neighbors) - {q1, q2} if neighbor in active_qubits]

            # Apply depolarizing error through 'rz' gate to each active neighboring qubit
            for neighbor in neighbors:
                neighbor_error = depolarizing_error(neighbor_error_prob, 1)
                noise_model.add_quantum_error(neighbor_error, ['id', 'x', 'h', 'rz', 'sx'], [neighbor])

    return noise_model, neighbors

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


#Adds noisy identity gates on the neighboring qubits when multiple CNOT (cx) gates are present.
def create_circ_with_ids(transpiled_circ, cmap, n_qubits):
    qubits_used = []

    circ_data = list(transpiled_circ.data)

    qubits_used = []
    for n in range(n_qubits):
        for j in range(len(list(transpiled_circ.data))):
            if list(transpiled_circ.data)[j][1][0]._index == n and list(transpiled_circ.data)[j][
                0].name != 'barrier' and list(transpiled_circ.data)[j][0].name != 'measure':
                qubits_used.append(n)

    active_qubits = list(collections.Counter(qubits_used).keys())
    print(active_qubits)
    # index_max_gates=np.argmax(list(collections.Counter(qubits_used).values()))
    # qubit_max_gates=list(collections.Counter(qubits_used).keys())[index_max_gates]

    circ_dict = {}
    step = dict((el, 0) for el in active_qubits)
    for j in range(len(circ_data)):
        if circ_data[j][0].name != 'measure' and circ_data[j][0].name != 'barrier':
            if j == 0:
                step[circ_data[j][1][0]._index] += 1

                if circ_data[j][0].name == 'rz':
                    angle = float(circ_data[j][0].params[0])
                else:
                    angle = 0

                circ_dict[circ_data[j][1][0]._index, step[circ_data[j][1][0]._index]] = [circ_data[j][0].name,
                                                                                         circ_data[j][1][0]._index,
                                                                                         angle]

            if j > 0:
                if circ_data[j][0].num_qubits == 1:
                    step[circ_data[j][1][0]._index] += 1
                    if circ_data[j][0].name == 'rz':
                        angle = float(circ_data[j][0].params[0])
                    else:
                        angle = 0

                    circ_dict[circ_data[j][1][0]._index, step[circ_data[j][1][0]._index]] = [circ_data[j][0].name,
                                                                                             circ_data[j][1][
                                                                                                 0]._index, angle]

                if circ_data[j][0].num_qubits == 2:
                    id1 = circ_data[j][1][0]._index
                    id2 = circ_data[j][1][1]._index
                    max_step = max([step[id1], step[id2]])
                    step[id1] = max_step + 1
                    step[id2] = max_step + 1

                    if circ_data[j][0].name == 'rz':
                        angle = float(circ_data[j][0].params[0])
                    else:
                        angle = 0

                    circ_dict[id1, step[id1]] = [circ_data[j][0].name, id2, angle]

    all_steps = []
    for key in circ_dict.keys():
        all_steps.append(key[1])

    circ_schedule = dict()
    for count in range(1, max(all_steps) + 1):
        circ_schedule[count] = dict()
        circ_schedule[count]['gates'] = dict()
        circ_schedule[count]['qubits'] = dict()
        circ_schedule[count]['angles'] = dict()

        gates = []
        qubits = []
        angles = []

        for key in circ_dict.keys():
            if key[1] == count:
                gates.append(circ_dict[key][0])
                qubits.append([key[0], circ_dict[key][1]])
                angles.append(circ_dict[key][2])

        circ_schedule[count]['gates'] = gates
        circ_schedule[count]['qubits'] = qubits
        circ_schedule[count]['angles'] = angles

    qc = QuantumCircuit(n_qubits)

    for key in circ_schedule.keys():
        for j in range(len(circ_schedule[key]['qubits'])):
            qubit0 = circ_schedule[key]['qubits'][j][0]

            gate = circ_schedule[key]['gates'][j]
            qubit_list = []
            if gate == 'cx':
                qubit1 = circ_schedule[key]['qubits'][j][1]
                qc.cx(qubit0, qubit1)
                neibors = get_neighbors(qubit0, qubit1, cmap)
                for q in neibors:
                    qc.id(q)

            elif gate == 'sx':
                qc.sx(qubit0)
                neibors = get_neighbors(qubit0, qubit0, cmap)
                for q in neibors:
                    qc.id(q)
            elif gate == 'x':
                qc.x(qubit0)
                neibors = get_neighbors(qubit0, qubit0, cmap)
                for q in neibors:
                    qc.id(q)
            elif gate == 'rz':
                angle = circ_schedule[key]['angles'][j]
                qc.rz(angle, qubit0)
                neibors = get_neighbors(qubit0, qubit0, cmap)
                for q in neibors:
                    qc.id(q)

    qc.measure_all()

    return qc, circ_schedule

#Adds noisy identities when there are multiple CNOT gates in a step.
#Returns the new circuit and the circuit schedule.
def create_circ_with_ids_cxcx(transpiled_circ, cmap, n_qubits):
    qubits_used = []
    circ_data = list(transpiled_circ.data)

    qubits_used = []
    for n in range(n_qubits):
        for j in range(len(list(transpiled_circ.data))):
            if list(transpiled_circ.data)[j][1][0]._index == n and list(transpiled_circ.data)[j][
                0].name != 'barrier' and list(transpiled_circ.data)[j][0].name != 'measure':
                qubits_used.append(n)

    active_qubits = list(collections.Counter(qubits_used).keys())

    circ_dict = {}
    step = dict((el, 0) for el in active_qubits)
    for j in range(len(circ_data)):
        if circ_data[j][0].name != 'measure' and circ_data[j][0].name != 'barrier':
            if j == 0:
                step[circ_data[j][1][0]._index] += 1

                if circ_data[j][0].name == 'rz':
                    angle = float(circ_data[j][0].params[0])
                else:
                    angle = 0

                circ_dict[circ_data[j][1][0]._index, step[circ_data[j][1][0]._index]] = [circ_data[j][0].name,
                                                                                         circ_data[j][1][0]._index,
                                                                                         angle]

            if j > 0:
                if circ_data[j][0].num_qubits == 1:
                    step[circ_data[j][1][0]._index] += 1
                    if circ_data[j][0].name == 'rz':
                        angle = float(circ_data[j][0].params[0])
                    else:
                        angle = 0

                    circ_dict[circ_data[j][1][0]._index, step[circ_data[j][1][0]._index]] = [circ_data[j][0].name,
                                                                                             circ_data[j][1][0]._index,
                                                                                             angle]

                if circ_data[j][0].num_qubits == 2:
                    id1 = circ_data[j][1][0]._index
                    id2 = circ_data[j][1][1]._index
                    max_step = max([step[id1], step[id2]])
                    step[id1] = max_step + 1
                    step[id2] = max_step + 1

                    if circ_data[j][0].name == 'rz':
                        angle = float(circ_data[j][0].params[0])
                    else:
                        angle = 0

                    circ_dict[id1, step[id1]] = [circ_data[j][0].name, id2, angle]

    all_steps = []
    for key in circ_dict.keys():
        all_steps.append(key[1])

    circ_schedule = dict()
    # count=1
    for count in range(1, max(all_steps) + 1):
        # print(count)
        circ_schedule[count] = dict()
        circ_schedule[count]['gates'] = dict()
        circ_schedule[count]['qubits'] = dict()
        circ_schedule[count]['angles'] = dict()
        circ_schedule[count]['cx_count'] = dict()
        circ_schedule[count]['cx_idxs'] = dict()

        gates = []
        qubits = []
        angles = []
        cx_count = 0
        cx_indices = []
        for key in circ_dict.keys():
            circ_dict[key][0]
            if key[1] == count:
                # print(key)
                gates.append(circ_dict[key][0])
                qubits.append([key[0], circ_dict[key][1]])
                angles.append(circ_dict[key][2])
                if circ_dict[key][0] == 'cx':
                    cx_count += 1
                    cx_indices.append([key[0], circ_dict[key][1]])

        circ_schedule[count]['gates'] = gates
        circ_schedule[count]['qubits'] = qubits
        circ_schedule[count]['angles'] = angles
        circ_schedule[count]['cx_count'] = cx_count
        circ_schedule[count]['cx_idxs'] = cx_indices

    # create new circuit with added noisy identities for extra cross-talk noise
    qc = QuantumCircuit(n_qubits)

    for key in circ_schedule.keys():
        for j in range(len(circ_schedule[key]['qubits'])):
            qubit0 = circ_schedule[key]['qubits'][j][0]

            gate = circ_schedule[key]['gates'][j]

            if gate == 'cx':
                qubit1 = circ_schedule[key]['qubits'][j][1]
                qc.cx(qubit0, qubit1)


            elif gate == 'sx':
                qc.sx(qubit0)

            elif gate == 'x':
                qc.x(qubit0)

            elif gate == 'rz':
                angle = circ_schedule[key]['angles'][j]
                qc.rz(angle, qubit0)

        if circ_schedule[key]['cx_count'] >= 2:
            # noisy identity on qubits themselves
            for q in circ_schedule[key]['cx_idxs']:
                qc.id(q[0])
                qc.id(q[1])

                # noisy identity on neighboring qubits
                neibors = get_neighbors(q[0], q[1], cmap)
                for nq in neibors:
                    qc.id(nq)

    qc.measure_all()

    return qc, circ_schedule


def create_circ_with_ids_cxcx_(transpiled_circ, cmap, n_qubits, num_cx):
    circ_data = list(transpiled_circ.data)
    active_qubits = {q.index for instr, qargs, _ in circ_data for q in qargs}

    circ_dict = {}
    step = dict((el, 0) for el in active_qubits)
    for j in range(len(circ_data)):
        if circ_data[j][0].name != 'measure' and circ_data[j][0].name != 'barrier':
            if circ_data[j][0].num_qubits == 1:
                step[circ_data[j][1][0].index] += 1
                circ_dict[circ_data[j][1][0].index, step[circ_data[j][1][0].index]] = circ_data[j]
            elif circ_data[j][0].num_qubits == 2:
                id1 = circ_data[j][1][0].index
                id2 = circ_data[j][1][1].index
                max_step = max([step[id1], step[id2]])
                step[id1] = max_step + 1
                step[id2] = max_step + 1
                circ_dict[id1, step[id1]] = circ_data[j]
                circ_dict[id2, step[id2]] = circ_data[j]

    all_steps = set(step.values())
    circ_schedule = {count: [] for count in all_steps}
    print(circ_schedule)
    print(circ_dict)
    for key in circ_dict.keys():
        print(key)
        print(circ_schedule)
        print(circ_dict[key])
        circ_schedule[key[1]].append(circ_dict[key])

    # create new circuit with added noisy identities for extra cross-talk noise
    qc = QuantumCircuit(n_qubits)
    for count in circ_schedule:
        cx_indices = []
        for instr, qargs, _ in circ_schedule[count]:
            if instr.name == 'cx':
                qc.cx(qargs[0].index, qargs[1].index)
                cx_indices.append((qargs[0].index, qargs[1].index))
            else:
                qc.append(instr, [q.index for q in qargs])

        if len(cx_indices) >= num_cx:
            # noisy identity on qubits themselves
            for q in cx_indices:
                qc.id(q[0])
                qc.id(q[1])

                # noisy identity on neighboring qubits
                neighbors = get_neighbors(q[0], q[1], cmap)

                for nq in neighbors:
                    qc.id(nq)

    qc.measure_all()

    return qc, neighbors


def apply_crosstalk_qc_with_ids(base_error_prob, neighbor_error_prob):
    noise_model = NoiseModel()

    cx_error = depolarizing_error(base_error_prob, 2)
    noise_model.add_all_qubit_quantum_error(cx_error, ['cx', 'cz'])

    id_error = phase_damping_error(neighbor_error_prob)
    noise_model.add_all_qubit_quantum_error(id_error, ['id'])

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


def create_circ_with_ids_cxcx__(transpiled_circ, cmap, n_qubits):
    circ_data = list(transpiled_circ.data)
    active_qubits = {q.index for instr, qargs, _ in circ_data for q in qargs}

    circ_dict = {}
    step = dict((el, 0) for el in active_qubits)

    for j in range(len(circ_data)):
        if circ_data[j][0].name != 'measure' and circ_data[j][0].name != 'barrier':
            if circ_data[j][0].num_qubits == 1:
                step[circ_data[j][1][0].index] += 1
                circ_dict[circ_data[j][1][0].index, step[circ_data[j][1][0].index]] = circ_data[j]
            elif circ_data[j][0].num_qubits == 2:
                id1 = circ_data[j][1][0].index
                id2 = circ_data[j][1][1].index
                max_step = max([step[id1], step[id2]])
                step[id1] = max_step + 1
                step[id2] = max_step + 1
                circ_dict[id1, step[id1]] = circ_data[j]
                circ_dict[id2, step[id2]] = circ_data[j]

    all_steps = set(step.values())
    circ_schedule = {count: [] for count in all_steps}

    print("Active Qubits: ", active_qubits)
    print("Steps: ", step)
    print("Circ Dict: ", circ_dict)
    print("Initialized circ_schedule: ", circ_schedule)

    for key in circ_dict.keys():
        if key[1] not in circ_schedule:
            print(f"Missing key in circ_schedule: {key[1]}")
            circ_schedule[key[1]] = []
        circ_schedule[key[1]].append(circ_dict[key])

    print("Circ Schedule: ", circ_schedule)

    # create new circuit with added noisy identities for extra cross-talk noise
    qc = QuantumCircuit(n_qubits)
    neighbors = set()  # Initialize neighbors to avoid UnboundLocalError
    for count in circ_schedule:
        cx_indices = []
        for instr, qargs, _ in circ_schedule[count]:
            if instr.name == 'cx':
                qc.cx(qargs[0].index, qargs[1].index)
                cx_indices.append((qargs[0].index, qargs[1].index))
            else:
                qc.append(instr, [q.index for q in qargs])

        if len(cx_indices) >= 2:
            # noisy identity on qubits themselves
            for q in cx_indices:
                qc.id(q[0])
                qc.id(q[1])

                # noisy identity on neighboring qubits
                new_neighbors = get_neighbors(q[0], q[1], cmap)
                neighbors.update(new_neighbors)  # Update the set of neighbors

                for nq in new_neighbors:
                    qc.id(nq)

    qc.measure_all()

    return qc, list(neighbors)  # Convert neighbors to list before returning


