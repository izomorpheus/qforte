from qiskit import QuantumCircuit
from qforte import Circuit
import qforte as qf
import numpy as np

def qforte_to_qiskit(qforte_circuit: Circuit) -> QuantumCircuit:
    '''Takes a qForte circuit object and returns 
    an equivalent qiskit QuantumCircuit object'''
    qiskit_circuit = QuantumCircuit()
    gates = qforte_circuit.gates()

    Rzy_mat = 1/np.sqrt(2) * np.array([[1j, 1], [1, 1j]])

    gate_map_1qbit = {
        "X": qiskit_circuit.x,
        "Y": qiskit_circuit.y,
        "Z": qiskit_circuit.z,
        "H": qiskit_circuit.h,
        "V": qiskit_circuit.sx,
        "S": qiskit_circuit.s,
        "T": qiskit_circuit.t,
        "I": qiskit_circuit.id,
        "Rzy": qiskit_circuit.unitary
    }

    gate_map_2qbit = {
        "CNOT": qiskit_circuit.cx,
        "cX": qiskit_circuit.cx,
        "aCNOT": qiskit_circuit.cx,
        "acX": qiskit_circuit.cx,
        "cY": qiskit_circuit.cy,
        "cV": qiskit_circuit.csx,
        "SWAP": qiskit_circuit.swap
    }
    
    gate_map_1qbit_with_param = {
        "R": qiskit_circuit.r,
        "Rx": qiskit_circuit.rx,
        "Ry": qiskit_circuit.ry,
        "Rz": qiskit_circuit.rz,
        "rU1": qiskit_circuit.unitary
    }

    gate_map_2qubit_with_param = {
        "A": qiskit_circuit.unitary,
        "cR": qiskit_circuit.unitary,
        "cRz": qiskit_circuit.unitary,
        "rU2": qiskit_circuit.unitary
    }

    for gate in gates:
        if gate.target() == gate.control():
            if gate.gate_id() in gate_map_1qbit:
                if gate.gate_id() == "Rzy": #special case for Rzy
                    qiskit_circuit.unitary(Rzy_mat, [gate.target()], label=gate.gate_id())
                else:
                    gate_map_1qbit[gate.gate_id()](gate.target())
            elif gate.gate_id() in gate_map_1qbit_with_param:
                if gate.gate_id() == "rU1": #special case for rU1
                    rU1_mat = np.array([[]]) #TODO
                    qiskit_circuit.unitary(rU1_mat, [gate.target()], label=gate.gate_id())
                else:
                    gate_map_1qbit_with_param[gate.gate_id()](gate.target(), gate.param())
            else:
                raise ValueError(f"Unsupported one-qubit gate: {gate.gate_id()}")
        else:
            if gate.gate_id() in gate_map_2qbit:
                gate_map_2qbit[gate.gate_id()](gate.control(), gate.target())
            elif gate.gate_id() in gate_map_2qubit_with_param:
                custom_mat = np.array([]) #TODO
                qiskit_circuit.unitary(custom_mat, [gate.control(), gate.target()], label=gate.gate_id())
            else:
                raise ValueError(f"Unsupported two-qubit gate: {gate.gate_id()}")
    
    return qiskit_circuit
