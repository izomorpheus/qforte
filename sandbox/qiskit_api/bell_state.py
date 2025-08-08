from qforte import Circuit, gate
from qforte.qiskit_api.translators import qforte_to_qiskit
from qiskit.qasm3 import dumps

def test_bell_state_conversion():

    #construct the bell state in qforte
    qf_circ = Circuit()
    qf_circ.add_gate(gate("H", 0))
    qf_circ.add_gate(gate("CNOT", 1, 0))

    #convert to qiskit circuit
    qiskit_circ = qforte_to_qiskit(qf_circ, 2)

    #draw the circuit in the terminal
    print(qiskit_circ.draw())

    # Convert circuit to QASM string
    qasm3_str = dumps(qiskit_circ)

    # Save the QASM string to a file
    with open("circuit.qasm3", "w") as f:
        f.write(qasm3_str)

if __name__ == "__main__":
    test_bell_state_conversion()