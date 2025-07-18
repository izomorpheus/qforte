from qforte import Circuit, gate
from qforte.adapters.qiskit_adapters import qforte_to_qiskit_V2

from qiskit.visualization import circuit_drawer
from matplotlib import pyplot as plt

def test_bell_state_conversion():

    #construct the bell state in qforte
    qf_circ = Circuit()
    qf_circ.add_gate(gate("H", 0))
    qf_circ.add_gate(gate("CNOT", 1, 0))

    #convert to qiskit circuit
    qiskit_circ = qforte_to_qiskit_V2(qf_circ, nqubits=2)

    #draw the circuit as text
    print(qiskit_circ.draw())

    #draw circuit using matplotlib
    circuit_drawer(qiskit_circ, output='mpl')
    plt.show()

    import qiskit.qasm3

    qasm3_str = qiskit.qasm3.dumps(qiskit_circ)

    # To save to a file:
    with open("circuit.qasm3", "w") as f:
        f.write(qasm3_str)

if __name__ == "__main__":
    test_bell_state_conversion()