from qiskit.qasm3 import load
from qiskit import transpile, QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram, plot_state_city

qc = load("circuit.qasm3")
qc.measure_all()

simulator = AerSimulator()
transpiled_circ = transpile(qc, simulator)

result = simulator.run(transpiled_circ).result()
counts = result.get_counts(transpiled_circ)
plot_histogram(counts, title='Bell-State counts')
print(counts)
