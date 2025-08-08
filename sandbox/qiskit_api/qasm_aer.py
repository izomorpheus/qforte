from qiskit.qasm3 import load
from qiskit import transpile, QuantumCircuit
from qiskit_aer import AerSimulator, StatevectorSimulator
from qiskit.visualization import plot_histogram, plot_state_city

qc = load("circuit.qasm3")
qc.measure_all()

simulator = AerSimulator()
simulator.set_options(shots = 100000000)

# regular simulation
transpiled_circ = transpile(qc, simulator)
result = simulator.run(transpiled_circ).result()
counts = result.get_counts(transpiled_circ)
plot_histogram(counts, title='Bell-State counts')
print(counts)

## statevector simulation
statevector_simulator = StatevectorSimulator()
qc = load("circuit.qasm3")
transpiled_circ_sv = transpile(qc, statevector_simulator)
result_sv = statevector_simulator.run(transpiled_circ_sv).result()
vec = result_sv.get_statevector()
print(vec)
print(type(result_sv))