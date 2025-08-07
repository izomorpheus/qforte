from qiskit.qasm3 import load
from qiskit import transpile, QuantumCircuit
from qiskit_aer import AerSimulator, StatevectorSimulator
from qiskit.visualization import plot_histogram, plot_state_city

#qc = QuantumCircuit(2)
#qc.x(0)
#qc.x(1)
#qc.barrier()
#qc.compose(load("sandbox/time_evolution/circuit.qasm3"))
qc = load("sandbox/time_evolution/circuit.qasm3")
qc.measure_all()

simulator = AerSimulator()
statevector_simulator = AerSimulator(method='statevector')
simulator.set_options(shots = 100000)

## statevector simulation
#transpiled_circ_sv = transpile(qc, statevector_simulator)
#result_sv = statevector_simulator.run(transpiled_circ_sv).result()
#statevector = result_sv.get_statevector(transpiled_circ_sv)
#print("STATEVECTOR"  + str(statevector))

# regular simulation
transpiled_circ = transpile(qc, simulator)
result = simulator.run(transpiled_circ).result()
counts = result.get_counts(transpiled_circ)
plot_histogram(counts, title='Bell-State counts')
print(counts)

statevector_simulator = StatevectorSimulator()
qc = load("sandbox/time_evolution/circuit.qasm3")
transpiled_circ_sv = transpile(qc, statevector_simulator)
result_sv = statevector_simulator.run(transpiled_circ_sv).result()
vec = result_sv.get_statevector()
print(vec)