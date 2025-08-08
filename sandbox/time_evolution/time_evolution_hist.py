from qforte.qiskit_api.workflows import readQASM, writeQASM
from qforte.qiskit_api.dispatchers import QpuDispatcher, AerDispatcher, Dispatcher
from qiskit_aer import StatevectorSimulator
from qiskit.visualization import plot_histogram
from qiskit_ibm_runtime.fake_provider import FakeBrisbane
from matplotlib import pyplot as plt
import os
import pickle

# Helper to cache and load results
def get_or_cache(name, dispatcher, circuits, **kwargs):
    cache_file = f"{name}_result.pkl"
    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            return pickle.load(f)
    result = dispatcher.dispatch_sampler(circuits=circuits, **kwargs)
    with open(cache_file, 'wb') as f:
        pickle.dump(result, f)
    return result


def hist(circuit=None, shots=None):

    if not circuit:
    # Load the time evolution circuit from a QASM file
        qc = readQASM("circuit.qasm3")
        if not qc:
            raise ValueError("Failed to read the quantum circuit from the QASM file.")
    else:
        # make a copy so we don’t mutate the original (e.g. add measurements)
        qc = circuit.copy()

    qc.measure_all()

    # Create dispatchers for QPU and Aer
    qpu_dispatcher = QpuDispatcher()
    aer_dispatcher = AerDispatcher()
    aer_statevector_simulator = Dispatcher(StatevectorSimulator())
    fake_dispatcher = Dispatcher(FakeBrisbane())


    # Dispatch the circuits using QPU (cached)
    qpu_result = get_or_cache("qpu", qpu_dispatcher, circuits=[qc], shots=shots)

    # Dispatch the circuits using Aer (cached)
    aer_result = get_or_cache("aer", aer_dispatcher, circuits=[qc], shots=shots)

    aer_statevector_simulator_result = get_or_cache("aer_statevector", aer_statevector_simulator, circuits=[qc], shots=shots)

    fake_dispatcher_result = get_or_cache("fake", fake_dispatcher, circuits=[qc], shots=shots)

    # Extract data from the QPU results
    qpu_counts = qpu_result[0].data.meas.get_counts()

    # Extract data from the Aer results
    aer_counts = aer_result[0].data.meas.get_counts()

    # Extract data from the Aer statevector simulator results
    aer_statevector_counts = aer_statevector_simulator_result[0].data.meas.get_counts()
    fake_sherbrooke_counts = fake_dispatcher_result[0].data.meas.get_counts()

    # plot the qpu results
    plot_histogram(qpu_counts, title="Single Qubit Circuit Counts")

    # Plot the aer results
    plot_histogram(aer_counts, title="Aer Single Qubit Circuit Counts")

    # Plot the aer statevector simulator results
    plot_histogram(aer_statevector_counts, title="Aer Statevector Single Qubit Circuit Counts")

    plot_histogram(fake_sherbrooke_counts, title="Fake Sherbrooke Single Qubit Circuit Counts")

    plt.show()
