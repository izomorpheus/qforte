# from matplotlib import pyplot as plt
from qforte.qiskit_api.dispatchers import QpuDispatcher, AerDispatcher, Dispatcher
from qiskit.visualization import plot_histogram
from qiskit_ibm_runtime.fake_provider import FakeTorino
from matplotlib import pyplot as plt
import numpy as np
import os
import pickle

from qforte import Circuit
from qiskit.circuit import QuantumCircuit
from qforte.qiskit_api.translators import qforte_to_qiskit
from qforte.qiskit_api.dispatchers import Dispatcher
from qiskit.visualization import plot_histogram
from qiskit import qasm3

class SamplerHistFlowV0:
    def __init__(self, circuits, dispatchers):
        self.qforte_circuits = [circuit for circuit in circuits if isinstance(circuit, Circuit)]
        self.qiskit_circuits = [circuit for circuit in circuits if isinstance(circuit, QuantumCircuit)]
        self.dispatchers = [dispatcher for dispatcher in dispatchers if isinstance(dispatcher, Dispatcher)]
        self.results = []

        for dispatcher in self.dispatchers:
            dispatcher.dispatch_sampler(self.qiskit_circuits)


    def run(self):
        pass

    def get_data(self):
        pass
    def histogram(self, circuit, shots=1024):

        result = self.dispatcher.dispatch_sampler(circuits=[qc], shots=1024)
        return plt


def writeQASM(circuit, f=None):

    qasm_str = qasm3.dumps(circuit)
    if f:
        try:
            fd = open(f, "w")
            fd.write(qasm_str)
        except Exception as e:
            raise RuntimeError(f"Failed to write QASM to file: {e}")

    return qasm_str

def readQASM(f):
    if f:
        try:
            circuit = qasm3.load(f)
            return circuit
        except Exception as e:
            raise RuntimeError(f"Failed to read QASM from file: {e}")


# Add the dispatchers directory to the path
# script_dir = os.path.dirname(__file__)
# api_root = os.path.abspath(os.path.join(script_dir, '..'))
# sys.path.insert(0, api_root)

from qforte.qiskit_api.dispatchers import QpuDispatcher
from qforte.qiskit_api.dispatchers import AerDispatcher
from qiskit import QuantumCircuit

def main():

    # Create a trivial quantum circuit
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.measure_all()

    # Create a Bell state circuit
    qc_bell = QuantumCircuit(2, 2)
    qc_bell.h(0)
    qc_bell.cx(0, 1)
    qc_bell.measure_all()

    # Create dispatchers for QPU and Aer
    qpu_dispatcher = QpuDispatcher()
    aer_dispatcher = AerDispatcher()

    #dispatch the circuits using QPU
    qpu_result = qpu_dispatcher.dispatch_sampler(circuits=[qc], shots=1024)
    qpu_result_bell = qpu_dispatcher.dispatch_sampler(circuits=[qc_bell], shots=1024)

    #dispatch the circuits using Aer
    aer_result = aer_dispatcher.dispatch_sampler(circuits=[qc], shots=1024)
    aer_result_bell = aer_dispatcher.dispatch_sampler(circuits=[qc_bell], shots=1024)

    #extract data from the QPU results
    qpu_counts = qpu_result[0].data.meas.get_counts()
    qpu_counts_bell = qpu_result_bell[0].data.meas.get_counts()

    #extract data from the Aer results
    aer_counts = aer_result[0].data.meas.get_counts()
    aer_counts_bell = aer_result_bell[0].data.meas.get_counts()
    for count in qpu_counts:
        print(count)

    #plot the qpu results
    plot_histogram(qpu_counts, title="Single Qubit Circuit Counts")
    plot_histogram(qpu_counts_bell, title="Bell State Circuit Counts")

     #plot the aer results
    plot_histogram(aer_counts, title="Aer Single Qubit Circuit Counts")
    plot_histogram(aer_counts_bell, title="Aer Bell State Circuit Counts")
    plt.show()


class SamplerHistFlow:

    def __init__(self, computer, circuit, shots=1024):
        self.circuit = circuit
        self.computer = computer
        self.shots = shots
        self.results = []

    def plot(self, legend=True):
        return self.hist()

    def get_or_cache(self, name, dispatcher, circuits, **kwargs):
        cache_file = f"{name}_result.pkl"
        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
        result = dispatcher.dispatch_sampler(circuits=circuits, **kwargs)
        with open(cache_file, 'wb') as f:
            pickle.dump(result, f)
        return result


    def hist(self):

        if not self.circuit:
        # Load the time evolution circuit from a QASM file
            qc = readQASM("circuit.qasm3")
            if not qc:
                raise ValueError("Failed to read the quantum circuit from the QASM file.")
        else:
            # make a copy so we don’t mutate the original (e.g. add measurements)
            qc = self.circuit.copy()

        qc.measure_all()

        # Create dispatchers for QPU and Aer
        qpu_dispatcher = QpuDispatcher("ibm_torino")
        aer_dispatcher = AerDispatcher()
        # aer_statevector_simulator = Dispatcher(StatevectorSimulator())
        fake_dispatcher = Dispatcher(FakeTorino())


        # Dispatch the circuits using QPU (cached)
        qpu_result = self.get_or_cache("qpu", qpu_dispatcher, circuits=[qc], shots=self.shots)

        # Dispatch the circuits using Aer (cached)
        aer_result = self.get_or_cache("aer", aer_dispatcher, circuits=[qc], shots=self.shots)

        # aer_statevector_simulator_result = get_or_cache("aer_statevector", aer_statevector_simulator, circuits=[qc], shots=shots)

        fake_dispatcher_result = self.get_or_cache("fake", fake_dispatcher, circuits=[qc], shots=self.shots)

        # Extract data from the QPU results
        qpu_counts = qpu_result[0].data.meas.get_counts()

        # Extract data from the Aer results
        aer_counts = aer_result[0].data.meas.get_counts()

        # Extract data from the Aer statevector simulator results
        # aer_statevector_counts = aer_statevector_simulator_result[0].data.meas.get_counts()
        fake_brisbane_counts = fake_dispatcher_result[0].data.meas.get_counts()

        # plot the qpu results
        # plot_histogram(qpu_counts, title="Single Qubit Circuit Counts")

        # # Plot the aer results
        # plot_histogram(aer_counts, title="Aer Single Qubit Circuit Counts")

        # # Plot the aer statevector simulator results
        # plot_histogram(aer_statevector_counts, title="Aer Statevector Single Qubit Circuit Counts")

        # plot_histogram(fake_sherbrooke_counts, title="Fake Sherbrooke Single Qubit Circuit Counts")
        q_legend = ['IBM Heron r1', 'Simulated IBM Heron r1', 'Aer Simulator']
        q_dists = [qpu_counts, fake_brisbane_counts, aer_counts]

        if self.computer:
            coeffs = self.computer.get_coeff_vec()
            probs = np.abs(coeffs) ** 2
            probs = [int(np.round(p) * self.shots) for p in probs]
            n = int(np.log2(len(probs)))
            bitstrings = [format(i, f'0{n}b') for i in range(len(probs))]
            amp_probs = dict(zip(bitstrings, probs))
            c_legend = ['QForte Simulator']
            c_dist = [amp_probs]
            legend = q_legend + c_legend
            dists = q_dists + c_dist
            print(dists)
        else:
            legend = q_legend
            dists = q_dists


        ax = plot_histogram(dists,
                        title=r'Trotterized Time Evolution: $H_2 \to \hat{H} \to \prod_j\;  e^{-i H_j t}$',
                        legend=legend,
                        figsize=(10, 6),
                        bar_labels=False)
        plt.savefig("trotter_hist.png", dpi=300, bbox_inches='tight')
        plt.show()
        return dists