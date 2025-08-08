# from matplotlib import pyplot as plt
# from qiskit.visualization import circuit_drawer

from qforte import Circuit
from qiskit.circuit import QuantumCircuit
from qforte.qiskit_api.translators import qforte_to_qiskit
from qforte.qiskit_api.dispatchers import Dispatcher
from qiskit.visualization import plot_histogram
from qiskit import qasm3

class SamplerHistFlow:
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

if __name__ == "__main__":
    main()
