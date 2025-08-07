import matplotlib.pyplot as plt  # type: ignore
from qiskit.visualization import plot_histogram

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
    plot1 = plt.hist(qpu_counts.values(), bins=len(qpu_counts), alpha=0.5, label='QPU Single Qubit Circuit')
    plot2 = plt.hist(aer_counts.values(), bins=len(aer_counts), alpha=0.5, label='QPU Single Qubit Circuit')
    # plot_histogram(qpu_counts, title="Single Qubit Circuit Counts")
    # plot_histogram(qpu_counts_bell, title="Bell State Circuit Counts")

    # #plot the aer results
    # plot_histogram(aer_counts, title="Aer Single Qubit Circuit Counts")
    # plot_histogram(aer_counts_bell, title="Aer Bell State Circuit Counts")
    plt.show()

if __name__ == "__main__":
    main()
