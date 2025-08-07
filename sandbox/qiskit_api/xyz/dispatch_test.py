import os
import sys
import matplotlib.pyplot as plt  # type: ignore

# Add the dispatchers directory to the path
script_dir = os.path.dirname(__file__)
api_root = os.path.abspath(os.path.join(script_dir, '..'))
sys.path.insert(0, api_root)

from dispatchers.dispatch_qpu import QPUDispatcher
from qiskit import QuantumCircuit
from qiskit_ibm_runtime import QiskitRuntimeService


def main():


    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.measure_all()

    qc_bell = QuantumCircuit(2, 2)
    qc_bell.h(0)
    qc_bell.cx(0, 1)
    qc_bell.measure_all()

    dispatcher = QPUDispatcher()
    result = dispatcher.dispatch_sampler(circuits=[qc])
    result_bell = dispatcher.dispatch_sampler(circuits=[qc_bell])

    counts = result[0].data.meas.get_counts()
    counts_bell = result_bell[0].data.meas.get_counts()

    plt.hist(counts)
    plt.hist(counts_bell)
    plt.show()

if __name__ == "__main__":
    main()
