from qiskit_ibm_runtime import QiskitRuntimeService, Sampler, Estimator

class Dispatcher:
    """
    Base IBM Runtime dispatcher for Sampler and Estimator primitives.
    """
    def __init__(self, circuit, backend_name, service=None):
        self.circuit = circuit
        self.backend_name = backend_name
        self.service = service or QiskitRuntimeService()
        self.job_options = {}
        self.sampler = None
        self.estimator = None

    def configure_job(self, shots: int = 1024, **options):
        """Configure job-level options such as number of shots and runtime options."""
        self.job_options = {'shots': shots, **options}

    def configure_sampler(self):
        """Instantiate the Sampler primitive for the configured backend."""
        # SamplerV2: pass backend mode (name or object)
        self.sampler = Sampler(mode=self.backend_name)

    def dispatch_sampler(self):
        """Run the circuit on the sampler and return the results."""
        if self.sampler is None:
            self.configure_sampler()
        job = self.sampler.run([self.circuit], **self.job_options)
        return job.result()

    def configure_estimator(self):
        """Instantiate the Estimator primitive for the configured backend."""
        # EstimatorV2: pass backend mode (name or object)
        self.estimator = Estimator(mode=self.backend_name)

    def dispatch_estimator(self, observables):
        """Run the circuit on the estimator with given observables and return the results."""
        if self.estimator is None:
            self.configure_estimator()
        job = self.estimator.run([self.circuit], [observables], **self.job_options)
        return job.result()


class SimDispatcher(Dispatcher):
    """
    Simulator dispatcher preconfigured to use the aer_simulator backend.
    """
    def __init__(self, circuit, service=None):
        super().__init__(circuit, backend_name='aer_simulator', service=service)

if __name__ == "__main__":
    from qiskit import QuantumCircuit
    import matplotlib.pyplot as plt  # type: ignore

    # Construct a simple 1-qubit circuit
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.measure_all()

    # Run using the simulator dispatcher
    dispatcher = SimDispatcher(qc)
    dispatcher.configure_job(shots=1024)
    result = dispatcher.dispatch_sampler()

    # Extract probability distribution
    dist = result.quasi_dists[0]  # type: ignore

    # Plot histogram of probabilities
    plt.bar(dist.keys(), dist.values())
    plt.xlabel("Bitstring")
    plt.ylabel("Probability")
    plt.title("SimDispatcher Sampler Results")
    plt.show()
