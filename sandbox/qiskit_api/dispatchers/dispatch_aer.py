from qiskit_ibm_runtime import QiskitRuntimeService

class AerDispatcher:
    def __init__(self, circuit, backend_name: str = 'aer_simulator', service=None):
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
        """Instantiate the Sampler primitive using the runtime service on aer_simulator."""
        self.sampler = self.service.sampler(backend=self.backend_name)  # type: ignore

    def dispatch_sampler(self):
        """Run the circuit on the sampler and return the results."""
        if self.sampler is None:
            self.configure_sampler()
        job = self.sampler.run([self.circuit], **self.job_options)  # type: ignore
        return job.result()

    def configure_estimator(self):
        """Instantiate the Estimator primitive using the runtime service on aer_simulator."""
        self.estimator = self.service.estimator(backend=self.backend_name)  # type: ignore

    def dispatch_estimator(self, observables):
        """Run the circuit on the estimator with given observables and return the results."""
        if self.estimator is None:
            self.configure_estimator()
        job = self.estimator.run([self.circuit], [observables], **self.job_options)  # type: ignore
        return job.result()
