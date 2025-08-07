from qiskit_ibm_runtime import QiskitRuntimeService, Sampler, Estimator
from qiskit import transpile
from qiskit.transpiler import generate_preset_pass_manager

class QPUDispatcher:
    def __init__(self, backend_name=None):
        self.backend_name = backend_name
        self.service = QiskitRuntimeService()
        self.pass_manager = None
        self.sampler = None
        self.estimator = None
        self.backend = self.set_backend_from_name()

    def set_backend_from_name(self):
        try:
            self.backend = self.service.least_busy() if not self.backend_name else self.service.backend(backend_name)
            self.pass_manager = generate_preset_pass_manager(backend=self.backend)
        except Exception as e:
            raise ValueError(f"Could not set backend: {e}")

    def get_backend(self):
        return self.backend

    def set_backend(self, backend_name=None):
        self.backend_name = backend_name
        self.set_backend_from_name()

    def dispatch_sampler(self, circuits, shots=None):
        # Ensure backend is initialized
        if self.backend is None:
            self.set_backend_from_name()
        # Transpile to match hardware connectivity
        transpiled = transpile(circuits, backend=self.backend)
        # Instantiate and run sampler
        self.sampler = Sampler(mode=self.backend)
        job = self.sampler.run(transpiled, shots=shots)
        return job.result()

    def dispatch_estimator(self, circuits, observables, precision=None):
        if self.backend is None:
            self.set_backend_from_name()
        transpiled = transpile(circuits, backend=self.backend)
        pubs = zip(circuits, observables)
        if self.estimator is None:
            self.estimator = Estimator(mode=self.backend)
        job = self.estimator.run(pubs, precision)  # type: ignore
        return job.result()

class AerDispatcher:
    def __init__(self, backend_name=None):
        self.backend_name = backend_name
        self.service = QiskitRuntimeService()
        self.backend = self.service.least_busy() if not backend_name else self.service.backend(backend_name)
        self.sampler = None
        self.estimator = None

    def set_backend_from_name(self):
        try:
            self.backend = self.service.backend(self.backend_name) if self.backend_name else self.service.least_busy()
        except Exception as e:
            raise ValueError(f"Could not set backend: {e}")

    def get_backend(self):
        return self.backend

    def set_backend(self, backend_name=None):
        self.backend_name = backend_name
        self.set_backend_from_name()

    def dispatch_sampler(self, circuits, shots=None):
        # Ensure backend is initialized
        if self.backend is None:
            self.set_backend_from_name()
        # Transpile to match hardware connectivity
        transpiled = transpile(circuits, backend=self.backend)
        # Instantiate and run sampler
        self.sampler = Sampler(mode=self.backend)
        job = self.sampler.run(transpiled, shots=shots)
        return job.result()

    def dispatch_estimator(self, circuits, observables, precision=None):
        if self.backend is None:
            self.set_backend_from_name()
        transpiled = transpile(circuits, backend=self.backend)
        pubs = zip(circuits, observables)
        if self.estimator is None:
            self.estimator = Estimator(mode=self.backend)
        job = self.estimator.run(pubs, precision)  # type: ignore
        return job.result()