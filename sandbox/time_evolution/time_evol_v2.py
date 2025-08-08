import qforte as qf
from qforte.qiskit_api.workflows import SamplerHistFlow
from qforte.qiskit_api.translators import qforte_to_qiskit
from qforte.evolution.trotter import trotter_evolve, hartree_fock

geom = [
    ('H', (0., 0., 1.0)), 
    ('H', (0., 0., 2.0)),
     ]

mol = qf.system_factory(
    build_type='psi4',  #use psi4
    mol_geometry=geom,  #use the geom we specified
    basis='sto-3g',     #use a minimal basis
    df_icut=1.0e-6)     #???

circuit, computer = trotter_evolve(mol, dt=0.1, N=10, r=1, order=1)
qiskit_circuit = hartree_fock(mol, qforte_to_qiskit(circuit, computer.get_nqubit()))
workflow = SamplerHistFlow(computer,
                qiskit_circuit,
                shots=10000)
workflow.plot()