import qforte as qf
import numpy as np
from qforte.helper.df_ham_helper import *
from qforte.utils.exponentiate import exponentiate_pauli_string
from qforte.adapters.qiskit_adapters import qforte_to_qiskit_V1, qforte_to_qiskit_V2
from qiskit.visualization import circuit_drawer
from matplotlib import pyplot as plt

##############
# INIT TIMER #
##############

timer = qf.local_timer()
timer.reset()

#################
# INIT MOLECULE #
#################

#molecular geometry
geom = [
    ('H', (0., 0., 1.0)), 
    ('H', (0., 0., 2.0)),
    # ('H', (0., 0., 3.0)), 
    # ('H', (0., 0., 4.0)),
     ]

#molecule object that contains both the fermionic and qubit Hamiltonians.
mol = qf.system_factory(
    build_type='psi4',  #use psi4
    mol_geometry=geom,  #use the geom we specified
    basis='sto-3g',     #use a minimal basis
    build_qb_ham = True,# VERY STRANGE BEHAVIOR!!
    run_fci=1,          #fci energy
    store_mo_ints=1,    #???
    build_df_ham=0,     #double factorized hamiltonian
    df_icut=1.0e-6)     #???

#time the molecule instantiation
timer.record('Run Psi4 and Initialize')

#######################
# INIT TIME EVOLUTION #
#######################

dt = 0.1  #time step
N = 10    #num of steps
r = 1     #trotter number
order = 1 #trotter order

ref = mol.hf_reference   #hartree fock reference state
nel = sum(ref)           #num of electrons
nqubits = len(ref)
norb = int(len(ref) / 2) #num of spatial orbitals
sz = 0                   #z component of spin

gphase = np.exp(-1.0j*dt*mol.nuclear_repulsion_energy) #not sure why we do this?


# print initialization params
print(f"dt:    {dt}")
print(f"r:     {r}") 
print(f"order: {order}")

print(f"nel {nel} norb {norb}")

print(f"\nmol.nuclear_repulsion_energy: {mol.nuclear_repulsion_energy}")
print(f"\nexp(-i*dt*nre): {np.exp(-1.0j*dt*mol.nuclear_repulsion_energy)}")



################################
# TROTTERIZED HAMILTONIAN INIT #
################################

#get the hamiltonian
sqham = mol.sq_hamiltonian 

#extract the hermitian pairs
hermitian_pairs = qf.SQOpPool()
hermitian_pairs.add_hermitian_pairs(1.0, sqham)

# Get the 2nd quantized hamiltonian
sqham = mol.sq_hamiltonian
print("SQ HAM")
print(sqham)
print("\nQB HAM")
print(mol.hamiltonian)
print(mol.nuclear_repulsion_energy)

# Get the hermitian pairs from the hamiltonian
hermitian_pairs = qf.SQOpPool()
hermitian_pairs.add_hermitian_pairs(1.0, sqham)

# Initialize a circuit to hold the pauli strings and a float to hold the global phase
trotter_circ = qf.Circuit()
trotter_phase = 1.0
gphase2 = np.exp((-1.0j*dt)*hermitian_pairs.terms()[0][1].terms()[0][0]*2)
for pair in hermitian_pairs.terms():
    #unpack the pair
    coeff = pair[0]
    sqop = pair[1]

    #transform the hermitian pair to a linear combination of pauli strings
    qbop = sqop.jw_transform()
    #reapply the coefficient (doublecheck you understand this)
    qbop.mult_coeffs(coeff)

    #exponentiate each pauli string and multiply them all together
    for term in qbop.terms():
        #exponentiate the string (as the time evolution operator)
        exp_op = exponentiate_pauli_string(-1j*dt*term[0], term[1])
        #left-multiply the final circuit to create the product of exponentiated pauli strings
        trotter_circ.add_circuit(exp_op[0])
        #extract the phases to a single prefactor
        trotter_phase *= exp_op[1]
        
print(trotter_circ)
    
######################
# FOCK COMPUTER INIT #
######################

#initialize fock computer
c = qf.Computer(nqubits)

#and set it to the hartree fock state
c.hartree_fock(nel) 

#####################
# FCI COMPUTER INIT #
#####################

#initialize fci computers
fc1 = qf.FCIComputer(nel=nel, sz=sz, norb=norb)
fc2 = qf.FCIComputer(nel=nel, sz=sz, norb=norb)

#and set them to the hartree fock state
fc1.hartree_fock()
fc2.hartree_fock()

#######################
# TIME EVOLUTION MAIN #
#######################
print(f"trotter phase: {trotter_phase}, global phase: {gphase}, global phase2: {gphase2}")
#for each time step
for i in range(N):

    #do time evolution on the fock computer
    c.apply_circuit(trotter_circ)
    c.scale(gphase)

    #do time evolution with the taylor series on fci computer
    fc1.evolve_op_taylor(
        sqham,
        dt,
        1.0e-15,
        30,
        False)

    #do time evolution with the trotter approximation on fci computer
    fc2.evolve_pool_trotter(
        hermitian_pairs,
        dt,
        r,
        order,
        antiherm=False,
        adjoint=False)
    
    #not sure what this does ???
    # fc2.scale(gphase)

    #get the energy by taking the expectation value of the hamiltonian

    E3 = np.real(c.direct_op_exp_val(mol.hamiltonian))

    E1 = np.real(fc1.get_exp_val(sqham))
    E2 = np.real(fc2.get_exp_val(sqham))

    #get the difference between the state vectors
    C1 = fc1.get_state_deep()
    dC2 = fc2.get_state_deep()
    dC2.subtract(C1)
    
    #print the energies for the taylor expansion and the trotter approximation, and the norm of the difference in the state vectors
    print(f"t {(i+1)*dt:6.6f} |dC2| {dC2.norm():6.6f} {E1:6.6f} {E2:6.6f} {E3:6.6f}")

print(fc2.str(print_complex=True))
print(c)
print(type(mol.hamiltonian))

#######################
# QISKIT TESTING CODE #
#######################

# Convert the circuit to Qiskit format

try: # Try both versions of the conversion
    qiskit_trotter_circ_v1 = None
    qiskit_trotter_circ_v2 = None
    qiskit_trotter_circ_v1 = qforte_to_qiskit_V1(trotter_circ, nqubits)
    qiskit_trotter_circ_v2 = qforte_to_qiskit_V2(trotter_circ, nqubits)
except Exception as e:
    print(f"Error converting circuit to Qiskit: {e}")

if qiskit_trotter_circ_v1:
    print(qiskit_trotter_circ_v1.draw())
    #draw circuit using matplotlib
    circuit_drawer(qiskit_trotter_circ_v1, output='mpl')
    plt.show()
if qiskit_trotter_circ_v2: print(qiskit_trotter_circ_v2.draw())

# run the resulting circuit

import qiskit.qasm3

qasm3_str = qiskit.qasm3.dumps(qiskit_trotter_circ_v1)

# To save to a file:
with open("circuit.qasm3", "w") as f:
    f.write(qasm3_str)