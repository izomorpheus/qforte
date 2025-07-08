import qforte as qf
import numpy as np
from qforte.helper.df_ham_helper import *
from qforte.utils.exponentiate import exponentiate_pauli_string

# NOTE(Nick): this sandbox file compares the evolution of the HF wfn
# under trotterized evoltuion of the double factorized hamiltonain  
# to the exact time evolution (via taylor expansion of e^-itH)


def t_diff(Tqf, npt, name, print_both=False):
    print(f"\n  ===> {name} Tensor diff <=== ")
    Tnp = qf.Tensor(shape=np.shape(npt), name='Tnp')
    Tnp.fill_from_nparray(npt.ravel(), np.shape(npt))
    if(print_both):
        print(Tqf)
        print(Tnp)
    Tnp.subtract(Tqf)
    print(f"  ||dT||: {Tnp.norm()}")
    if(Tnp.norm() > 1.0e-12):
        print(Tnp)

geom = [
    ('H', (0., 0., 1.0)), 
    ('H', (0., 0., 2.0)),
    # ('H', (0., 0., 3.0)), 
#     ('H', (0., 0., 4.0)),
#     ('H', (0., 0., 5.0)), 
#     ('H', (0., 0., 6.0)),
#     ('H', (0., 0., 7.0)), 
#     ('H', (0., 0., 8.0)),
#     # ('H', (0., 0., 9.0)), 
#     # ('H', (0., 0.,10.0)),
#     # ('H', (0., 0.,11.0)), 
#     # ('H', (0., 0.,12.0))
     ]

#
# geom = [
#     ('H', (0., 0., 1.0)), 
#     ('Be', (0., 0., 2.0)),
#     ('H', (0., 0., 3.0)), 
#     ]

# geom = [('Li', [0.0, 0.0, 0.0]), ('H', [0.0, 0.0, 1.45])]


timer = qf.local_timer()

timer.reset()
# Get the molecule object that now contains both the fermionic and qubit Hamiltonians.
mol = qf.system_factory(
    build_type='psi4', 
    mol_geometry=geom, 
    basis='sto-3g', 
    build_qb_ham = True, ## VERY STRANGE BEHAVIOR!!
    run_fci=1,
    store_mo_ints=1,
    build_df_ham=0,
    df_icut=1.0e-6)

timer.record('Run Psi4 and Initialize')


## ====> Set up Time Step and number of steps <==== ##
dt = 0.1
N = 10
r = 1 #trotter number
order = 1 #trotter order

## ====> Set up DF and Trotter Stuff <==== ##

sqham = mol.sq_hamiltonian
print()
hermitian_pairs = qf.SQOpPool()
hermitian_pairs.add_hermitian_pairs(1.0, sqham)





##########################################
## BEGIN FOCK TROTTERIZED TIME EVOLUTION #
##########################################

# Initialize a Fock Computer for trotterized time evolution
#trotter_computer = qf.Computer()

# Get the 2nd quantized hamiltonian
sqham = mol.sq_hamiltonian

# Get the hermitian pairs from the hamiltonian
hermitian_pairs = qf.SQOpPool()
hermitian_pairs.add_hermitian_pairs(1.0, sqham)

# Initialize a circuit to hold the pauli strings and a float to hold the global phase
trotter_circ = qf.Circuit()
trotter_phase = 1.0

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
    
#######################################
# END FOCK TROTTERIZED TIME EVOLUTION #
#######################################


































# print(f"sqham: {sqham}")
# print(f"hermitian pairs: {hermitian_pairs}")
#
# print("")
# print(f"len(sqham.terms()):  {len(sqham.terms())} ")
# print(f"len(hermitian_pairs.terms()):  {len(hermitian_pairs.terms())} ")
# print("")


## ====> set up FCIComputers <==== ##
ref = mol.hf_reference

nel = sum(ref)
sz = 0
norb = int(len(ref) / 2)

print(f"nel {nel} norb {norb}")

fc1 = qf.FCIComputer(nel=nel, sz=sz, norb=norb)
fc2 = qf.FCIComputer(nel=nel, sz=sz, norb=norb)

fc1.hartree_fock()
fc2.hartree_fock()





## ====> DF Evolution it a 32-bit integer? it a 32-bit integer? <==== ##
print(f"dt:    {dt}")
print(f"r:     {r}") 
print(f"order: {order}")

print("\n\n")

# print(sqham)

print(f"\nmol.nuclear_repulsion_energy: {mol.nuclear_repulsion_energy}")
print(f"\nexp(-i*dt*nre): {np.exp(-1.0j*dt*mol.nuclear_repulsion_energy)}")
print("")

gphase = np.exp(-1.0j*dt*mol.nuclear_repulsion_energy)

for i in range(N):

    fc1.evolve_op_taylor(
        sqham,
        dt,
        1.0e-15,
        30,
        False)

    fc2.evolve_pool_trotter(
        hermitian_pairs,
        dt,
        r,
        order,
        antiherm=False,
        adjoint=False)
    
    # fc2.scale(gphase)
    

    E1 = np.real(fc1.get_exp_val(sqham))
    E2 = np.real(fc2.get_exp_val(sqham))

    C1 = fc1.get_state_deep()
    dC2 = fc2.get_state_deep()

    dC2.subtract(C1)
    
    print(f"t {(i+1)*dt:6.6f} |dC2| {dC2.norm():6.6f} {E1:6.6f} {E2:6.6f} ")
