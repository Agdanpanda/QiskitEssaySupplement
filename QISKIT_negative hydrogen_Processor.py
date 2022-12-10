
import numpy as np
# Importing standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, Aer, IBMQ
from qiskit import execute
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from ibm_quantum_widgets import *
from qiskit.exceptions import QiskitError
from qiskit import ClassicalRegister, QuantumRegister
from qiskit.circuit import quantumcircuit
from qiskit.circuit.controlledgate import ControlledGate
from qiskit.providers.ibmq import least_busy
# Loading your IBM Quantum account(s)
#provider = IBMQ.load_account()
from scipy.optimize import minimize
def Z1(theta):
    q = QuantumRegister(2)
    #Create a Classical Register called "c" with 3 bits
    c = ClassicalRegister(2)
    qc = QuantumCircuit(q,c)
    qc.u1(theta[0] , q[0])
    qc.u3(theta[1], -np.pi/2, np.pi/2, q[0])
    qc.u1(theta[2], q[0])
    qc.cx(q[0],q[1])
    qc.u1(theta[3], q[1])
    qc.u3(theta[4], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[5],q[1])
    qc.cx(q[1], q[0])
    qc.u1(theta[6],q[0])
    qc.u1(theta[7],q[1])
    qc.u3(theta[8], -np.pi/2, np.pi/2, q[0])
    qc.u3(theta[9], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[10], q[0])
    qc.u1(theta[11], q[1])
    qc.z(q[0])
    qc.measure(q[0], c[0])
    qc.measure(q[1], c[1])
    shots = T
    max_credits = 3
    job_hpc = execute(qc, backend,  shots = shots, max_credits = max_credits )
    print(job_hpc.job_id())
    try:        
        result_hpc = job_hpc.result()
    except:
        print('Error found: ',job_hpc.error_message() )
    counts12 = result_hpc.get_counts(qc)
    Z = 0 
    if '00' in list(counts12):
        Z = Z + counts12['00']/T
    if '01' in list(counts12):
        Z = Z + counts12['01']/T
    if '10' in list(counts12):
        Z = Z - counts12['10']/T
    if '11' in list(counts12):
        Z = Z - counts12['11']/T
    return Z
def Z2(theta):
    q = QuantumRegister(2)
    c = ClassicalRegister(2)
    qc = QuantumCircuit(q,c)
    qc.u1(theta[0] , q[0])
    qc.u3(theta[1], -np.pi/2, np.pi/2, q[0])
    qc.u1(theta[2], q[0])
    qc.cx(q[0],q[1])
    qc.u1(theta[3], q[1])
    qc.u3(theta[4], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[5],q[1])
    qc.cx(q[1], q[0])
    qc.u1(theta[6],q[0])
    qc.u1(theta[7],q[1])
    qc.u3(theta[8], -np.pi/2, np.pi/2, q[0])
    qc.u3(theta[9], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[10], q[0])
    qc.u1(theta[11], q[1])
    #qc.z(q[0])
    qc.z(q[1])
    qc.measure(q[0], c[0])
    qc.measure(q[1], c[1])
    shots = T
    max_credits = 3
    job_hpc = execute(qc, backend,  shots = shots, max_credits = max_credits )
    print(job_hpc.job_id())
    try:        
        result_hpc = job_hpc.result()
    except:
        print('Error found: ',job_hpc.error_message() )
    counts12 = result_hpc.get_counts(qc)
    Z = 0 
    if '00' in list(counts12):
        Z = Z + counts12['00']/T
    if '01' in list(counts12):
        Z = Z - counts12['01']/T
    if '10' in list(counts12):
        Z = Z + counts12['10']/T
    if '11' in list(counts12):
        Z = Z - counts12['11']/T  
    return Z
def Z3(theta):
    q = QuantumRegister(2)
    c = ClassicalRegister(2)
    qc = QuantumCircuit(q,c)
    qc.u1(theta[0] , q[0])
    qc.u3(theta[1], -np.pi/2, np.pi/2, q[0])
    qc.u1(theta[2], q[0])
    qc.cx(q[0],q[1])
    qc.u1(theta[3], q[1])
    qc.u3(theta[4], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[5],q[1])
    qc.cx(q[1], q[0])
    qc.u1(theta[6],q[0])
    qc.u1(theta[7],q[1])
    qc.u3(theta[8], -np.pi/2, np.pi/2, q[0])
    qc.u3(theta[9], -np.pi/2, np.pi/2, q[1])
    qc.u1(theta[10], q[0])
    qc.u1(theta[11], q[1])
    qc.z(q[0])
    qc.z(q[1])
    qc.measure(q[0], c[0])
    qc.measure(q[1], c[1])
    shots = T
    max_credits = 3
    job_hpc = execute(qc, backend,  shots = shots, max_credits = max_credits )
    print(job_hpc.job_id())
    try:        
        result_hpc = job_hpc.result()
    except:
        print('Error found: ',job_hpc.error_message() )
    counts12 = result_hpc.get_counts(qc)

    Z = 0 
    if '00' in list(counts12):
        Z = counts12['00']/T
    if '01' in list(counts12):
        Z = Z - counts12['01']/T
    if '10' in list(counts12):
        Z = Z - counts12['10']/T
    if '11' in list(counts12):
        Z = Z + counts12['11']/T  
    return Z
def totalEnergy(theta):
    z1Val = Z1(theta)
    z2Val = Z2(theta)
    z3Val = Z3(theta)
    energy = ( -1.0/4.0*( 1.0 - z1Val ) - 1.0/4.0*( 1.0 - z2Val ) + 5.0/32.0*( 1.0 - z1Val - z2Val + z3Val ) )
    print('Current Energy: ' ,energy)
    print('Current Z1 {:.3f}, Z2 {:.3f}, Z3 {:.3f}'.format(z1Val,z2Val,z3Val))
    return energy
apitoken = ''
try :
    IBMQ.enable_account(token=apitoken , hub='')
except:
    pass
provider = IBMQ.get_provider(project='main')
backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits == 5 and not x.configuration().simulator and x.status().operational==True))
print('Backend:',backend)

T = 8192
theta0 = [0, np.pi/2, 0, 0, np.pi/2, 0, 0, 0, np.pi/2, np.pi/2, 0, 0 ]
res = minimize( totalEnergy, theta0, method = 'COBYLA', options = {'xtol':1e-6, 'disp': True})
print('totalEnergy:',res.fun)
print('Theta: ',res.x)


