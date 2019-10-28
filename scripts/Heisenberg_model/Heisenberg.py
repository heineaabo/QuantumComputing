import qiskit as qk
import numpy as np
import matplotlib.pylab as plt


class QCheisenberg:
	"""
	Class to implement the pairing hamiltonian and estimate its eigenvalues.
	"""
	def __init__(self,n_work,n_simulation,h0=1,dt=0.005,Emax=500):
            """
            Input:
                n_work (int) - Number of work-qubits
                n_simulation (int) - Number of simulation qubits
                h0 (float) - parameter in the one-body hamiltonian
                dt (float) - timestep in the phase estimation algorithm
                Emax (float) - subtracted from the hamiltonian to yield the whole eigenvalue sprectrum.

            """
            self.n_work = n_work
            self.n_simulation = n_simulation
            self.n_qubits = n_work + n_simulation + 1
            self.qb = qk.QuantumRegister(self.n_qubits)
            self.cb = qk.ClassicalRegister(self.n_qubits)
            self.qz = qk.QuantumCircuit(self.qb,self.cb)
            self.h0 = h0
            self.dt = dt
            self.Emax = Emax
	def set_dt(self,dt):
            self.dt = dt
	def set_Emax(self,Emax):
            self.Emax = Emax
	def set_h0(self,h0):
            self.h0 = h0
	
	def H0(self,dt,control_qubit):
            """
            Implements the one-body part of the parining hamiltonian
            Input:
                    dt (float) - timestep in the Phase estimation algorithm
                    control_qubit (int) - The conditional qubit
            """
            h0=self.h0
            n_work = self.n_work
            n_simulation=self.n_simulation
            n_qubits=self.n_qubits
            Emax=self.Emax
            qb = self.qb
            cb = self.cb
            qz = self.qz
            for q_state in range(0,n_simulation):
                qz.crz(dt*h0,qb[control_qubit],qb[q_state+n_work])

            qz.cu1(Emax*dt,qb[control_qubit],qb[n_work])
            qz.x(qb[n_work])
            qz.cu1(Emax*dt,qb[control_qubit],qb[n_work])
            qz.x(qb[n_work])

            self.qb = qb
            self.cb = cb
            self.qz = qz

	def H1(self,dt,control_qubit):
            """
            Implements the two-body part of the pairing hamiltonian
            Input:
                    dt (float) - timestep in the Phase estimation algorithm
                    control_qubit (int) - The conditional qubit
            """
            n_work = self.n_work
            n_simulation = self.n_simulation
            n_qubits = self.n_qubits
            ancilla = n_qubits - 1
            qb = self.qb
            cb = self.cb
            qz = self.qz
            
            for i in range(n_work,ancilla-1):
                # Fourth sum
                qz.cx(qb[i],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.crz(dt,qb[control_qubit],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.cx(qb[i],qb[ancilla])
            #for i in range(n_work,ancilla-1):
                # Third sum
                qz.crz(np.pi/2,qb[i+1],qb[ancilla])
                qz.crz(np.pi/2,qb[i],qb[ancilla]) 
                qz.ch(qb[i],qb[ancilla])
                qz.ch(qb[i+1],qb[ancilla])
                qz.cx(qb[i],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.crz(dt,qb[control_qubit],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.cx(qb[i],qb[ancilla])
                qz.ch(qb[i+1],qb[ancilla])
                qz.ch(qb[i],qb[ancilla])
                qz.crz(-np.pi/2,qb[i],qb[ancilla]) 
                qz.crz(-np.pi/2,qb[i+1],qb[ancilla])
            #for i in range(n_work,ancilla-1):
                # Second sum
                qz.ch(qb[i],qb[ancilla])
                qz.ch(qb[i+1],qb[ancilla])
                qz.cx(qb[i],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.crz(dt,qb[control_qubit],qb[ancilla])
                qz.cx(qb[i+1],qb[ancilla])
                qz.cx(qb[i],qb[ancilla])
                qz.ch(qb[i+1],qb[ancilla])
                qz.ch(qb[i],qb[ancilla])

            self.qb = qb
            self.cb = cb
            self.qz = qz


	def PhaseEstimation(self,t):
            """
            Phase estimation algorithm
            t (float) - how long to apply the hamiltonian on the simulation qubits (t/dt iterations)
            """
            self.t = t
            n_work = self.n_work
            n_simulation = self.n_simulation
            n_qubits = self.n_qubits
            dt = self.dt
            h0 = self.h0
            qb = self.qb
            cb = self.cb
            qz = self.qz

            #Initialize simulation qubits to superposition

            #for i in range(0,n_simulation,2):
                #qz.h(qb[n_work+i])
                #qz.cx(qb[n_work+i],qb[n_work+i+1])
            
            for cq in range(n_work):
                qz.h(qb[cq])
            for cq in range(n_work,n_qubits-1):
                qz.rx(np.pi/2,qb[cq])

            for cq in range(n_work):
                for j in range(int(t/dt)):
                    self.H0((2**cq)*dt,cq)
                    self.H1((2**cq)*dt,cq)

            self.qb = qb
            self.cb = cb
            self.qz = qz

	def inverse_Fourier(self):
            """
            Inverse Quantum Fourier algorithm
            """
            dt = self.dt
            h0=self.h0
            n_work = self.n_work
            n_simulation=self.n_simulation
            n_qubits=self.n_qubits
            qb = self.qb
            cb = self.cb
            qz = self.qz
            for cq in range(int(n_work/2)):
                qz.swap(qb[cq],qb[n_work-cq-1])
            for cq in range(n_work):
                for i in range(cq):
                    qz.cu1(-2*np.pi/(2**(1+cq-i)),qb[i],qb[cq])
                qz.h(qb[cq])
            self.qb = qb
            self.cb = cb
            self.qz = qz

	def solve(self,t=None):
            """
            Call this method to implement and output the qiskit circuit
            input:
                    t (float) - time the hamiltonian acts on the system. (t/dt iterations) t = dt if not set.

            output:
                    qz - qiskit circuit
                    qb - qiskit quantum bits
                    cb - qiskit classical bits
            """
            if t == None:
                t = self.dt
            self.PhaseEstimation(t)
            self.inverse_Fourier()
            return(self.qz,self.qb,self.cb)

	def run_simulation(self,t = None, shots = 1000):
            """
            Runs the simulation shots times and returns
            Input:
                    t (float) - time hamiltonian is applied to the quantum system (t/dt iterations)
                    shots (int) - number of times to run the simulation
            Output:
                    measurements (array) - array containing energy, state and times measured
            """

            qz,qb,cb = self.solve(t)
            print('Implemented')
            print('Measuring...')
            self.qz.measure(self.qb,self.cb)
            job = qk.execute(self.qz, backend = qk.Aer.get_backend('qasm_simulator'), shots=shots)
            result = job.result()
            result = result.get_counts(self.qz)
            print('Measured')
            measurements = []
            for key,value in result.items():
                key_ = key[self.n_simulation+1:]
                eigenstate = key[1:(self.n_simulation+1)]
                eigenstate = eigenstate[::-1]
                decimal = 0
                for i,bit in enumerate(key_):
                    decimal += int(bit)*2**(-i-1)
                if value != 0:
                    measurements.append(np.array([eigenstate, self.Emax-decimal*2*np.pi/self.t, value]))
            
            measurements = np.array(measurements)
            return(measurements)


	def statEig(self,measurements,min_measure=15):
            """
            Finds the estimated eigenvalue and variance by averaging the peaks found with run_simulation
            input:
                    measurements (array) - output from run_simulation
                    min_measure (int) - Minimum measurements of state before it is considered
                    for eigenvalue estimation.
            output:
                    eigenvalues (list) - Estimated eigenvalues
                    varEigs (list) - Estimated variance of eigenvalue approximation
            """

            sumxiyi = 0
            sumyi = 0
            xi_list = []
            eigenvalues = []
            varEigs = []
            minMeasBool = False
            x = measurements[:,1].astype(np.float)
            idx = np.argsort(x)
            y = measurements[:,2].astype(np.int)
            x = x[idx]
            y = y[idx]
            energy_dict = {}

            for xi,yi in zip(x,y):
                energy_dict[xi] = 0
            for xi,yi in zip(x,y):
                energy_dict[xi] += yi

            x = np.array(list(energy_dict.keys()))
            y = np.array(list(energy_dict.values()))
            idx = np.argsort(x)
            x = x[idx]
            y = y[idx]

            for xi, yi in zip(x,y):
                if yi >= min_measure:
                    minMeasBool = True
                    sumxiyi += xi*yi
                    sumyi += yi
                    xi_list.append(xi)
                if minMeasBool and yi < min_measure:
                    minMeasBool = False
                    mu = sumxiyi/sumyi
                    eigenvalues.append(mu)
                    sumxiyi=0
                    sumyi = 0
                    var = 0
                    for val in xi_list:
                        var += (val - mu)**2
                    var/= len(xi_list)
                    varEigs.append(var)
                    xi_list = []
            return(eigenvalues,varEigs)

	def sort_measurements(self,measurements):
            x = measurements[:,1]
            y = measurements[:,2]
            idx = np.argsort(x)
            x = x[idx]
            y = y[idx]
            eigdict = {}
            for xi in x:
                eigdict[xi] = 0
            for xi, yi in zip(x,y):
                eigdict[xi] += int(yi)

            x = np.array(list(eigdict.keys())).astype(np.float)
            idx = np.argsort(x)
            y = np.array(list(eigdict.values())).astype(np.int)
            x = x[idx]
            y = y[idx]
            return(x,y)

