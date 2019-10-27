import qiskit as qk
import numpy as np
import matplotlib.pylab as plt
qk.IBMQ.load_account()

class QCPairing:
    """
    Class to implement the pairing hamiltonian and estimate its eigenvalues.
    """
    def __init__(self,n_work,n_simulation,delta=1,g=1,dt=0.005,Emax=500):
        """
        Input:
            input:
                n_work (int) - Number of work-qubits
                n_simulation (int) - Number of simulation qubits
                delta (float) - parameter in the one-body hamiltonian
                g (float) - parameter in the two-body hamiltonian
                dt (float) - timestep in the phase estimation algorithm
                Emax (float) - subtracted from the hamiltonian to yield the whole eigenvalue sprectrum.
        """
        self.n_work = n_work
        self.n_simulation = n_simulation
        self.n_qubits = n_work + n_simulation + 1
        self.qb = qk.QuantumRegister(self.n_qubits)
        self.cb = qk.ClassicalRegister(self.n_qubits)
        self.qz = qk.QuantumCircuit(self.qb,self.cb)
        self.delta = delta
        self.g = g
        self.dt = dt
        self.Emax = Emax
    def set_dt(self,dt):
        self.dt = dt
    def set_Emax(self,Emax):
        self.Emax = Emax
    def set_delta(self,delta):
        self.delta = delta
    def set_g(self,g):
        self.g = g

    def H0(self,dt,control_qubit):
        """
        Implements the one-body part of the parining hamiltonian
        Input:
            dt (float) - timestep in the Phase estimation algorithm
            control_qubit (int) - The conditional qubit
        """
        g = self.g
        delta=self.delta
        n_work = self.n_work
        n_simulation=self.n_simulation
        n_qubits=self.n_qubits
        Emax=self.Emax
        qb = self.qb
        cb = self.cb
        qz = self.qz
        s_state = 0
        for q_state in range(0,n_simulation):
            if q_state % 2 == 0:
                s_state += 1
            qz.crz(dt*delta*(s_state - 1),qb[control_qubit],qb[q_state+n_work])

        qz.cu1(-dt*delta*(1/8)*(n_simulation-2)*n_simulation,qb[control_qubit],qb[n_work])
        qz.x(qb[n_work])
        qz.cu1(-dt*delta*(1/8)*(n_simulation-2)*n_simulation,qb[control_qubit],qb[n_work])
        qz.x(qb[n_work])

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
        g = self.g
        delta=self.delta
        n_work = self.n_work
        n_simulation=self.n_simulation
        n_qubits=self.n_qubits
        qb = self.qb
        cb = self.cb
        qz = self.qz
        for p in range(1,n_simulation,2):
            for q in range(p,n_simulation,2):
                if p == q:
                    theta = -2*(1/8)*g*dt
                    qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
                    qz.x(qb[n_work])
                    qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
                    qz.x(qb[n_work])

                    qz.crz(theta,qb[control_qubit],qb[p-1+n_work])
                    qz.crz(theta,qb[control_qubit],qb[p+n_work])

                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                else:
                    theta = -2*(1/16)*g*dt
                    #FIRST TERM:
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    ############
                    #SECOND TERM:
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(np.pi/2,qb[q+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(-np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.rz(-np.pi/2,qb[q+n_work])
                    ###########
                    #THIRD TERM:
                    qz.h(qb[p-1+n_work])
                    qz.rz(np.pi/2,qb[p+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(np.pi/2,qb[q+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(-np.pi/2,qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.rz(-np.pi/2,qb[q+n_work])
                    ###########
                    #FOURTH TERM
                    qz.h(qb[p-1+n_work])
                    qz.rz(np.pi/2,qb[p+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(-np.pi/2,qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(-np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    ###########
                    #FIFTH TERM
                    qz.rz(np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(np.pi/2,qb[q+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.rz(-np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.rz(-np.pi/2,qb[q+n_work])
                    ##########
                    #SIXTH TERM:
                    qz.rz(np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.rz(-np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(-np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    #######################
                    #SEVENTH TERM
                    qz.rz(np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p-1+n_work])
                    qz.rz(np.pi/2,qb[p+n_work])
                    qz.h(qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.rz(-np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(-np.pi/2,qb[p+n_work])
                    qz.h(qb[q-1+n_work])

                    qz.h(qb[q+n_work])
                    ##############
                    #EIGTH TERM:
                    qz.rz(np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p-1+n_work])
                    qz.rz(np.pi/2,qb[p+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(np.pi/2,qb[q+n_work])
                    qz.h(qb[q+n_work])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
                    qz.cx(qb[p-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[p+n_work],qb[n_qubits-1])
                    qz.cx(qb[q-1+n_work],qb[n_qubits-1])
                    qz.cx(qb[q+n_work],qb[n_qubits-1])
                    qz.h(qb[p-1+n_work])
                    qz.rz(-np.pi/2,qb[p-1+n_work])
                    qz.h(qb[p+n_work])
                    qz.rz(-np.pi/2,qb[p+n_work])
                    qz.h(qb[q-1+n_work])
                    qz.rz(-np.pi/2,qb[q-1+n_work])
                    qz.h(qb[q+n_work])
                    qz.rz(-np.pi/2,qb[q+n_work])
        self.qb = qb
        self.cb = cb
        self.qz = qz


    def PhaseEstimation(self,t):
        """
        Phase estimation algorithm
        t (float) - how long to apply the hamiltonian on the simulation qubits (t/dt iterations)
        """
        self.t = t
        dt = self.dt
        g = self.g
        delta=self.delta
        n_work = self.n_work
        n_simulation=self.n_simulation
        n_qubits=self.n_qubits
        qb = self.qb
        cb = self.cb
        qz = self.qz

        #Initialize simulation qubits to superposition of bell states

        for i in range(0,n_simulation,2):
            qz.h(qb[n_work+i])
            qz.cx(qb[n_work+i],qb[n_work+i+1])

        for cq in range(n_work):
            qz.h(qb[cq])

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
        g = self.g
        delta=self.delta
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
        self.qz.measure(self.qb,self.cb)
        job = qk.execute(self.qz, backend = qk.Aer.get_backend('qasm_simulator'), shots=shots)
        result = job.result()
        result = result.get_counts(self.qz)
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