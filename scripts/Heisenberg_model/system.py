import numpy as np
import qiskit as qk

class Model:
    def __init__(self,n_work,n_simulation,n_ancilla,Emax,*argv):
        self.w = n_work
        self.s = n_simulation
        self.N = n_work + n_simulation + n_ancilla
        self.qb = qk.QuantumRegister(self.N)
        self.cb = qk.ClassicalRegister(self.N)
        self.qc = qk.QuantumCircuit(self.qb,self.cb)
        self.Emax = Emax
        self.args = argv
        self.plotX = np.zeros(1)
        self.plotY = np.zeros(1)
        self.result = {}
    
    def toPhase(self,bitstr):
        phase = 0
        for i,c in enumerate(bitstr):
            phase += (2**(-i-1))*int(c)
        return phase
        
    def measure(self,t):
        self.qc.measure(self.qb,self.cb)
        job = qk.execute(self.qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=1024)
        #qk.tools.monitor.job_monitor(job)
        result = job.result().get_counts(self.qc)
        x = [] # phase
        y = [] # hits
        for key,val in result.items():
            eigenstate = key[:self.N-self.w]
            phi = key[self.N-self.w:]
            phi = phi[::-1]
            x.append(self.Emax - 2*np.pi*self.toPhase(phi)/t)
            y.append(val)
        x = np.array(x)
        y = np.array(y)
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]
        ## Check if same phase measured for different eigenstates
        x_ = []
        y_ = []
        for i,xi in enumerate(x):
            if i > 0:
                if xi == x_[-1]:
                    y_[-1] += y[i]
                else:
                    x_.append(xi)
                    y_.append(y[i])

            else:
                x_.append(xi)
                y_.append(y[i])
        self.plotX = x_
        self.plotY = y_
        self.result = result
        return self

#####################
### Pairing model ###
#####################
class Pairing(Model):
    def __init__(self,*args):
        super().__init__(*args)
        self.delta = self.args[0]
        self.g = self.args[1]
        self.ansatz()
    
    def __call__(self,control_qubit,dt):
        n_work = self.w
        n_simulation = self.s
        n_qubits = self.N
        g = self.g
        delta = self.delta
        Emax = self.Emax

        qz = self.qc
        qb = self.qb
        
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
        self.qc = qz
        self.qb = qb
        return self
    
    def ansatz(self):
        for i in range(0,self.s,2):
            self.qc.h(self.qb[self.w+i])
            self.qc.cx(self.qb[self.w+i],self.qb[self.w+i+1])
        return None