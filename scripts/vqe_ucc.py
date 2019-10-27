import qiskit as qk
import numpy as np
from scipy.optimize import minimize
from Qoperator import *
qk.IBMQ.load_accounts()



class VQE_UCC:
    def __init__(self,n_qubits,n_fermi,circuit_list,shots=1000):
        self.n_qubits = n_qubits
        self.n_fermi = n_fermi
        self.shots=shots
        self.delta=delta
        self.circuit_list = circuit_list

    def cluster_term(self,i,j,a,b,m,qc,qb,cb):
        n_qubits = self.n_qubits + 1
        if m == 0:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])
            qc.cx(qb[i],qb[n_qubits])
            qc.cx(qb[j],qb[n_qubits])
            qc.cy(qb[a],qb[n_qubits])
            qc.cx(qb[b],qb[n_qubits])
            return(1)
        if m == 1:
            for k in range(i+1,j-1):
            qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cy(qb[i],qb[n_qubits])
            qc.cx(qb[j],qb[n_qubits])
            qc.cy(qb[a],qb[n_qubits])
            qc.cy(qb[b],qb[n_qubits])
            return(1)
        if m == 2:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cx(qb[i],qb[n_qubits])
            qc.cy(qb[j],qb[n_qubits])
            qc.cy(qb[a],qb[n_qubits])
            qc.cy(qb[b],qb[n_qubits])
            return(1)
        if m == 3:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cx(qb[i],qb[n_qubits])
            qc.cx(qb[j],qb[n_qubits])
            qc.cx(qb[a],qb[n_qubits])
            qc.cy(qb[b],qb[n_qubits])
            return(1)
        if m == 4:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cy(qb[i],qb[n_qubits])
            qc.cx(qb[j],qb[n_qubits])
            qc.cx(qb[a],qb[n_qubits])
            qc.cx(qb[b],qb[n_qubits])
            return(-1)
        if m == 5:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cx(qb[i],qb[n_qubits])
            qc.cy(qb[j],qb[n_qubits])
            qc.cx(qb[a],qb[n_qubits])
            qc.cx(qb[b],qb[n_qubits])
            return(-1)
        if m == 6:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cy(qb[i],qb[n_qubits])
            qc.cy(qb[j],qb[n_qubits])
            qc.cy(qb[a],qb[n_qubits])
            qc.cx(qb[b],qb[n_qubits])
            return(-1)
        if m == 7:
            for k in range(i+1,j-1):
                qc.cz(qb[k],qb[n_qubits])

            for l in range(a+1,b-1):
                qc.cz(qb[l],qb[n_qubits])

            qc.cy(qb[i],qb[n_qubits])
            qc.cy(qb[j],qb[n_qubits])
            qc.cx(qb[a],qb[n_qubits])
            qc.cy(qb[b],qb[n_qubits])
            return(-1)


    def ansatz_term(self,i,j,a,b,qc,qb,cb,t):
        n_qubits = self.n_qubits
        #TERM 1
        qc.h(qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.rz(np.pi/2,qb[b])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.h(qb[b])
        qc.rz(-np.pi/2,qb[b])

        #TERM2

        qc.h(qb[i])
        qc.h(qb[j])
        qc.rz(np.pi/2,qb[a])
        qc.h(qb[a])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.rz(-np.pi/2,qb[a])
        qc.h(qb[b])


        #TERM3

        qc.h(qb[i])
        qc.rz(np.pi/2,qb[j])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.h(qb[j])
        qc.rz(-np.pi/2,qb[j])
        qc.h(qb[a])
        qc.h(qb[b])

        #TERM4

        qc.h(qb[i])
        qc.rz(np.pi/2,qb[j])
        qc.h(qb[j])
        qc.rz(np.pi/2,qb[a])
        qc.h(qb[a])
        qc.rz(np.pi/2,qb[b])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])


        qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.h(qb[j])
        qc.rz(-np.pi/2,qb[j])
        qc.h(qb[a])
        qc.rz(-np.pi/2,qb[a])
        qc.h(qb[b])
        qc.rz(-np.pi/2,qb[b])

        #TERM5

        qc.rz(np.pi/2,qb[i])
        qc.h(qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.rz(-np.pi/2,qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.h(qb[b])

        #TERM6

        qc.rz(np.pi/2,qb[i])
        qc.h(qb[i])
        qc.h(qb[j])
        qc.rz(np.pi/2,qb[a])
        qc.h(qb[a])
        qc.rz(np.pi/2,qb[b])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(-t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.h(qb[i])
        qc.rz(-np.pi/2,qb[i])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.rz(-np.pi/2,qb[a])
        qc.h(qb[b])
        qc.rz(-np.pi/2,qb[b])

        #TERM7

        qc.rz(np.pi/2,qb[i])
        qc.h(qb[i])
        qc.rz(np.pi/2,qb[j])
        qc.h(qb[j])
        qc.h(qb[a])
        qc.rz(np.pi/2,qb[b])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.rz(-np.pi/2,qb[i])
        qc.h(qb[j])
        qc.rz(-np.pi/2,qb[j])
        qc.h(qb[a])
        qc.h(qb[b])
        qc.rz(-np.pi/2,qb[b])

        #TERM8

        qc.rz(np.pi/2,qb[i])
        qc.h(qb[i])
        qc.rz(np.pi/2,qb[j])
        qc.h(qb[j])
        qc.rz(np.pi/2,qb[a])
        qc.h(qb[a])
        qc.h(qb[b])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        qc.rz(t[i,j,a-n_fermi,b-n_fermi]/8,qb[n_qubits])

        qc.cx(qb[i],qb[n_qubits])
        qc.cx(qb[j],qb[n_qubits])
        qc.cx(qb[a],qb[n_qubits])
        qc.cx(qb[b],qb[n_qubits])

        for k in range(i+1,j-1):
            qc.cx(qb[k],qb[n_qubits])

        for l in range(a+1,b-1):
            qc.cx(qb[l],qb[n_qubits])


        qc.h(qb[i])
        qc.rz(-np.pi/2,qb[i])
        qc.h(qb[j])
        qc.rz(-np.pi/2,qb[j])
        qc.h(qb[a])
        qc.rz(-np.pi/2,qb[a])
        qc.h(qb[b])

    def wavefunction_ansatz(self,qc,qb,cb,t,measure=False):

        n_qubits = self.n_qubits
        n_fermi = self.n_fermi
        t = t.reshape(self.n_fermi,self.n_fermi,self.n_qubits - self.n_fermi,self.n_qubits - self.n_fermi)

        for i in range(n_fermi,n_qubits):
            qc.x(qb[i])

        for i in range(n_fermi):
            for j in range(i+1,n_fermi):
                for a in range(n_fermi,n_qubits):
                    for b in range(a+1,n_qubits):
                        self.ansatz_term(i,j,a,b,qc,qb,cb,t)

        if measure:
            qc.measure(qb,cb)
            job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
            result = job.result()
            self.result = result.get_counts(qc)



    def energy(self,t):
        n_qubits = self.n_qubits
        circuit_list = self.circuit_list
        E = 0
        for circuit in circuit_list:

            qb = qk.QuantumRegister(n_qubits+1)
            cb = qk.ClassicalRegister(n_qubits+1)
            qc = qk.QuantumCircuit(qb,cb)
            self.wavefunction_ansatz(qc,qb,cb,t)
            indices = []
            I_count = 0
            for i in range(n_qubits):
                operation = circuit.get(i).op
                if operation == '':
                    I_count += 1
                    continue
                eval('qc.{}(qb[{}])'.format(operation,i))

                if operation == 'x':
                    qc.h(qb[i])
                if operation == 'y':
                    qc.u1(-np.pi/2,qb[i])
                    qc.h(qb[i])
                indices.append(i)

            if I_count < n_qubits:
                qc.measure(qb,cb)
                job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
                result = job.result()
                self.result = result.get_counts(qc)
                H0 = 0
                for key, value in self.result.items():
                    key1 = key[::-1]
                    eigenval = 1
                    for ind in indices:
                        e =  1 if key1[ind] == '0' else -1
                        eigenval *= e
                    H0 += eigenval*value
                H0 /= self.shots

                E += H0*circuit.factor.real
            else:
                E += circuit.factor.real

        qb = qk.QuantumRegister(n_qubits+1)
        cb = qk.ClassicalRegister(n_qubits+1)
        qc = qk.QuantumCircuit(qb,cb)
        self.wavefunction_ansatz(qc,qb,cb,t,measure=True)
        print('----------')
        print(self.result)
        print(E)

        return(E)


    def non_gradient_optimization(self):
        t = np.random.rand(self.n_fermi*self.n_fermi*(self.n_qubits - self.n_fermi)*(self.n_qubits - self.n_fermi))

        res = minimize(self.energy, t, method='COBYLA', options={'disp': True},tol=1e-12)
        print(res.x)

        print(self.result)

    def energy_gradient(self,t):
        n_qubits = self.n_qubits
        n_fermi = self.n_fermi
        circuit_list = self.circuit_list
        E = 0

        t_grad = np.zeros(t.shape)

        for k in range(n_fermi):
            for l in range(k+1,n_fermi):
                for m in range(n_fermi,n_qubits):
                    for n in range(m+1,n_qubits):

                        for circuit in circuit_list:
                            t_grad_temp = 0
                            for term in range(8):
                                qb = qk.QuantumRegister(n_qubits+2)
                                cb = qk.ClassicalRegister(n_qubits+2)
                                qc = qk.QuantumCircuit(qb,cb)

                                qc.h(qb[n_qubits+1])
                                for i in range(n_fermi,n_qubits):
                                    qc.x(qb[i])

                                for i in range(n_fermi):
                                    for j in range(i+1,n_fermi):
                                        for a in range(n_fermi,n_qubits):
                                            for b in range(a+1,n_qubits):
                                                self.ansatz_term(i,j,a,b,qc,qb,cb,t)
                                                if i == k and l == j and m == a and n == b:
                                                    fact = self.cluster_term(i,j,a,b,term,qc,qb,cb)

                                qc.h(qb[n_qubits+1])
                                qc.rx(np.pi/2,qb[n_qubits+1])

                                for q_bit in range(n_qubits):
                                    operation = circuit.get(q_bit).op
                                    if operation == '':
                                        continue
                                    eval('qc.c{}(qb[{}],qb[n_qubits+1])'.format(operation,q_bit))

                                qc.measure(qb,cb)
                                job = qk.execute(qc, backend = qk.Aer.get_backend('qasm_simulator'), shots=self.shots)
                                result = job.result()
                                self.result = result.get_counts(qc)
                                shots = 0
                                for key, value in self.result.items():
                                    key1 = key[::-1]
                                    if key1[-1] == '0':
                                        shots += value

                                prob = shots/self.shots

                                prob = -2*prob + 1

                                t_grad_temp += prob*fact/8
                            t_grad[k,l,m-n_fermi,n-n_fermi] += circuit.factor.real*t_grad_temp
        t_grad *= 2
        return(t_grad)

    def gradient_descent(self,iters=1000,learning_rate=1e-2):
        t = np.random.rand(self.n_fermi,self.n_fermi,(self.n_qubits - self.n_fermi),(self.n_qubits - self.n_fermi))	
        self.energy(t)
        print('-----')	
        for i in range(iters):
            grads = self.energy_gradient(t)
            t = t - learning_rate*grads
            self.energy(t)
            print('--------')


if __name__ == '__main__':
    n_qubits = 4
    n_fermi = 2
    delta = 1
    g = 1
    g /= 4

    h_pq = np.identity(n_qubits)
    for p in range(n_qubits):
        h_pq[p,p] *= delta*(p - (p%2))/2

    h_pqrs = np.zeros((n_qubits,n_qubits,n_qubits,n_qubits))
    for p in range(0,n_qubits-1,2):
        for r in range(0,n_qubits-1,2):
            h_pqrs[p,p+1,r,r+1] = -0.5*g

    Pairing = Hamiltonian(n_qubits)

    circuit_list = Pairing.get_circuits(h_pq,h_pqrs)
    for oplist in circuit_list:
        print(oplist.factor)
        for i in range(n_qubits):
            print('qb[{}]'.format(i),oplist.get(i).op)

    test = VQE_UCC(n_qubits=n_qubits,n_fermi=n_fermi,circuit_list=circuit_list)
    test.gradient_descent(learning_rate=1e-1)