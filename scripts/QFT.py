import qiskit as qk
import numpy as np
qk.IBMQ.load_account()

def QFT(Qcircuit, inverse=False):
    """ _________________________
        
        Quantum Fourier Transform
        _________________________
        
        Input: 
        
            Qcircuit = [qc,qr,cr,n]
                - qc -> Quantum circuit object
                - qr -> Quantum register object
                - cr -> Classical register object
                - n  -> Number of qubits
                
            inverse:
                True,False
   
        Output:
        
            Qcircuit
    """
    
    qc       =  Qcircuit[0]
    qr       =  Qcircuit[1]
    cr       =  Qcircuit[2]
    n_qubits =  Qcircuit[3]
    if not inverse:
        for i in range(n_qubits):
            qc.h(qr[i])
            for j in range(i+1,n_qubits):
                qc.cu1(np.pi/2**(j-i),qr[j],qr[i])

        for i in range(int(n_qubits/2)):
            qc.swap(qr[i],qr[-(i+1)])
    else:
        for i in range(int(n_qubits/2)):
            qc.swap(qr[i],qr[-(i+1)])
            
        for i in range(n_qubits):
            for j in range(i):
                qc.cu1(-np.pi/2**(i-j),qr[j],qr[i])
            qc.h(qr[i])    
    
    return [qc,qc,cr,n_qubits] 