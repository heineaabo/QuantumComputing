import numpy as np
import qiskit as qk
#from systems import Model

def InvFourierTransform(model):
    qc = model.qc
    qb = model.qb
    w = model.w
    
    for i in range(int(w/2)):
        qc.swap(qb[i],qb[w-i-1])   
    for i in range(w):
        for j in range(i):
            qc.cu1(-(2*np.pi)/(2**(i+1-j)),qb[j],qb[i])
        qc.h(qb[i])

    model.qc = qc
    model.qb = qb
    return model

def PhaseEstimation(model,t=0.5,dt=0.005):
    s = model.s
    w = model.w
    
    # Initialize / Create superposition
    for i in range(w):
        model.qc.h(model.qb[i])
    
    # model ansatz in model.__init__
        
    # Apply controlled-U operations
    for i in range(model.w):
        for n in range(int(t/dt)):
            model = model(i,(2**i)*dt)
    
    # Inverse Quantum Fourier Transform
    model = InvFourierTransform(model)
    
    return model