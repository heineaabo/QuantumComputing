import matplotlib.pylab as plt
from Heisenberg import QCheisenberg
import qiskit as qk
import numpy as np
#Outout from FCI:
#numpy eig:  [-0.61803399  1.61803399] 1pair-4basis eigenvalues with diagonalization
#numpy eig:  [1.] 					   2pair-4basis 	-----||------


n_simulation = 3
n_work = 8
Emax= 7
dt = 0.005
h0=1
t= 100*dt
model = QCheisenberg(n_work,n_simulation,h0=h0,dt=dt,Emax=Emax) #Initialization
measurements = model.run_simulation(t = t,shots = 1000) #Run simulation
x,y = model.sort_measurements(measurements) #Calculates how many times a single eigenvalue is measured
plt.plot(x,y)
plt.xlabel('Eigenvalue')
plt.ylabel('Times measured')
plt.show()
    

