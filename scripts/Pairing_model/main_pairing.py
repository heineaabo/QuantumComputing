import matplotlib.pylab as plt
from qcPairingModel import QCPairing
import qiskit as qk
import numpy as np
#Outout from FCI:
#numpy eig:  [-0.61803399  1.61803399] 1pair-4basis eigenvalues with diagonalization
#numpy eig:  [1.] 					   2pair-4basis 	-----||------


for n_work in [4,6,8]:
	n_simulation = 4
	Emax= 2
	dt = 0.005
	g = 1
	delta=1
	t= 100*dt
	Pairing1 = QCPairing(n_work,n_simulation,delta=delta,g=g,dt=dt,Emax=Emax) #Initialization
	measurements = Pairing1.run_simulation(t = t,shots = 1000) #Run simulation
	x,y = Pairing1.sort_measurements(measurements) #Calculates how many times a single eigenvalue is measured
	plt.plot(x,y)
	plt.xlabel('Eigenvalue')
	plt.ylabel('Times measured')
	plt.show()
	











