from gates import *
from qubit import *

class Circuit:

    def __init__(self,n):
        self.mat = [Qubit() for i in range(n)]

    def apply(self,O,i=None):
        if isinstance(O,OneQubitGate):
            if self.mat[i].act(O):
                for j in range(len(self.mat)):
                    if j != i:
                        self.mat[j].Identity()
        if isinstance(O,ControlGate):
            ctrl,targ = O.get_connections()
            for j in range(len(self.mat)):
                self.mat[j].Identity()
            self.mat[ctrl].circ[-1] = CTRL(targ)
            self.mat[targ].circ[-1] = TARG(ctrl)
            #self.control_gate = 

            

    def optimize(self):
        # Remove double action of same gate
        for i in range(len(self.mat)):
            self.mat[i].single_gate_optimization()
        # Remove unessecary Id gates
        id_list = []
        for i in range(len(self.mat[0].circ)):
            if type(self.mat[0].circ[i]) == type(Id()):
                b = True
                for j in range(1,len(self.mat)):
                    if type(self.mat[j].circ[i]) != type(Id()):
                        b = False
                if b:
                    id_list.append(i)
        for i in range(len(self.mat)):
            self.mat[i].remove(id_list)
                

    def __str__(self):
        return str(self.mat)

