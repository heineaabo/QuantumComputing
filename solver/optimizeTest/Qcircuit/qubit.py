from gates import *

class Qubit:
    def __init__(self):
        self.circ = []

    def act(self,O):
        appended = False
        if len(self.circ) == 0:
            self.circ.append(O)
            appended = True
            return appended
        elif self.all_Identity():
            self.circ[0] = O
            return appended
        else:
            for i in reversed(range(len(self.circ))):
                if isinstance(self.circ[i],Id):
                    continue
                else:
                    if i == len(self.circ)-1:
                        self.circ.append(O)
                        appended = True
                        break
                    else:
                        self.circ[i+1] = O
                        break
            return appended

    def Identity(self):
        self.circ.append(Id())

    def all_Identity(self):
        b = True
        for i in range(len(self.circ)):
            if not isinstance(self.circ[i],Id):
                b = False
        return b
    
    def single_gate_optimization(self):
        for i in range(len(self.circ)-1):
            if type(self.circ[i]) == type(self.circ[i+1]) and\
               type(self.circ[i]) != type(Id()):
                self.circ[i] = Id()
                self.circ[i+1] = Id()

    def remove(self,operations):
        for o in reversed(operations):
            self.circ.pop(o)

    def __repr__(self):
        return str(self.circ)
