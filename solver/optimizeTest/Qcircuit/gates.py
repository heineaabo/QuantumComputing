

class Gate:
    def __init__(self):
        self.char = ''

    def __repr__(self):
        return self.char

class OneQubitGate(Gate):
    def __init__(self):
        super().__init__()

class TwoQubitGate(Gate):
    def __init__(self):
        super().__init__()

class Id(OneQubitGate):
    def __init__(self):
        super().__init__()
        self.char = 'I'


class X(OneQubitGate):
    def __init__(self):
        super().__init__()
        self.char = 'X'

class Y(OneQubitGate):
    def __init__(self):
        super().__init__()
        self.char = 'Y'

class Z(OneQubitGate):
    def __init__(self):
        super().__init__()
        self.char = 'Z'

class H(OneQubitGate):
    def __init__(self):
        super().__init__()
        self.char = 'H'


class ControlGate(Gate):
    def __init__(self,gate,ctrl,targ):
        self.connection = None
        self.gate = gate
        self.char = 'C'+str(gate)
        self.targ = targ
        self.ctrl = ctrl

    def get_connections(self):
        return self.ctrl, self.targ
        

#class CNOT(TwoQubitGate):
#    def __init__(self,ctrl,targ):
#        super().__init__()
#        self.char = 'CX'
#        self.ctrl = CTRL(targ)
#        self.targ = TARG(ctrl)

class CTRL(Gate):
    def __init__(self,targ):
        self.connection = targ
        self.char = '\u2A00 '

class TARG(Gate):
    def __init__(self,ctrl):
        self.connection = ctrl
        self.char = '\u2A01 '


