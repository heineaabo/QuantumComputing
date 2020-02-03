from circuit import *
from qubit import *
from gates import *

I = Id()
x = X()
y = Y()
z = Z()

Q = Circuit(4)
Q.apply(x,1)
Q.apply(y,1)
Q.apply(y,2)
Q.apply(z,2)
Q.apply(z,2)

print(Q)
Q.optimize()
print(Q)
Q.apply(ControlGate(X(),1,3))
print(Q)
