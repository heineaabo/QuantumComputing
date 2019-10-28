import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as spst

def basisVisualizer(L,psi):
    '''Given psi=(#)_10, outputs the state in arrows'''
    #ex: |↓|↑↓|↑|↑|
    psi_2 = bin(psi)[2:]
    N  = len(psi_2)
    up = (L-N)*'0'+psi_2
    configStr = "|"
    uparrow   = '\u2191'
    downarrow = '\u2193'
    for i in range(L):
        blank = True
        if up[i] == '1':
            configStr+=uparrow
            blank = False
        if up[i] == '0':
            configStr+=downarrow
            blank = False
        if blank:
            configStr+="_"
        configStr +="|"
    print(configStr)
def countBits(x):
    '''Counts number of 1s in bin(n)'''
    #From Hacker's Delight, p. 66
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    x = x + (x >> 8)
    x = x + (x >> 16)
    return x & 0x0000003F
#helper function to print binary numbers
def binp(num, length=4):
    '''print a binary number without python 0b and appropriate number of zeros'''
    return format(num, '#0{}b'.format(length + 2))[2:]





def makeH(SzList,L,Jxy,Jz,h):
    '''Make a 1D Heisenberg chain of length L with Jxy,Jz and magnetic field h out of an SzList of states'''
    
    basisMap = {}
    stateID  = 0 
    #generate an ordering
    for state in SzList:
        #print(state) #basisVisualizer(L,state)
        basisMap[state] = stateID
        stateID+=1
    nH = stateID
    H = np.zeros([nH,nH])
    #now fill H
    for state in SzList:
        idxA = basisMap[state]
        H[idxA,idxA] += -h*countBits(state) # B field
        for i in range(L):
            j = (i+1)%L #nearest neighbors are hard coded here
            if (state>>i & 1) == (state>>j & 1):#matching bit check
                H[idxA,idxA] += -Jz/4
            else:
                H[idxA,idxA] -= -Jz/4
                mask = 2**(i)+2**j
                stateB= state^mask #this flips the bits at i,j
                idxB = basisMap[stateB]
                H[idxA,idxB]+= -Jxy/2
    #print(np.all(H==H.T)) #check if Hermitian and is coded properly; very slow
    return H



if __name__ == '__main__':
    #syst = makeH(_,2,1,1,1)
    print("The binary number for |up up down> is: ",int('110',2))
    print("In binary the state is:",binp(int('110',2),3))
    print("The state",int('110',2),"can be visualized as:")
    #use helper function to visualize
    basisVisualizer(L=3,psi=int('110',2))
