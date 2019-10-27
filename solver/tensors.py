import numpy as np

def tensor_prod(A,B):
    assert A.shape == B.shape
    n,m = A.shape
    T = np.zeros((n*n,m*m),dtype=complex)
    for i in range(n*n):
        for j in range(m*m):
            T[i,j] = B[i%2,j%2]

    for i in range(n*n):
        for j in range(m*m):
            if i < n and j < m:
                a = A[0,0]
            if i < n and j >= m:
                a = A[0,1]
            if i >= n and j < m:
                a = A[1,0]
            if i >= n and j >= m:
                a = A[1,1]
            T[i,j] = a*T[i,j]
    return(T)

def tprint(tensor):
    n,m = tensor.shape
    for i in range(n):
        elems = []
        for j in range(m):
            a = tensor[i,j]
            if a.real == 0 and a.imag == 0:
                elems.append(' 0.0 ')
            elif a.real != 0 and a.imag == 0:
                if a.real > 0:
                    b = ' ' + str(a.real)
                else:
                    b = str(a.real)
                b += ' '
                elems.append(b)
            elif a.real == 0 and a.imag != 0:
                if a.imag > 0 and a.imag == 1.0:
                    b = '  i  '
                if a.imag > 0 and a.imag != 1.0:
                    b = ' ' + str(a.imag)+'i'
                if a.imag < 0 and a.imag == -1.0:
                    b = ' -i  '
                if a.imag < 0 and a.imag != -1.0:
                    b = str(a.imag)+'i'
                elems.append(b)
            else:
                b = str(a)
                elems.append(b)
        if n == 4:
            print('| {:<4} {:<4} {:<4} {:<4} |'.format(*elems))
        if n == 2:
            print('| {:<4} {:<4} |'.format(*elems))

