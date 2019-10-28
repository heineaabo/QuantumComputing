import numpy as np

def not_eye(n):
    M = np.zeros((n,n),dtype=complex)
    for i in range(n):
        M[i, n-i-1] = 1.0
    return M

def pauli(n,i = np.complex(0,1)):
    if n == 'x':
        r = np.array([[0,1],[1,0]],dtype=complex)
    if n == 'y':
        r = np.array([[0,-i],[i,0]],dtype=complex)
    if n == 'z':
        r = np.array([[1,0],[0,-1]],dtype=complex)
    return r

class MatrixOperator:
    def __init__(self,matrix=np.zeros((1,0))):
        self.shape = [2,2]
        self.matrix = np.zeros(self.shape,dtype=complex)
        self.i = np.complex(0,1)
        if matrix.shape != (1,0):
            self.set_matrix(matrix)

    def set_matrix(self,A):
        n,m = A.shape
        self.shape = [n,m]
        self.matrix = A

    def __getitem__(self,i):
        return self.matrix[i[0],i[1]]

    def __neg__(self):
        new = MatrixOperator(matrix=-self.matrix)
        return new

    def __pos__(self):
        return self

    def __eq__(self,other):
        assert isinstance(other,MatrixOperator)
        return np.array_equal(self.matrix,other.matrix)

    def __mul__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            new = MatrixOperator(matrix=np.matmul(self.matrix,other.matrix))
        elif type(other) == type(self.matrix):
            new = MatrixOperator(matrix=np.matmul(self.matrix,other))
        elif type(other) in [int,float]:
            new = MatrixOperator(matrix=other*self.matrix)
        return new 

    def __rmul__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            new = MatrixOperator(matrix=np.matmul(other.matrix,self.matrix))
        elif type(other) == type(self.matrix):
            new = MatrixOperator(matrix=np.matmul(other,self.matrix))
        elif type(other) in [int,float]:
            new = MatrixOperator(matrix=other*self.matrix)
        return new 

    def __add__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            new = MatrixOperator(matrix=self.matrix + other.matrix)
        elif type(other) == type(self.matrix):
            new = MatrixOperator(matrix=self.matrix + other)
        return new

    def __radd__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            new = MatrixOperator(matrix=other.matrix + self.matrix)
        elif type(other) == type(self.matrix):
            new = MatrixOperator(matrix=other + self.matrix)
        return new

    def __sub__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            new = MatrixOperator(matrix=self.matrix - other.matrix)
        elif type(other) == type(self.matrix):
            new = MatrixOperator(matrix=self.matrix - other)
        return new

    def dagger(self):
        return MatrixOperator(matrix = self.matrix.conj().T)

    def __matmul__(self,other):
        new = None
        if  isinstance(other,MatrixOperator):
            assert self.shape == other.shape
            n,m = self.shape
            T = np.zeros((n*n,m*m),dtype=complex)
            for i in range(n*n):
                for j in range(m*m):
                    T[i,j] = other.matrix[i%n,j%m]

            for i in range(n*n):
                for j in range(m*m):
                    if i < n and j < m:
                        a = self.matrix[0,0]
                    if i < n and j >= m:
                        a = self.matrix[0,1]
                    if i >= n and j < m:
                        a = self.matrix[1,0]
                    if i >= n and j >= m:
                        a = self.matrix[1,1]
                    T[i,j] = a*T[i,j]
            #self.shape = [n*n,m*m]
            #self.matrpix = T
            new = MatrixOperator(matrix=T)
        elif type(other) == type(self.matrix):
            assert self.matrix.shape == other.shape
            n,m = self.shape
            T = np.zeros((n*n,m*m),dtype=complex)
            for i in range(n*n):
                for j in range(m*m):
                    T[i,j] = other[i%n,j%m]

            for i in range(n*n):
                for j in range(m*m):
                    if i < n and j < m:
                        a = self.matrix[0,0]
                    if i < n and j >= m:
                        a = self.matrix[0,1]
                    if i >= n and j < m:
                        a = self.matrix[1,0]
                    if i >= n and j >= m:
                        a = self.matrix[1,1]
                    T[i,j] = a*T[i,j]
            #self.shape = [n*n,m*m]
            #self.matrpix = T
            new = MatrixOperator(matrix=T)
        return new

    def __str__(self):
        strings, longest = self._convert_matrix()
        n,m = self.shape
        full = '\n'
        upline = chr(9122)
        for i in range(n):
            if i == 0:
                full += '{} '.format(chr(9121))
            elif i == n-1:
                full += '{} '.format(chr(9123))
            else:
                full += '{} '.format(chr(9122))
            for j in range(m):
                if longest > 5:
                    full += '{:^7}'.format(strings[i][j])
                elif longest > 3:
                    full += '{:^5}'.format(strings[i][j])
                else:
                    full += '{:^3}'.format(strings[i][j])
            if i == 0:
                full += ' {}\n'.format(chr(9124))
            elif i == n-1:
                full += ' {}\n'.format(chr(9126))
            else:
                full += ' {}\n'.format(chr(9125))
        return full

    def _convert(self,number):
        real = number.real
        im = number.imag
        string = ''
        if np.isclose(real,0) != True and np.isclose(im,0) != True:
            return str(number)
        elif np.isclose(real,0) == False and np.isclose(im,0):
            if real < 0:
                string += '-'
            string += str(abs(np.around(real,4)))
            return string
        elif np.isclose(real,0) and np.isclose(im,0) == False:
            if im < 0:
                string += '-'
            if np.isclose(abs(im),1) == False and np.isclose(abs(im),0) == False:
                string += str(abs(im))
            string += 'i'
            return string
        else:
            return str(0)

    def _convert_matrix(self):
        strings = []
        longest = 0
        n,m = self.shape
        fact,new = self._check_fact()
        if fact == 0:
            for i in range(n):
                strings.append([])
                for j in range(m):
                    elem = self._convert(self.matrix[i,j])
                    strings[i].append(elem)
                    if len(elem) > longest:
                        longest = len(elem)
        else:
            print('{} x'.format(fact))
            for i in range(n):
                strings.append([])
                for j in range(m):
                    elem = self._convert(new)
                    strings[i].append(elem)
                    if len(elem) > longest:
                        longest = len(elem)
        return strings, longest

    def _check_fact(self):
        n,m = self.shape
        x = self._get_fact(self.matrix[0,0])
        new = np.zeros_like(self.matrix)
        for i in range(n):
            for j in range(m):
                x2 = self._get_fact(self.matrix[i,j])
                if np.isclose(x,x2):
                    if self.matrix[i,j] < 0:
                        new[i,j] = -1
                    else:
                        new[i,j] = 1
                    continue
                else:
                    x = 0
                    new = []
                    break
            if len(new) == 0:
                break
        return x,new

    def _get_fact(self,x):
        if x.real != 0 and x.imag != 0:
            if np.isclose(x.real,x.imag):
                return x.real
            else:
                return 1
        if x.real !=0 and x.imag == 0:
            return x.real
        if x.real ==0 and x.imag != 0:
            return x.imag
        if x.real == 0 and x.imag == 0:
            return 0

    

class X(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [2,2]
        self.matrix = np.array([[0,1],[1,0]],dtype=complex)

class Y(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [2,2]
        self.matrix = np.array([[0,-self.i],[self.i,0]],dtype=complex)

class Z(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [2,2]
        self.matrix = np.array([[1,0],[0,-1]],dtype=complex)

class R(MatrixOperator):
    def __init__(self,axis,theta):
        super().__init__()
        self.axis = axis
        self.theta = theta
        self.shape = [2,2]
        self.matrix = np.eye(2,dtype=complex)*np.cos(0.5*self.theta) - pauli(self.axis)*self.i*np.sin(0.5*self.theta)

class I(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [2,2]
        self.matrix = np.eye(2,dtype=complex)

class H(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [2,2]
        self.matrix = np.sqrt(0.5)*np.array([[1,1],[1,-1]])

class CNOT(MatrixOperator):
    def __init__(self):
        super().__init__()
        self.shape = [4,4]
        self.matrix = np.eye(4)
        self.matrix[2:,2:] = not_eye(2)

if __name__ == '__main__':
    x = X()
    z = Z()
    print(z@x)
    print(x@z)
