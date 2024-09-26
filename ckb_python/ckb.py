import scipy
import numpy
import networkx

class Ckb:
    def __init__(self, H):
        self.n, _ = numpy.shape(H)
        self._D = numpy.diag(H)
        self.D = numpy.zeros((self.n,self.n))
        numpy.fill_diagonal(self.D, numpy.exp(self._D))
        self._H = numpy.copy(H)
        self.H = numpy.copy(H)
        numpy.fill_diagonal(self.H, 0)
        self._color()
        self._expAs()
        self._asmatrix()
    
    def _color(self):
        # strategy: to color edges, make a new graph with edges as vertices
        # and vertex color that graph
        self.G = networkx.from_numpy_array(self.H)
        self.LG = networkx.line_graph(self.G)
        self.coloring = networkx.coloring.greedy_color(self.LG)
        self.colors = set(self.coloring.values())
        self.mono_As = {color: numpy.zeros((self.n,self.n), dtype=int) for color in self.colors}
        for key, value in self.coloring.items():
            color = value
            i, j = key
            self.mono_As[color][i, j], self.mono_As[color][j, i] = 1, 1
        self.As = [self.H * mono_A for mono_A in self.mono_As.values()]

    def _asmatrix(self):
        self.M = numpy.identity(self.n)
        for mono_A in self.expAs:
            self.M = mono_A @ self.M
        self.M = self.D @ self.M

    def _expAs(self):
        self.expAs = [scipy.linalg.expm(A) for A in self.As]

    def multiply(self, matrix):
        M = numpy.copy(matrix)
        for A in self.expAs:
            M = A @ M
        M = self.D @ M
        return M
    
    def exacterror(self):
        return scipy.linalg.norm(self.M - scipy.linalg.expm(self._H), 2)
    
    def saveckb(self):
        error = self.exacterror()
        ncolors = len(self.expAs)
        n = self.n
        print(n)
        print(ncolors)
        print(error)
        for k, expA in enumerate(self.expAs):
            print()
            I, J = expA.nonzero()
            for i, j in zip(I, J):
                print(i, j, expA[i, j])