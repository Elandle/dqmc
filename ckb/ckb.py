import numpy
import scipy
import networkx

# Assumption: no diagonal, can easily be implemented or worked around
class ckb:
    def __init__(self, K, dtau=1):
        self.T = dtau * numpy.copy(K)
        self.n_T, _ = numpy.shape(self.T)
        self._colour()
        self.n = len(self.colours)
        self._constructckb()
        self.exact = scipy.linalg.expm(self.T)
        self._constructapprox()
        self.error = scipy.linalg.norm(self.exact - self.approx, ord=2)
        

    def _colour(self):
        self.G = networkx.from_numpy_array(self.T)
        self.LG = networkx.line_graph(self.G)
        self.colouring = networkx.coloring.greedy_color(self.LG)
        self.colours = set(self.colouring.values())
        self.mono_As = {colour: numpy.zeros((self.n_T, self.n_T), dtype=int) for colour in self.colours}
        for key, value in self.colouring.items():
            colour = value
            i, j = key
            self.mono_As[colour][i, j], self.mono_As[j, i] = 1, 1
        self.As = [self.T * mono_A for mono_A in self.mono_As.values()]

    def _constructckb(self):
        self.ckbcolours = self.n * [None]
        for k, colour in enumerate(self.colours):
            ijs = [ij for ij, col in self.colouring.items() if col == colour]
            pairs = [sympair(self.T[i, j], self.T[j, i], i, j) for i, j in ijs]
            self.ckbcolours[k] = ckbcolour(pairs)

    def _constructapprox(self):
        self.approx = numpy.identity(self.n_T)
        self.right_mult(self.approx)

    
    def right_mult(self, A):
        for colour in self.ckbcolours:
            colour.right_mult(A)

    def left_mult(self, A):
        for colour in reversed(self.ckbcolours):
            colour.left_mult(A)
    
    def saveckb(self, filename="ckb.txt"):
        with open(filename, "w") as file:
            print(self.n, file=file)
            for colour in self.colours:
                print(file=file)
                ijs = [ij for ij, col in self.colouring.items() if col == colour]
                print(len(ijs), file=file)
                for i, j in ijs:
                    print(i+1, j+1, self.T[i, j], self.T[j, i], file=file)

    def is_bipartite(self):
        return networkx.algorithms.is_bipartite(self.G)
    
    def savebipartite(self, filename="bip.txt", sgn=1):
        coloring = networkx.algorithms.bipartite.color(self.G)
        with open(filename, "w") as file:
            for vertex, color in coloring.items():
                vertex = vertex + 1
                if (color == 0):
                    color = sgn * 1
                else:
                    color = sgn * (-1)
                print(vertex, color, file=file)


class ckbcolour:
    def __init__(self, pairs):
        n = len(pairs)
        self.pairs = n * [None]
        for i, sym in enumerate(pairs):
            self.pairs[i] = sym

    def right_mult(self, A):
        for pair in self.pairs:
            pair.right_mult(A)

    def left_mult(self, A):
        for pair in reversed(self.pairs):
            pair.left_mult(A)


    

    

class sympair:
    def __init__(self, a, b, i, j):
        self.i = i
        self.j = j
        if (a * b > 0):
            x = numpy.sqrt(a * b)
            self.d  = numpy.cosh(x)
            self.ij = a * numpy.sinh(x) / x
            self.ji = b * numpy.sinh(x) / x
        elif (a * b < 0):
            x = numpy.sqrt(-a * b)
            self.d = numpy.cos(x)
            self.ij = a * numpy.sin(x) / x
            self.ji = b * numpy.sin(x) / x
        elif (b == 0):
            self.d = 1
            self.ij = a
            self.ji = 0
        elif (a == 0):
            self.d = 1
            self.ij = 0
            self.ji = b

    def right_mult(self, A):
        work = numpy.copy(A[:, self.i])
        A[:, self.i] = self.d * A[:, self.i] + self.ji * A[:, self.j]
        A[:, self.j] = self.d * A[:, self.j] + self.ij * work

    def left_mult(self, A):
        work = numpy.copy(A[self.i, :])
        A[self.i, :] = self.d * A[self.i, :] + self.ij * A[self.j, :]
        A[self.j, :] = self.d * A[self.j, :] + self.ji * work
