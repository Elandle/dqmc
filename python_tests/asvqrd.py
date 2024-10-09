import numpy
import scipy

# Python implementation of the ASvQRD algorithm without
# care about speed or memory - for testing

def perm(P: numpy.ndarray[int]):
    n = len(P)
    A = numpy.zeros((n, n))
    for i in range(n):
        j = P[i]
        A[j, i] = 1
    return A

def invperm(P: numpy.ndarray[int]):
    return perm(P).T

# A @ perm(P) = Q @ R
# A           = Q @ R @ invperm(P)
# Q, R, P = scipy.linalg.qr(A, pivoting=True)

def firstchainiteration(B):
    n, _ = numpy.shape(B)
    # B = Q @ R @ invperm(P)
    Q, R, P = scipy.linalg.qr(B, pivoting=True)
    D = numpy.diag(numpy.diag(R))
    T = scipy.linalg.inv(D) @ R @ invperm(P)
    return Q, D, T

def chainiteration(B, Q, D, T):
    n, _ = numpy.shape(B)
    C = (B @ Q) @ D
    Q, R, P = scipy.linalg.qr(C, pivoting=True)
    D = numpy.diag(numpy.diag(R))
    T = (scipy.linalg.inv(D) @ R @ invperm(P)) @ T
    return Q, D, T

def DbDs(D):
    n, _ = numpy.shape(D)
    Db = numpy.zeros((n, n))
    Ds = numpy.zeros((n, n))
    for i in range(n):
        if (abs(D[i, i]) > 1):
            Db[i, i] = D[i, i]
            Ds[i, i] = 1
        else:
            Db[i, i] = 1
            Ds[i, i] = D[i, i]
    return Db, Ds

# B: shape (n, n, L)
# Returns B[:, :, L-1] * ... * B[:, :, 0]
def asvqrdmult(B):
    n, _, L = numpy.shape(B)

    # First iteration
    j = 0
    Bj = B[:, :, j]
    Q, D, T = firstchainiteration(Bj)
    for j in range(1, L):
        Bj = B[:, :, j]
        Q, D, T = chainiteration(Bj, Q, D, T)
    Db, Ds = DbDs(D)
    return (Q @ Db @ Ds) @ T


