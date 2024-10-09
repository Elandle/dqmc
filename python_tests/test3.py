import numpy
import scipy

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
n = 2000
A = 1000 * numpy.random.rand(n, n)

Q, R, P = scipy.linalg.qr(A, pivoting=True)

print(scipy.linalg.norm(A - Q @ R @ invperm(P)))

