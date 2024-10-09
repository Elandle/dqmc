import simulationsetup
import equalgreens
import numpy
import scipy

def norm(A):
    return scipy.linalg.norm(A)

nx = 2
ny = 2
N = nx * ny
tx = 1
ty = 1

Kx = numpy.diag((nx-1) * [tx], k=1) + numpy.diag((nx-1) * [tx], k=-1)
Kx[0, nx-1] = 1; Kx[nx-1, 0] = 1
Ky = numpy.diag((ny-1) * [ty], k=1) + numpy.diag((ny-1) * [ty], k=-1)
Ky[0, ny-1] = 1; Ky[ny-1, 0] = 1
Idx = numpy.identity(nx)
Idy = numpy.identity(ny)
K = numpy.kron(Idy, Kx) + numpy.kron(Ky, Idx)

L = 6
U = 2
temp = 5

S = simulationsetup.Simulation(N, L, K, U, temp)

l = 1
Gup_exact = equalgreens.newGup_exact(S, l)
Gup_unstable_exact = equalgreens.newGup_unstable_exact(S, 1)

print(Gup_exact)
print(Gup_unstable_exact)

print(norm(Gup_unstable_exact - Gup_exact) / norm(Gup_unstable_exact))