import asvqrd
import numpy
import simulationsetup
import bmult
import scipy

nx = 10
ny = 10
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

L = 4
U = 2
temp = 5
dtau = 1

S = simulationsetup.Simulation(N, L, K, U, temp, dtau=dtau)

B = numpy.zeros((S.N, S.N, S.L))
Bexact = numpy.zeros((S.N, S.N, S.L))
Bckbexact = numpy.zeros((S.N, S.N, S.L))
for j in range(S.L):
    B[:, :, j] = bmult.make_B(S, j, 1)
    Bexact[:, :, j] = bmult.make_Bexact(S, j, 1)
    Bckbexact[:, :, j] = bmult.make_Bckbexact(S, j, 1)

X = asvqrd.asvqrdmult(B)
Xexact = asvqrd.asvqrdmult(Bexact)
Xckbexact = asvqrd.asvqrdmult(Bckbexact)

id = numpy.identity(S.N)
Y = numpy.copy(id)
Yexact = numpy.copy(id)
Yckbexact = numpy.copy(id)
for j in range(S.L):
    Y = B[:, :, j] @ Y
    Yexact = B[:, :, j] @ Yexact
    Yckbexact = B[:, :, j] @ Yckbexact

print(scipy.linalg.norm(X - Y))
print(scipy.linalg.norm(Xexact - Yexact))
print(scipy.linalg.norm(Xckbexact - Yckbexact))
print(scipy.linalg.norm(X - Xexact))
print(scipy.linalg.norm(X - Xckbexact))
print(scipy.linalg.norm(Xckbexact - Xexact))

print(scipy.linalg.norm(S.ckb.right_mult(id) - S.ckb.left_mult(id)))
print(scipy.linalg.norm(S.expT - S.ckb.left_mult(id)) / scipy.linalg.norm(S.expT))
print(scipy.linalg.norm(S.expT - S.ckb.right_mult(id)) / scipy.linalg.norm(S.expT))
