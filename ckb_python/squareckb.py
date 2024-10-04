import numpy
import networkx
import scipy
import ckb

# 2d PERL
# Periodic rectangular lattice
# nx sites in the horizontal direction, and ny sites in the vertical direction
# n = nx * ny total sites
# tx hopping in the horizontal direction, and ty hopping in the vertical direction

nx = 10
ny = 10

n = nx * ny

tx = 1
ty = 1

Kx = numpy.diag((nx-1) * [tx], k=1) + numpy.diag((nx-1) * [tx], k=-1)
Kx[0, nx-1] = 1; Kx[nx-1, 0] = 1

Ky = numpy.diag((ny-1) * [ty], k=1) + numpy.diag((ny-1) * [ty], k=-1)
Ky[0, ny-1] = 1; Ky[ny-1, 0] = 1

Idx = numpy.identity(nx)
Idy = numpy.identity(ny)

K = numpy.kron(Idy, Kx) + numpy.kron(Ky, Idx)


check = ckb.Ckb(K)
print(len(check.colors))
for color in check.colors:
    print()
    ijs = [ij for ij, col in check.coloring.items() if col == color]
    print(len(ijs))
    for i, j in ijs:
        print(i+1, j+1, K[i, j], K[j, i])