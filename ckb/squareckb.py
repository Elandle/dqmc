import numpy
import ckb

# 2d PERL
# Periodic rectangular lattice
# nx sites in the horizontal direction, and ny sites in the vertical direction
# n = nx * ny total sites
# tx hopping in the horizontal direction, and ty hopping in the vertical direction

nx = 2
ny = 2

n = nx * ny

tx = 1
ty = 1

Kx = numpy.diag((nx-1) * [tx], k=1) + numpy.diag((nx-1) * [tx], k=-1)
Kx[0, nx-1] = tx; Kx[nx-1, 0] = tx

Ky = numpy.diag((ny-1) * [ty], k=1) + numpy.diag((ny-1) * [ty], k=-1)
Ky[0, ny-1] = ty; Ky[ny-1, 0] = ty

Idx = numpy.identity(nx)
Idy = numpy.identity(ny)

K = numpy.kron(Idy, Kx) + numpy.kron(Ky, Idx)


check = ckb.ckb(K)
check.saveckb("squareckbtwo.txt")
