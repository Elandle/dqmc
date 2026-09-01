import numpy
import ckb

# must both be even
nx = 8
ny = 8

t = 1.0
delta = 1.0

ckbname = f"2dhatano_{nx}x{ny}_d{delta:.3f}_ckb.txt"
bipname = f"2dhatano_{nx}x{ny}_d{delta:.3f}_bip.txt"


N = nx * ny

tp = t                
tm = t * (1.0 - delta)


def indx(x, y):
    return (y % ny) * nx + (x % nx)


def sublattice(x, y):
    # True: first sublattice
    # False: second sublattice
    return (x + y) % 2 == 0


K = numpy.zeros((N, N), dtype=float)

for y in range(ny):
    for x in range(nx):
        i = indx(x, y)

        j = indx(x + 1, y)
        if sublattice(x, y):
            K[i, j] = tp
            K[j, i] = tm
        else:
            K[j, i] = tp
            K[i, j] = tm

        j = indx(x, y + 1)
        if sublattice(x, y):
            K[j, i] = tp
            K[i, j] = tm
        else:
            K[i, j] = tp
            K[j, i] = tm

K = K.T

check = ckb.ckb(K)
check.saveckb(ckbname)
check.savebipartite(bipname)
