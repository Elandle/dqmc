import numpy
import checkerboard
import scipy

class Simulation:
    def __init__(self, N, L, K, U, temp, mu=0, dtau=0.125, h=None, nstab=4, nbin=32, nmeassweep=32*1000,
                 nskip=5, nequil=1000):
        self.N = N
        self.L = L
        self.U = U
        self.temp = temp
        self.mu = mu
        self.dtau = dtau
        self.nstab = nstab
        self.nbin = nbin
        self.nmeassweep = nmeassweep
        self.nskip = nskip
        self.nequil = nequil

        self.binsize = nmeassweep / nbin
        self.ntotal = nequil + (nmeassweep - 1 ) * nskip + nmeassweep
        self.beta = 1 / temp
        self.alpha = numpy.acosh(numpy.exp(dtau * U / 2))
        self.aldmu = self.alpha * numpy.exp(dtau * U / 2)
        self.aldmuinv = 1 / self.aldmu

        if h is None:
            self.h = numpy.random.choice([-1, 1], (N, L))
        else:
            self.h = numpy.copy(h)

        self.Gup = numpy.zeros((N, N))
        self.Gdn = numpy.zeros((N, N))
        self.R = 0
        self.deltaup = 0
        self.deltadn = 0
        self.upstabi = 0
        self.dnstabi = 0
        self.upsgn = 0
        self.dnsgn = 0
        self.sgn = 0
        self.ckb = checkerboard.checkerboard(K, dtau)
        self.ckbinv = checkerboard.checkerboard(K, -dtau)
        self.T = dtau * numpy.copy(K)
        self.expT = scipy.linalg.expm(self.T)
        self.expTinv = scipy.linalg.inv(self.expT)
        self.expTckb = self.ckb.approx
        self.expTckbinv = scipy.linalg.inv(self.expTckb)

def periodicsquare(nx, ny, tx, ty):
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
    return K, N
