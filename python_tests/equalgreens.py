import numpy
import scipy
import bmult
import asvqrd
import simulationsetup
import typing



def upflipupdate(S, i):
    x = S.deltaup / S.Rup
    ei = numpy.zeros(S.N)
    ei[i] = 1
    M = S.Gup + x * numpy.outer(S.Gup[:, i] - ei, S.Gup[i, :])
    return M

def dnflipupdate(S, i):
    x = S.deltaup / S.Rup
    ei = numpy.zeros(S.N)
    ei[i] = 1
    M = S.Gdn + x * numpy.outer(S.Gdn[:, i] - ei, S.Gdn[i, :])
    return M

def getdeltaup(S, l, i):
    sigma = 1
    deltaup = deltaup = numpy.exp(sigma * S.alpha * (-S.h[i, l] - S.h[i, l])) - 1
    return deltaup

def getRup(S, l, i):
    deltaup = getdeltaup(S, l, i)
    Rup = 1 + deltaup * (1 - S.Gup[i, i])
    return Rup


def upflipupdateloop(S, l, i):
    sigma = 1
    deltaup = getdeltaup(S, l, i)
    Rup = getRup(S, l, i)
    X = numpy.identity(S.N) - S.Gup
    M = numpy.zeros((S.N, S.N))
    for j in range(S.N):
        for k in range(S.N):
            M[j, k] = S.Gup[j, k] - (1 / Rup) * S.Gup[j, i] * deltaup * X[i, k]
    return M

def dnflipupdateloop(S, l, i):
    sigma = -1
    deltadn = numpy.exp(sigma * S.alpha * (-S.h[i, l] - S.h[i, l])) - 1
    Rdn = 1 + deltadn * (1 - S.Gdn[i, i])
    X = numpy.identity(S.N) - S.Gdn
    M = numpy.zeros((S.N, S.N))
    for j in range(S.N):
        for k in range(S.N):
            M[j, k] = S.Gdn[j, k] - (1 / Rdn) * deltadn * X[i, k]
    return M



def upflipdupatematrix(S, l, i):
    sigma = 1
    deltaup = numpy.exp(sigma * S.alpha * (-S.h[i, l] - S.h[i, l])) - 1
    Rup = 1 + deltaup * (1 - S.Gup[i, i])
    delta = numpy.zeros((S.N, S.N))
    delta[i, i] = deltaup
    M = S.Gup - (1 / Rup) * S.Gup @ delta @ (numpy.identity(S.N) - S.Gup)
    return M

def dnflipdupatematrix(S, l, i):
    sigma = -1
    deltadn = numpy.exp(sigma * S.alpha * (-S.h[i, l] - S.h[i, l])) - 1
    Rdn = 1 + deltadn * (1 - S.Gdn[i, i])
    delta = numpy.zeros((S.N, S.N))
    delta[i, i] = deltadn
    M = S.Gdn - (1 / Rdn) * S.Gdn @ delta @ (numpy.identity(S.N) - S.Gdn)
    return M



def upwrap(S, l):
    M = bmult.right_Binvmult(S, S.Gup, l, 1)
    M = bmult.left_Bmult(S, M, l, 1)
    return M

def dnwrap(S, l):
    M = bmult.right_Binvmult(S, S.Gdn, l, 1)
    M = bmult.left_Bmult(S, M, l, 1)
    return M



def upwrap_exact(S, l):
    M = bmult.right_Binvmult_exact(S, S.Gup, l, 1)
    M = bmult.left_Bmult_exact(S, M, l, 1)
    return M

def dnwrap_exact(S, l):
    M = bmult.right_Binvmult_exact(S, S.Gdn, l, 1)
    M = bmult.left_Bmult_exact(S, M, l, 1)
    return M



def upwrap_ckbexact(S, l):
    M = bmult.right_Binvmult_ckbexact(S, S.Gup, l, 1)
    M = bmult.left_Bmult_ckbexact(S, M, l, 1)
    return M

def dnwrap_ckbexact(S, l):
    M = bmult.right_Binvmult_ckbexact(S, S.Gdn, l, 1)
    M = bmult.left_Bmult_ckbexact(S, M, l, 1)
    return M



def newGupchain(S: simulationsetup.Simulation, l: int):
    sigma = 1
    j = 0
    Bj = bmult.make_B(S, j, sigma)
    Q, D, T = asvqrd.firstchainiteration(Bj)
    for j in range(1, S.L):
        Bj = bmult.make_B(S, getj(j, S.L, l), sigma)
        Q, D, T = asvqrd.chainiteration(Bj, Q, D, T)
    return Q, D, T

def newGup(S, l):
    Q, D, T = newGupchain(S, l)
    Db, Ds = asvqrd.DbDs(D)
    Dbinv = scipy.linalg.inv(Db)
    M = scipy.linalg.inv(Dbinv @ Q.T + Ds @ T) @ Dbinv @ Q.T
    return M



def newGupchain_exact(S: simulationsetup.Simulation, l: int) -> tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    sigma = 1
    j = 0
    print(getj(j, S.L, l))
    Bj = bmult.make_Bexact(S, j, sigma)
    Q, D, T = asvqrd.firstchainiteration(Bj)
    for j in range(1, S.L):
        print(getj(j, S.L, l))
        Bj = bmult.make_Bexact(S, getj(j, S.L, l), sigma)
        Q, D, T = asvqrd.chainiteration(Bj, Q, D, T)
    return Q, D, T

def newGup_exact(S, l):
    Q, D, T = newGupchain_exact(S, l)
    Db, Ds = asvqrd.DbDs(D)
    Dbinv = scipy.linalg.inv(Db)
    M = scipy.linalg.inv(Dbinv @ Q.T + Ds @ T) @ Dbinv @ Q.T
    #M = scipy.linalg.inv(numpy.identity(S.N) + Q @ D @ T)
    return M



def newGupchain_ckbexact(S, l):
    sigma = 1
    j = 0
    Bj = bmult.make_Bckbexact(S, j, sigma)
    Q, D, T = asvqrd.firstchainiteration(Bj)
    for j in range(1, S.L):
        Bj = bmult.make_Bexact(S, getj(j, S.L, l), sigma)
        Q, D, T = asvqrd.chainiteration(Bj, Q, D, T)
    return Q, D, T

def newGup_ckbexact(S, l):
    Q, D, T = newGupchain_ckbexact(S, l)
    Db, Ds = asvqrd.DbDs(D)
    Dbinv = scipy.linalg.inv(Db)
    M = scipy.linalg.inv(Dbinv @ Q.T + Ds @ T) @ Dbinv @ Q.T
    return M



def newGup_unstable(S, l):
    sigma = 1
    M = numpy.identity(S.N)
    for j in range(l, S.L):
        Bj = bmult.make_B(S, j, sigma)
        M = Bj @ M
    for j in range(l):
        Bj = bmult.make_B(S, j, sigma)
        M = Bj @ M
    M = numpy.identity(S.N) + M
    M = scipy.linalg.inv(M)
    return M

def newGup_unstable_exact(S, l):
    sigma = 1
    M = numpy.identity(S.N)
    for j in range(l, S.L):
        print(j)
        Bj = bmult.make_Bexact(S, j, sigma)
        M = Bj @ M
    for j in range(l):
        print(j)
        Bj = bmult.make_Bexact(S, j, sigma)
        M = Bj @ M
    M = numpy.identity(S.N) + M
    M = scipy.linalg.inv(M)
    return M

def newGup_unstable_ckbexact(S, l):
    sigma = 1
    M = numpy.identity(S.N)
    for j in range(l, S.L):
        Bj = bmult.make_Bckbexact(S, j, sigma)
        M = Bj @ M
    for j in range(l):
        Bj = bmult.make_Bckbexact(S, j, sigma)
        M = Bj @ M
    M = numpy.identity(S.N) + M
    M = scipy.linalg.inv(M)
    return M



# Fortran:
# B(l) * ... * B(1) * B(L) * ... * B(l+1)
# Python:
# B(l-1) * ... * B(0) * B(L-1) * ... * B(l)

def getj(i, N, l):
    if (i < N - l):
        j = l + i
    else:
        j = l - N + i
    return j