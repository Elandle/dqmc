import equalgreens
import simulationsetup
import numpy
import bmult
import asvqrd
import scipy

def greens_Rup(S: simulationsetup.Simulation, i: int, l: int):
    deltaup = numpy.exp(-2 * S.alpha * S.h[i, l]) - 1
    Rup = 1 + (1 - S.Gup[i, i]) * deltaup
    return deltaup, Rup

def greens_Rdn(S: simulationsetup.Simulation, i: int, l: int):
    deltadn = numpy.exp(2 * S.alpha * S.h[i, l]) - 1
    Rdn = 1 + (1 - S.Gdn[i, i]) * deltadn
    return deltadn, Rdn

def upweight_definition(S: simulationsetup.Simulation):
    sigma = 1
    B = numpy.zeros((S.N, S.N, S.L))
    for l in range(S.L):
        B[:, :, l] = bmult.make_B(S, l, sigma)
    M = asvqrd.asvqrdmult(B)
    M = M + numpy.identity(S.N)
    Rup = scipy.linalg.det(M)
    return Rup