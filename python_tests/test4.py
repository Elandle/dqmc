import numpy
import scipy
import simulationsetup
import equalgreens

def norm(A):
    return scipy.linalg.norm(A)

l = 4
i = 12

K, N = simulationsetup.periodicsquare(5, 5, 1, 1)
S = simulationsetup.Simulation(N, 100, K, 4, 2, dtau=0.125)

S.Gup = equalgreens.newGup(S, l)
S.deltaup = equalgreens.getdeltaup(S, l, i)
S.Rup = equalgreens.getRup(S, l, i)

a = equalgreens.upflipupdate(S, i)
b = equalgreens.upflipdupatematrix(S, l, i)
c = equalgreens.upflipupdateloop(S, l, i)
S.h[i, l] = -S.h[i, l]
S.Gup = equalgreens.newGup_exact(S, l)
print(norm(a))
print(norm(b))
print(norm(c))
print(norm(S.Gup))
print(norm(a-b))
print(norm(a-c))
print(norm(b-c))
print(norm(a-S.Gup))
print(norm(b-S.Gup))
print(norm(c-S.Gup))
print(norm(a-S.Gup) / norm(S.Gup))
print(norm(b-S.Gup) / norm(S.Gup))
print(norm(c-S.Gup) / norm(S.Gup))