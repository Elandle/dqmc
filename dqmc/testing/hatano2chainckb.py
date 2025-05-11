import numpy
import ckb
import scipy

delta = 0.4
t = 1
v = 1

ckbfilename = "delta00.txt"


n = 8
tp = t + delta
tm = t - delta

K = numpy.zeros((n, n))

K[0, 1], K[1, 2], K[2, 3], K[3, 0] = tp, tp, tp, tp
K[1, 0], K[2, 1], K[3, 2], K[0, 3] = tm, tm, tm, tm

K[0, 4], K[1, 5], K[2, 6], K[3, 7] = v, v, v, v
K[4, 0], K[5, 1], K[6, 2], K[7, 3] = v, v, v, v

K[4, 5], K[5, 6], K[6, 7], K[7, 4] = tm, tm, tm, tm
K[5, 4], K[6, 5], K[7, 6], K[4, 7] = tp, tp, tp, tp


# K[0, 4], K[1, 5], K[2, 6], K[3, 7] = t, t, t, t
# K[4, 0], K[5, 1], K[6, 2], K[7, 3] = t, t, t, t



#check = ckb.ckb(K)
#check.saveckb(ckbfilename)
beta = 1 * 0.08
M = numpy.identity(n) + scipy.linalg.expm(-beta * K)
M = scipy.linalg.inv(M)


eigenvalues, eigenvectors = numpy.linalg.eig(M)
for a in eigenvalues:
    print(a)