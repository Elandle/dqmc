import numpy
import ckb

# 2d Hatano Nelson
#   tp-->    tp-->    tp-->   tp-->
#   <--tm    <--tm    <--tm   <--tm
# 0 ------ 1 ------ 2 ------ 3
# |        |        |        |
# |v       |v       |v       |v
# |        |        |        |
# 4 ------ 5 ------ 6 ------ 7
#   tm-->    tm-->    tm-->   tm-->
#   <--tp    <--tp    <--tp   <--tp
#
# t(i, j) = hopping j to i
#

# only even m for now
n = 64
m = n//2

t = 1
v = 0.9
delta = 0.8
tp = t + delta
tm = t - delta

K = numpy.zeros((n, n))

for i in range(m-1):
    K[i, i+1], K[i+1, i] = tp, tm
K[0, m-1], K[m-1, 0] = tm, tp

for i in range(m, n-1):
    K[i, i+1], K[i+1, i] = tm, tp
K[m, n-1], K[n-1, m] = tp, tm

for i in range(m):
    K[i, i+m], K[i+m, i] = v, v

K = K.T

check = ckb.ckb(K)
check.saveckb("hatanockb08.txt")
check.savebipartite()