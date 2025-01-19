import numpy
import ckb

# 2d Hatano Nelson

n = 8

t = 1
v = 1
delta = 0.1
tp = t + delta
tm = t - delta

K = numpy.zeros((n, n))

K[0, 1], K[1, 0] = tp, tm
K[1, 2], K[2, 1] = tp, tm
K[2, 3], K[3, 2] = tp, tm
K[0, 3], K[3, 0] = tm, tp

K[4, 5], K[5, 4] = tm, tp
K[5, 6], K[6, 5] = tm, tp
K[6, 7], K[7, 6] = tm, tp
K[4, 7], K[7, 4] = tp, tm

K[0, 4], K[4, 0] = v, v
K[1, 5], K[5, 1] = v, v
K[2, 6], K[6, 2] = v, v
K[3, 7], K[7, 3] = v, v


check = ckb.ckb(K)
check.saveckb("hatanockb.txt")
