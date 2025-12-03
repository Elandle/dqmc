import numpy
import ckb

#
#      0  ----   1   ----   2  ----  3  --  ...  --    n-2  ----  n-1
#      |         |          |        |                  |          |
#  --  n  ----  n+1  ----  n+2 ---- n+3 --  ...  --   2*n-2 ---- 2*n-1 --
#      |         |          |        |                  |          |
#  -- 2*n ---- 2*n+1 ---- 2*n+2 ---- 2*n+3 -- ... --  3*n-2 ---- 3*n-2 --
#
#
#


m = 8
n = 8
N = m*n

tx = 1
deltax = 0.1
txp = tx + deltax
txm = tx - deltax

ty = 1
deltay = 0.0
typ = ty + deltay
tym = ty - deltay


K = numpy.zeros((N, N))

for i in range(m):
    for j in range(n-1):
        K[j+i*n, j+i*n+1] = txp

for i in range(m):
    for j in range(1, n):
        K[j+i*n, j+i*n-1] = txm

for i in range(m-1):
    for j in range(n):
        K[j+i*n, j+i*n+n] = tym

for i in range(1, m):
    for j in range(n):
        K[j+i*n, j+i*n-n] = typ

# Periodic y
for j in range(n):
    K[j, N-n+j] = typ
    K[N-n+j, j] = tym

# Periodic x
# for i in range(m):
#     K[i*n, i*n+n-1] = txm
#     K[i*n+n-1, i*n] = txp

K = K.T
check = ckb.ckb(K)
check.saveckb("d010_x8_y8_obcx_pbcy_ckb.txt")
check.savebipartite()
print(K)