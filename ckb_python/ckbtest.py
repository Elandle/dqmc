import numpy
import networkx
import scipy
import matplotlib.pyplot
import ckb

# 4 honeycombs
#
#     15 - - - - 16
#     /           \
#    /             \
#   /               \
# 12                 13 - - - - 14
#  \                /             \
#   \              /               \
#    \            /                 \
#     9 - - - - 10                  11     
#    /           \                 /                  1 + delta
#   /             \               /           +     +
#  /               \             /             \   /
# 6                 7 - - - - - 8               \ / 
#  \               /            \         -  - -   - - +
#   \             /              \              /  \
#    \           /                \            /    \
#     3 - - - - 4                  5          -      -
#                \                /       1 - delta
#                 \              /
#                  \            /
#                   1 - - - - - 2
#

delta = 0.5
plus = 1 + delta
minus = 1 - delta
n = 16
H = numpy.zeros((n, n))

H[1-1, 2-1], H[1-1, 4-1] = plus, plus
H[2-1, 1-1], H[2-1, 5-1] = minus, plus
H[3-1, 4-1], H[3-1, 6-1] = plus, plus
H[4-1, 1-1], H[4-1, 3-1], H[4-1, 7-1] = minus, minus, plus
H[5-1, 2-1], H[5-1, 8-1] = minus, plus
H[6-1, 3-1], H[6-1, 9-1] = minus, plus
H[7-1, 4-1], H[7-1, 8-1], H[7-1, 10-1] = minus, plus, plus
H[8-1, 5-1], H[8-1, 7-1], H[8-1, 11-1] = minus, minus, plus
H[9-1, 6-1], H[9-1, 10-1], H[9-1, 12-1] = minus, plus, plus
H[10-1, 7-1], H[10-1, 9-1], H[10-1, 13-1] = minus, minus, plus
H[11-1, 8-1], H[11-1, 14-1] = minus, plus
H[12-1, 9-1], H[12-1, 15-1] = minus, plus
H[13-1, 10-1], H[13-1, 14-1], H[13-1, 16-1] = minus, plus, plus
H[14-1, 11-1], H[14-1, 13-1] = minus, minus
H[15-1, 12-1], H[15-1, 16-1] = minus, plus
H[16-1, 13-1], H[16-1, 15-1] = minus, minus


check = ckb.Ckb(H)

# Python print version
#for color in check.colors:
#    print(color)
#    ijs = [ij for ij, col in check.coloring.items() if col == color]
#    for i, j in ijs:
#        print(i, j, H[i, j], H[j, i])
#    print()




print(len(check.colors))
for color in check.colors:
    print()
    ijs = [ij for ij, col in check.coloring.items() if col == color]
    print(len(ijs))
    for i, j in ijs:
        print(i+1, j+1, H[i, j], H[j, i])