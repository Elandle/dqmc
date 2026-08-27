import numpy
import ofile
import matplotlib.pyplot
import fnmatch
import glob
import contextlib
matplotlib.rcParams["figure.dpi"] = 400

deltas = ["000", "002", "005", "010", ]
Ls     = ["6", "12", "24", "48", "96", "120"]
xbcs   = ["obcx", "pbcx"]

#outputs = dict()
#for d in deltas:
#    for L in Ls:
#        for xbc in xbcs:
#            name = "d" + d + "_L" + L + "_x8_y8_" + xbc + "_pbcy_out" + ".txt"
#            outputs[name] = ofile.OutputFile(name)




def lrindices(name, left, right):
    indl = name.find(left)
    indr = name.find(right)
    return indl, indr

def betweenindices(name, indl, indr):
    return name[indl+1:indr]

def betweentext(name, left, right):
    indl, indr = lrindices(name, left, right)
    return betweenindices(name, indl, indr)

def todecimal(name, point):
    name = str(name)
    if point == 0:
        return float(name)
    return float(name[:point] + "." + name[point:])


outputs = {}
for name in glob.glob("*_out.txt"):
    outputs[name] = ofile.OutputFile(name)

dtau = 0.08

for key, val in outputs.items():
    delta = betweentext(key, "d", "_")
    delta = todecimal(delta, 1)
    L = betweentext(key, "L", "_x")
    L = int(L)
    xbc = betweentext(key, "y8_", "_pbcy")[2:]
    if (xbc == "pbcx"):
        xbc = True
    else:
        xbc = False
    val.delta = delta
    val.L = L
    val.pbcx = xbc
    val.beta = L * dtau
    val.temp = 1 / val.beta
    
def printthings_delta_pbcx(outputs, delta, pbcx, nx, ny):
    n = nx*ny
    print(f"delta = {delta}, pbcx = {pbcx}")
    print("beta1 KE2 KEerr3 PE4 PEerr5 CHEE6 CHEEerr7 E8 Eerr9 sign10 signerr11 totaldensity12 totaldensityerr13 upden14 updenerr15 dnden16 dndenerr17 SScorrlr18 SScorrlrerr19 SScorrself20 SScorrselferr21 SScorrnext22 SScorrnexterr23 dbocc24 dboccerr25")
    for key, val in outputs.items():
        if val.delta == delta:
            if val.pbcx is pbcx:
                print(f"{val.beta:.5f}", end=" ")
                print(f"{val.KE/n:.5f} {val.KEerr/n:.5f}", end=" ")
                print(f"{val.PE/n:.5f} {val.PEerr/n:.5f}", end=" ")
                print(f"{val.CHEE/n:.5f} {val.CHEEerr/n:.5f}", end=" ")
                print(f"{val.E/n:.5f} {val.Eerr/n:.5f}", end=" ")
                print(f"{val.sign:.5f} {val.signerr:.5f}", end=" ")
                print(f"{val.total_density:.5f} {val.total_densityerr:.5f}", end=" ")
                print(f"{val.upden:.5f} {val.updenerr:.5f}", end=" ")
                print(f"{val.dnden:.5f} {val.dndenerr:.5f}", end=" ")
                print(f"{val.SScorr[0, 1]:.5f} {val.SScorrerr[0, 1]:.5f}", end=" ")
                print(f"{val.SScorr[0, 0]:.5f} {val.SScorrerr[0, 0]:.5f}", end=" ")
                print(f"{val.SScorr[0, 2]:.5f} {val.SScorrerr[0, 2]:.5f}", end=" ")
                v = 1 / val.dbocc_fullerr**2
                vsum = numpy.sum(v)
                dbocc = numpy.sum(v*val.dbocc_full) / vsum
                dboccerr = numpy.sqrt(1/vsum)
                print(f"{dbocc:.5f} {dboccerr:.5f}", end=" ")
                print("")


nx = 8 ; ny = 8
delta = 0.0 ; pbcx = True

pbcs = [("pbcx", True), ("obcx", False)]
deltas = [("d000", 0.0), ("d002", 0.02), ("d005", 0.05), ("d010", 0.10)]

for pbcname, pbcx in pbcs:
    for deltaname, delta in deltas:
        fname = deltaname + "_" + pbcname + "_data.txt"
        with open(fname, "w") as f:
            with contextlib.redirect_stdout(f):
                printthings_delta_pbcx(outputs, delta, pbcx, nx, ny)