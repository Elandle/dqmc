import numpy
import ofile
import matplotlib.pyplot
matplotlib.rcParams["figure.dpi"] = 400

deltas = ["00", "01"]
vs     = ["050", "075", "100", "125", "150", "175", "200"]
Ls     = ["50", "76", "100"]

outputs = dict()
for d in deltas:
    for v in vs:
        for L in Ls:
            name = "d" + d + "v" + v + "L" + L + ".txt"
            outputs[name] = ofile.OutputFile(name)

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

def nameparameters(name):
    delta = betweentext(name, "d", "v")
    v = betweentext(name, "v", "L")
    L = betweentext(name, "L", ".")
    delta = todecimal(delta, 1)
    v = todecimal(v, 1)
    L = int(L)
    return delta, v, L

outputs = {nameparameters(key): val for key, val in outputs.items()}
dtau = 0.08
colors = ["red", "blue"]
n = int(next(iter(outputs.values())).N / 2)


for L in [50, 76, 100]:
    figure, (axisSSlr, axisSSud) = matplotlib.pyplot.subplots(2, 1, sharex=True, figsize=(5, 9), constrained_layout=True)

    for i, delta in enumerate([0.0, 0.1]):
        vs      = []
        SSlr    = []
        SSlrerr = []
        SSud    = []
        SSuderr = []

        for (deltakey, vkey, Lkey), file in outputs.items():
            if (deltakey == delta and Lkey == L):
                vs.append(vkey)
                SSlr.append(file.SScorr[0, 1])
                SSlrerr.append(file.SScorrerr[0, 1])
                SSud.append(file.SScorr[0, n])
                SSuderr.append(file.SScorrerr[0, n])

        axisSSlr.errorbar(vs, SSlr, yerr=SSlrerr, fmt="*", capsize=1.5, elinewidth=1.5, markersize=2.5, label=f"δ = {deltakey}", color=colors[i])
        axisSSud.errorbar(vs, SSud, yerr=SSuderr, fmt="*", capsize=1.5, elinewidth=1.5, markersize=2.5, label=f"δ = {deltakey}", color=colors[i])


    axisSSlr.set_ylabel("SS (left-right)", fontsize=12)
    axisSSlr.tick_params(axis="both", which="major", labelsize=12)
    axisSSlr.legend()
    axisSSlr.set_aspect("equal", adjustable="box")
    axisSSlr.grid(True)
    axisSSud.set_ylabel("SS (up-down)", fontsize=12)
    axisSSud.tick_params(axis="both", which="major", labelsize=12)
    axisSSud.legend()
    axisSSud.set_aspect("equal", adjustable="box")
    axisSSud.grid(True)
    figure.suptitle(f"β = {L*dtau} left-right SS correlation (24 sites)", fontsize=14, y=0.98)
    figure.supxlabel("v", fontsize=12)
    # matplotlib.pyplot.show()

figure.savefig("abc.png")

