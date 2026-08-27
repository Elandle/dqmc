import ifile

input = ifile.InputFile()

input.N           = 64
input.L           = 40
input.nstab       = 2
input.north       = 2
input.nbin        = 32
input.nmeassweep  = 64000
input.nskip       = 2
input.nequil      = 2000
input.dtau        = 0.08
input.U           = 4
input.mu          = 0
input.ckbfilename = "ckb.txt"
input.outfilename = "out.txt"
input.debfilename = "deb.txt"

Ls = [6, 12, 24, 48, 96, 120]
ds = ["000", "002", "005", "010", "020", "030"]
mus = ["p000", "p020", "p040", "p060", "p080", "p100", "p120", "p140", "p160", "p180", "p200",
       "n020", "n040", "n060", "n080", "n100", "n120", "n140", "n160", "n180", "n200"]
xbcs = ["pbcx"]
Us = [4]

for xbc in xbcs:
    for d in ds:
        for L in Ls:
            for mu in mus:
                for U in Us:
                    input.L = L
                    L = str(L)
                    muinput = mu
                    if muinput[0] == "p":
                        muinput = float(muinput[1:]) * pow(10, -2)
                    elif muinput[0] == "n":
                        muinput = -float(muinput[1:]) * pow(10, -2)
                    input.mu = muinput
                    input.U = U
                    U = str(U)
                    name =  "U" + U + "_" + "L" + L + "_" + "d" + d + "_" + "mu" + mu + "_x8_y8_" + xbc + "_pbcy"
                    ckbname = "d" + d + "_x8_y8_" + xbc + "_pbxy"
                    input.ckbfilename = ckbname + "_ckb.txt"
                    input.outfilename = name + "_out.txt"
                    input.debfilename = name + "_deb.txt"
                    input.print(name + "_input.txt")
                    print(name + "_input.txt")