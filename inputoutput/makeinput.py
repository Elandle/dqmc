import ifile

input = ifile.InputFile()

input.N           = 64
input.L           = 40
input.nstab       = 2
input.north       = 2
input.nbin        = 32
input.nmeassweep  = 64000
input.nskip       = 2
input.nequil      = 200
input.dtau        = 0.08
input.U           = 4
input.mu          = 0
input.ckbfilename = "ckb.txt"
input.outfilename = "out.txt"
input.debfilename = "deb.txt"

Ls = [6, 12, 24, 48, 96, 120]
ds = ["000", "002", "005", "010"]
xbcs = ["obcx", "pbcx"]

for xbc in xbcs:
    for d in ds:
        for L in Ls:
            input.L = L
            L = str(L)
            basename = "d" + d + "_" + "x8_y8_" + xbc + "_pbcy"
            input.ckbfilename = basename + "_ckb.txt"
            input.outfilename = basename + "_out.txt"
            input.debfilename = basename + "_deb.txt"
            input.print(basename + "_input.txt")
            print(basename + "_input.txt")