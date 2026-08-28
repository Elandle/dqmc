import ifile

input = ifile.InputFile()

input.nstab = 10
input.nbin = 32
input.nmeassweep = 64000
input.nskip = 2
input.dtau = 0.08
input.mu = 0.0

sizes = [8, 10, 12]
Us = [1, 2, 3, 4]
Ls = [24, 48, 72, 102, 126, 150]
ds = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.10]
nequils = {8: 4000, 10: 6000, 12: 8000}

count = 0

for size in sizes:
    for U in Us:
        for L in Ls:
            for d in ds:

                input.summary = f"{size} x {size} delta = {d}"
                input.N = size * size
                input.L = L
                input.U = U
                input.north = 6
                input.nequil = nequils[size]

                beta = L * input.dtau

                dstr = str(int(round(d * 100))).zfill(3)
                sizestr = str(size)
                Ustr = str(U)
                Lstr = str(L)
                betastr = str(round(beta, 2))

                input.ckbfilename = "d" + dstr + "_x" + sizestr + "_y" + sizestr + "_pbcx_pbcy_ckb.txt"
                input.outfilename = "U" + Ustr + "_L" + Lstr + "_b" + betastr + "_d" + dstr + "_x" + sizestr + "_y" + sizestr + "_pbcx_pbcy_out.txt"
                input.debfilename = "U" + Ustr + "_L" + Lstr + "_b" + betastr + "_d" + dstr + "_x" + sizestr + "_y" + sizestr + "_pbcx_pbcy_deb.txt"
                input.bipfilename = "d" + dstr + "_x" + sizestr + "_y" + sizestr + "_pbcx_pbcy_bip.txt"

                inputfilename = "U" + Ustr + "_L" + Lstr + "_b" + betastr + "_d" + dstr + "_x" + sizestr + "_y" + sizestr + "_pbcx_pbcy_input.txt"

                input.print(inputfilename)
                print(inputfilename)

                count = count + 1

print(count, "total input files")