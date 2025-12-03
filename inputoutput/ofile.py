import linecache
import numpy

class OutputFile():
    def __init__(self, ofile):
        self.ofile         = None
        self.N             = None
        self.sign          = None; self.signerr          = None
        self.total_density = None; self.total_densityerr = None
        self.upden         = None; self.updenerr         = None
        self.dnden         = None; self.dndenerr         = None
        self.KE            = None; self.KEerr            = None
        self.CHEE          = None; self.CHEEerr          = None
        self.PE            = None; self.PEerr            = None
        self.E             = None; self.Eerr             = None
        self.AFSF          = None; self.AFSFerr          = None
        self.Gup           = None; self.Guperr           = None
        self.Gdn           = None; self.Gdnerr           = None
        self.SDcorr        = None; self.SDcorrerr        = None
        self.SScorr        = None; self.SScorrerr        = None
        self.upden_full    = None; self.upden_fullerr    = None
        self.dnden_full    = None; self.dnden_fullerr    = None
        self.dbocc_full    = None; self.dbocc_fullerr    = None
        self.magmom_full   = None; self.magmom_fullerr   = None
        self.read_ofile(ofile)

    def find_similar_line(self, lines, line, linenum=False, startnum=0):
        for i, l in enumerate(lines[startnum:]):
            if line in l:
                if linenum:
                    return l, i+startnum
                else:
                    return l
    
    def scalar_avg_err(self, line):
        _, _, right = line.partition("=")
        left, _, right = right.partition("+-")
        return float(left), float(right)
    
    def matrix_avg_err(self, lines, linenum, nrows):
        matrix = numpy.array([lines[linenum+i].split() for i in range(nrows)], dtype=numpy.float64)
        matrixerr = numpy.array([lines[linenum+nrows+1+i].split() for i in range(nrows)], dtype=numpy.float64)
        return matrix, matrixerr


    def read_ofile(self, ofile):
        self.ofile = ofile
        file = open(ofile, "r")
        lines = file.read().splitlines()
        file.close()

        line = self.find_similar_line(lines, "sign")
        self.sign, self.signerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "total density")
        self.total_density, self.total_densityerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "upden")
        self.upden, self.updenerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "dnden")
        self.dnden, self.dndenerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "KE")
        self.KE, self.KEerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "CHEE")
        self.CHEE, self.CHEEerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "PE")
        self.PE, self.PEerr = self.scalar_avg_err(line)

        line = self.find_similar_line(lines, "E")
        self.E, self.Eerr = self.scalar_avg_err(line)
        
        line = self.find_similar_line(lines, "antiferromagnetic")
        self.AFSF, self.AFSFerr = self.scalar_avg_err(line)

        _, linenum1 = self.find_similar_line(lines, "Gup", True)
        _, linenum2 = self.find_similar_line(lines, "+-", True, linenum1)
        self.N = linenum2 - linenum1 - 1
        self.Gup, self.Guperr = self.matrix_avg_err(lines, linenum1+1, self.N)

        _, linenum = self.find_similar_line(lines, "Gdn", True)
        self.Gdn, self.Gdnerr = self.matrix_avg_err(lines, linenum+1, self.N)

        _, linenum = self.find_similar_line(lines, "Spin density", True)
        self.SDcorr, self.SDcorrerr = self.matrix_avg_err(lines, linenum+1, self.N)

        _, linenum = self.find_similar_line(lines, "Spin spin", True)
        self.SScorr, self.SScorrerr = self.matrix_avg_err(lines, linenum+1, self.N)


        _, linenum = self.find_similar_line(lines, "Average upden (full)", True)
        self.upden_full, self.upden_fullerr = self.matrix_avg_err(lines, linenum+1, 1)

        _, linenum = self.find_similar_line(lines, "Average dnden (full)", True)
        self.dnden_full, self.dnden_fullerr = self.matrix_avg_err(lines, linenum+1, 1)

        _, linenum = self.find_similar_line(lines, "Average double occupancy (full)", True)
        self.dbocc_full, self.dbocc_fullerr = self.matrix_avg_err(lines, linenum+1, 1)

        _, linenum = self.find_similar_line(lines, "Average magnetic moment (full)", True)
        self.magmom_full, self.magmom_fullerr = self.matrix_avg_err(lines, linenum+1, 1)
