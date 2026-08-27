class InputFile():
    def __init__(self):
        self.N           = None
        self.L           = None
        self.nstab       = None
        self.north       = None
        self.nbin        = None
        self.nmeassweep  = None
        self.nskip       = None
        self.nequil      = None
        self.dtau        = None
        self.U           = None
        self.mu          = None
        self.ckbfilename = None
        self.outfilename = None
        self.debfilename = None
    def print(self, ofile=None, force=False):
        allset = self.checkset()
        if (not allset):
            print("Not all parameters are set for printing")
            if (not force):
                print("Stopping printing")
                return
            else:
                print("Printing has been forced. Proceeding to print")
        self.readyparameters()
        if ofile is not None:
            file = open(ofile, "w")
        else:
            file = None
        print(f"N           = {self.N}"          , file=file)
        print(f"L           = {self.L}"          , file=file)
        print(f"nstab       = {self.nstab}"      , file=file)
        print(f"north       = {self.north}"      , file=file)
        print(f"nbin        = {self.nbin}"       , file=file)
        print(f"nmeassweep  = {self.nmeassweep}" , file=file)
        print(f"nskip       = {self.nskip}"      , file=file)
        print(f"nequil      = {self.nequil}"     , file=file)
        print(f"dtau        = {self.dtau}"       , file=file)
        print(f"U           = {self.U}"          , file=file)
        print(f"mu          = {self.mu}"         , file=file)
        print(f"ckbfilename = {self.ckbfilename}", file=file)
        print(f"outfilename = {self.outfilename}", file=file)
        print(f"debfilename = {self.debfilename}", file=file)
        if ofile is not None:
            file.close()
    def readyparameters(self):
        self.N           = int(self.N)           if (self.N is not None)           else None
        self.L           = int(self.L)           if (self.L is not None)           else None
        self.nstab       = int(self.nstab)       if (self.nstab is not None)       else None
        self.north       = int(self.north)       if (self.north is not None)       else None
        self.nbin        = int(self.nbin)        if (self.nbin is not None)        else None
        self.nmeassweep  = int(self.nmeassweep)  if (self.nmeassweep is not None)  else None
        self.nskip       = int(self.nskip)       if (self.nskip is not None)       else None
        self.nequil      = int(self.nequil)      if (self.nequil is not None)      else None
        self.dtau        = float(self.dtau)      if (self.dtau is not None)        else None
        self.U           = float(self.U)         if (self.U is not None)           else None
        self.mu          = float(self.mu)        if (self.mu is not None)          else None
        self.ckbfilename = str(self.ckbfilename) if (self.ckbfilename is not None) else None
        self.outfilename = str(self.outfilename) if (self.outfilename is not None) else None
        self.debfilename = str(self.debfilename) if (self.debfilename is not None) else None
    def checkset(self):
        allset = True
        if (self.N is None):
            print(f"N           is not set")
            allset = False
        if (self.L is None):
            print(f"L           is not set")
            allset = False
        if (self.nstab is None):
            print(f"nstab       is not set")
            allset = False
        if (self.north is None):
            print(f"north       is not set")
            allset = False
        if (self.nbin is None):
            print(f"nbin        is not set")
            allset = False
        if (self.nmeassweep is None):
            print(f"nmeassweep  is not set")
            allset = False
        if (self.nskip is None):
            print(f"nskip       is not set")
            allset = False
        if (self.nequil is None):
            print(f"nequil      is not set")
            allset = False
        if (self.dtau is None):
            print(f"dtau        is not set")
            allset = False
        if (self.U is None):
            print(f"U           is not set")
            allset = False
        if (self.mu is None):
            print(f"mu          is not set")
            allset = False
        if (self.ckbfilename is None):
            print(f"ckbfilename is not set")
            allset = False
        if (self.outfilename is None):
            print(f"outfilename is not set")
            allset = False
        if (self.debfilename is None):
            print(f"debfilename is not set")
            allset = False
        return allset