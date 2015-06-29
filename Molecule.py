import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
import pickle
import tkFileDialog

from Green import *
from Graphics import *
from AtomParam import e0, tm

class Molecule (Green):
    """ Manages the description of the molecule
        Attributes
        Atom                    List with the atoms of the molecule
        lpL, lpR                List with the atoms that are coupled to the leads
        S                       Overlap matrix
        V                       Long range interaction matrix for PPP model
        
        Methods
        init
        Load(filename)          Loads a molecule from a file
        Save()                  Saves the current molecule to a file
        AddAtom(Atom, x, y, z)  Adds a certain atom at a specific position to the molecule
        UpdateAtom (iAtom, Atom)Changes the kind of atom at a specific position
        SetLead(iAtom, side)    Couples a specific atom to one of the leads
        SetLeads(lim)           Couples the atoms within a certain limit to the leads
        CreateHam()             Creates the tight binding Hamiltonian
        CreateOverlap()         Creates the overlap matrix
        CreateV()               Creates the long range interaction matrix for the PPP model
        GetCoord()              Returns vectors with the atom positions
        """
    def __init__ (self):
        self.Atom = []
        self.lpL = []
        self.lpR = []
        self.InitGreen()

    def __str__ (self):
        i = 0
        for item in self.Atom:
            print "#", i, ":", item[0], "on", item[1], ",", item [2], ",", item [3]
        return "\n"

    def Load (self, filename):
        f = open(filename, 'r')
        self.Atom = pickle.load(f)
        f.close()
        print "Loaded", filename

    def Save (self):
        filename = tkFileDialog.asksaveasfilename()
        f = open(filename, 'w')
        pickle.dump(self.Atom, f)
        f.close()
        print "Saved", filename

    def AddAtom (self, Atom, x, y, z):
        self.Atom.append((Atom, x, y, z))

    def UpdateAtom (self, iAtom, Atom):
        if Atom=='DEL':
            del self.Atom[iAtom]
        else:
            self.Atom[iAtom] = Atom

    def SetLead (self, iAtom, side):
        if side=='Left':
            self.lpL.append(iAtom)
        elif side=='Right':
            self.lpR.append(iAtom)
        else:
            return
        
    def SetLeads (self, lim=(0,0)):
        Atom, x, y, z = self.GetCoord()
        xmin = min(x)
        xmax = max(x)

        for i in range(len(x)):
            if x[i]-xmin<=lim[0]:
                self.SetLead(i, 'Left')
            if xmax-x[i]<=lim[1]:
                self.SetLead(i, 'Right')

    def CreateHam (self):
        self.N = len(self.Atom)
        L = 1.42 #Angstrom
  
        self.Ham = np.zeros((self.N, self.N), np.dtype('c16'))
        for i, (Atomi, xi, yi, zi) in enumerate(self.Atom):
            for j, (Atomj, xj, yj, zj) in enumerate(self.Atom):
                if i==j:
                    self.Ham[i,j] = e0[Atomi]
                elif sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))<=1.1*L:
                    if Atomi+Atomj in tm:
                        self.Ham[i,j] = tm[Atomi+Atomj]
        return self.Ham

    def CreateOverlap (self):
        L = 1.42 #Angstrom
        
        Stot = np.zeros((self.N, self.N), np.dtype('c16'))
        for i, (Atomi, xi, yi, zi) in enumerate(self.Atom):
            for j, (Atomj, xj, yj, zj) in enumerate(self.Atom):
                if i!=j and sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))<=1.1*L:
                    Stot[i,j] = 0.13
                elif i==j:
                    Stot[i,j] = 1
        EVal, EVec = la.eig(Stot)
        V = np.zeros((len(EVal),len(EVal)), dtype = np.dtype('c16'))
        K = np.zeros((len(EVal),len(EVal)), dtype = np.dtype('c16'))
        for i in range(len(EVal)):
            V[:,i] = EVec[i]
            K[i,i] = 1/sqrt(EVal[i])
        self.S = np.dot(K, np.transpose(np.conjugate(V)))
        self.S = np.dot(V, self.S)
       # print "S^(-1/2):"
       # print self.S

    def CreateV (self):
        self.V = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        for i in range(self.N):
            for j in range(self.N):
                if i!=j:
                    R2 =  (self.Atom[j][1]-self.Atom[i][1])**2+(self.Atom[j][2]-self.Atom[i][2])**2
                    self.V[i,j] += Hubbard/(2*sqrt(1+0.6117*R2))
                else:
                    self.V[i,j] = Hubbard
        return self.V

    def GetCoord (self):
        N = len(self.Atom)
        Atom = []
        x = np.zeros((N,1))
        y = np.zeros((N,1))
        z = np.zeros((N,1))
        for i, (Atomi, xi, yi, zi) in enumerate(self.Atom):
            Atom.append(Atomi)
            x[i] = xi
            y[i] = yi
            z[i] = zi
        return (Atom, x, y, z)
     
class Graphene (Molecule):
    """ Creates a sheet of graphene with width N and length M
        gtype can be 'Zigzag' or 'Armchair'

        Attributes
        
        Methods
        init                    Creates a GNR of type Zigzag or Armchair with width N and lengt M
        nDope()                 Dopes a GNR with Nitrogen
        """
    
    def __init__(self, N, M, gtype):
        self.Atom = []
        self.lpL = []
        self.lpR = []
        self.InitGreen()
        
        L = 1.42 #Angstrom
        
        if gtype == 'Zigzag':
            print "Zigzag"
    
            y = 0
            for i in range (N*2):
                if round(float(i)/2)%2==0:
                    bEven = True
                else:
                    bEven = False

                for j in range (M):
                    if j%2==bEven:
                        self.AddAtom ('C', (j)*L*sqrt(0.75), y, 0)
                    
                if i%2==0:
                    y+= 0.5*L
                else:
                    y+= 1.0*L
                    
        elif gtype == 'Armchair':
            print "Armchair"
            y = 0
            for i in range (N):
                if i%2==0:
                    bEven = True
                else:
                    bEven = False
                for j in range (M+2+int(math.trunc(float(M-1)/2))):
                    if (j+2-bEven)%3!=0:
                        self.AddAtom ('C', (2*j+bEven)*L/2, y, 0)
                    
                y+=L*sqrt(0.75)
        else:
            print "Type of graphene not recognized"

        self.SetLeads()


            
    def nDope (self):
        L = 1.42 #Angstrom
        n=0
        for i, (Atom, x, y, z) in enumerate(self.Atom):
            if y == L*sqrt(0.75):
                n+=1
                if n%2==0:
                    self.UpdateAtom(i, ('N', x, y, z))
        print n/2, "Atoms doped"



        
if __name__ == '__main__':
    from AtomParam import cAtom
  
    Mol1=Molecule()
    Mol1.Load('Molecules/L_Bend_Graphene')
    xmax = 7.37
    ymax = 5*1.42
    print len(Mol1.Atom)
    Num = []
    for i in range(len(Mol1.Atom)-1):
        if Mol1.Atom[i][1]>xmax and Mol1.Atom[i][2]>ymax:
            Num.append(i)
    for i in range (len(Num)):
        Mol1.UpdateAtom(Num[len(Num)-i-1], 'DEL')
    Mol1.UpdateAtom(140, 'DEL')
##    Mol1.Load('Molecules/L_Bend_Graphene')
##    Mol1.AddAtom('C',    0, 0, 0)
##    Mol1.AddAtom('C', 1.42, 0, 0)
   # Mol1.SetLeads((1,1))
    Graph = Graphics()
    Graph.SetMolecule(Mol1)

    plt.show()
