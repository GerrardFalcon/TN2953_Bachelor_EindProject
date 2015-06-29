import matplotlib.pyplot as plt

def InitiateMolecule (Mol, bOverlap=False):
    #Parameters
    bias = 0.1
    Emin, Emax = (-8.2, 8.2)

    
    Mol.CreateHam()
    Mol.SetOverlap(bOverlap)
    Mol.CreateOverlap()
    Mol.CreateV()
    Mol.UpdateGates(Emin, Emax)
    Mol.CreateHamExt (bias, 'Not self consistent')
    Mol.CalcGamma()
    Mol.UpdateGates(Emin, Emax)
    Mol.Density(bias)
    Mol.Transmission()
    Mol.Current(bias)
    Mol.CalcDOS()

def Graph (Mol1, Mol2):
    fig, (ax1) = plt.subplots(1,1)

    ax1.semilogy(Mol1.Gates, Mol1.T, 'k')
    ax1.semilogy(Mol2.Gates, Mol2.T, 'darkmagenta')

##    ax1.semilogy(Mol1.Gates, Mol1.DOS, 'darkorange', linestyle='--')
##    ax1.semilogy(Mol2.Gates, Mol2.DOS, 'darkorange')

    ax1.set_title('Seperate molecules')
#    ax1.legend(('Transmission 1', 'Transmission 2', 'DOS 1', 'DOS 2'))
    ax1.set_xlim([-8.2, 8.2])

##    ax2.plot(Mol1.Gates, Mol2.T/Mol1.T, 'g')
####    ax2.plot(Mol1.Gates, Mol2.DOS-Mol1.DOS, 'darkorange')
##    
##    ax2.set_title('Difference Molecules (2/1)')
##    ax2.set_xlabel('E-Ef (eV)')
##    ax2.legend(('Transmission', 'DOS'))
##    ax2.set_ylim([-1, 10])
    
    
if __name__ == '__main__':
    from Molecule import *
    from AtomParam import cAtom
    from Green import *

    Mol1 = Graphene(7,21, "Armchair")
    InitiateMolecule(Mol1)
    
    Mol2 = Graphene(7,21, "Armchair")
    atom = Mol2.Atom[23]
    atom = ('B', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (23, atom)
    
    atom = Mol2.Atom[27]
    atom = ('B', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (27, atom)
    
    atom = Mol2.Atom[31]
    atom = ('B', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (31, atom)
    
    atom = Mol2.Atom[34]
    atom = ('N', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (34, atom)
    
    atom = Mol2.Atom[38]
    atom = ('N', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (38, atom)
    
    atom = Mol2.Atom[42]
    atom = ('N', atom[1], atom[2], atom[3])
    Mol2.UpdateAtom (42, atom)
    InitiateMolecule(Mol2)
    


    Graph(Mol1, Mol2)

    plt.show()
