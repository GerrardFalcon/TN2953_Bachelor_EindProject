import matplotlib.pyplot as plt
import numpy as np



    
    
if __name__ == '__main__':
    from Molecule import *
    from AtomParam import cAtom
    from Green import *

    Emin, Emax = (-8.12, 8.12)
    Mol = Graphene(7,11, "Armchair")
    Mol.CreateHam()

    gL = np.linspace(0, 1.0, 100)
    TgLreg1 = np.zeros(len(gL))    
    TgLreg2 = np.zeros(len(gL))
    TgLreg3 = np.zeros(len(gL))

    reg1 = [426, 573]
    reg2 = [0, 228, 771, 999]
    reg3 = [229, 425, 574, 770 ]

    for i in range(len(gL)):
        Mol.SetG(gL[i], 0.2)
        Mol.CreateHamExt (0.1, 'Not self consistent')
        Mol.CalcGamma()
        Mol.UpdateGates(Emin, Emax)
        Mol.Transmission()
        TgLreg1[i] = np.mean(Mol.T[reg1[0]:reg1[1]])
        TgLreg2[i] = (np.mean(Mol.T[reg2[0]:reg2[1]])+np.mean(Mol.T[reg2[2]:reg2[3]]))/2
        TgLreg3[i] = (np.mean(Mol.T[reg3[0]:reg3[1]])+np.mean(Mol.T[reg3[2]:reg3[3]]))/2
        print i
        
    gR = gL
    TgRreg1 = np.zeros(len(gL))    
    TgRreg2 = np.zeros(len(gL))
    TgRreg3 = np.zeros(len(gL))

    for i in range(len(gL)):
        Mol.SetG(0.2, gR[i])
        Mol.CreateHamExt (0.1, 'Not self consistent')
        Mol.CalcGamma()
        Mol.UpdateGates(Emin, Emax)
        Mol.Transmission()
        TgRreg1[i] = np.mean(Mol.T[reg1[0]:reg1[1]])
        TgRreg2[i] = (np.mean(Mol.T[reg2[0]:reg2[1]])+np.mean(Mol.T[reg2[2]:reg2[3]]))/2
        TgRreg3[i] = (np.mean(Mol.T[reg3[0]:reg3[1]])+np.mean(Mol.T[reg3[2]:reg3[3]]))/2
        print i

    TgLRreg1 = np.zeros(len(gL))    
    TgLRreg2 = np.zeros(len(gL))
    TgLRreg3 = np.zeros(len(gL))

    for i in range(len(gL)):
        Mol.SetG(gL[i], gR[i])
        Mol.CreateHamExt (0.1, 'Not self consistent')
        Mol.CalcGamma()
        Mol.UpdateGates(Emin, Emax)
        Mol.Transmission()
        TgLRreg1[i] = np.mean(Mol.T[reg1[0]:reg1[1]])
        TgLRreg2[i] = (np.mean(Mol.T[reg2[0]:reg2[1]])+np.mean(Mol.T[reg2[2]:reg2[3]]))/2
        TgLRreg3[i] = (np.mean(Mol.T[reg3[0]:reg3[1]])+np.mean(Mol.T[reg3[2]:reg3[3]]))/2
        print i

    fig, (ax1, ax2, ax3) = plt.subplots(3,1)

    ax1.semilogy(gL, TgLreg1, 'g')
    ax1.semilogy(gL, TgLreg2, 'r')
    ax1.semilogy(gL, TgLreg3, 'y', linestyle='_')
    ax1.set_xlabel('gL')
    ax1.set_ylabel('Tranmission')
    
    ax2.semilogy(gL, TgRreg1, 'g')
    ax2.semilogy(gL, TgRreg2, 'r')
    ax2.semilogy(gL, TgRreg3, 'y', linestyle='_')
    ax1.set_xlabel('gR')
    ax1.set_ylabel('Tranmission')
    
    ax3.semilogy(gL, TgLRreg1, 'g')
    ax3.semilogy(gL, TgLRreg2, 'r')
    ax3.semilogy(gL, TgLRreg3, 'y', linestyle='_')
    ax1.set_xlabel('gL, gR')
    ax1.set_ylabel('Tranmission')

    

    plt.show()
