import numpy as np
import datetime
from Molecule import *
from AtomParam import cAtom

class TimeTest (object):
    def TestTotal(self, N):
        Natomss   = np.zeros(N)

        timeGreenDiag = np.zeros(N)
        timeGreenInv  = np.zeros(N)
        
        timeTransSeldenthuis = np.zeros(N)
        timeTransMatrix      = np.zeros(N)
        timeTransOrbitals    = np.zeros(N)
        
        timeCurrentMatDiscr10  = np.zeros(N)
        timeCurrentMatDiscr25  = np.zeros(N)
        timeCurrentMatDiscr50  = np.zeros(N)
        timeCurrentMatDiscr100 = np.zeros(N)
        timeCurrentMatAnal     = np.zeros(N)
        
        timeCurrentSeldenthuis = np.zeros(N)
        timeCurrentDiscr10  = np.zeros(N)
        timeCurrentDiscr25  = np.zeros(N)
        timeCurrentDiscr50  = np.zeros(N)
        timeCurrentDiscr100 = np.zeros(N)
        timeCurrentAnal  = np.zeros(N)

        for i in range(N):
            print "Start sequence", i
            Mol = Graphene(3,2+i, "Armchair")
            Mol.CreateHam()
            Mol.CreateHamExt(0.001, 'Not self consistent')
            Mol.CalcGamma()
            Natomss[i] = Mol.N

            #Green
            print "Start Green"
            ti = datetime.datetime.now()
            Mol.GreenDiag(0.1)
            te = datetime.datetime.now()
            delta = te-ti
            timeGreenDiag[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.GreenInv(0.1)
            te = datetime.datetime.now()
            delta = te-ti
            timeGreenInv[i] = delta.total_seconds()*1000
           
            #Transmission
            print "Start Transmission"
            ti = datetime.datetime.now()
            Mol.TransmissionFromSeldenthuis()
            te = datetime.datetime.now()
            delta = te-ti
            timeTransSeldenthuis[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.TransmissionFromMatrix()
            te = datetime.datetime.now()
            delta = te-ti
            timeTransMatrix[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.TransmissionFromOrbitals()
            te = datetime.datetime.now()
            delta = te-ti
            timeTransOrbitals[i] = delta.total_seconds()*1000         
      
            #CurrentMatrix
            print "Start Current Matrix"
            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(1, 0.01, 10)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr10[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(1, 0.01, 25)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr25[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(1, 0.01, 50)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr50[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(1, 0.01, 100)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr100[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentMatrixAnalytic(1, 0.01)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatAnal[i] = delta.total_seconds()*1000

            #Current
##            print "Start Current"
##            ti = datetime.datetime.now()
##            Mol.CurrentFromSeldenthuis(0.01)
##            te = datetime.datetime.now()
##            delta = te-ti
##            timeCurrentSeldenthuis[i] = delta.total_seconds()*1000
##            
##            ti = datetime.datetime.now()
##            Mol.CurrentFromMatrixDiscrete(0.01, 10)
##            te = datetime.datetime.now()
##            delta = te-ti
##            timeCurrentDiscr10[i] = delta.total_seconds()*1000
##
##            ti = datetime.datetime.now()
##            Mol.CurrentFromMatrixDiscrete(0.01, 25)
##            te = datetime.datetime.now()
##            delta = te-ti
##            timeCurrentDiscr25[i] = delta.total_seconds()*1000
##
##            ti = datetime.datetime.now()
##            Mol.CurrentFromMatrixDiscrete(0.01, 50)
##            te = datetime.datetime.now()
##            delta = te-ti
##            timeCurrentDiscr50[i] = delta.total_seconds()*1000
##
##            ti = datetime.datetime.now()
##            Mol.CurrentFromMatrixDiscrete(0.01, 100)
##            te = datetime.datetime.now()
##            delta = te-ti
##            timeCurrentDiscr100[i] = delta.total_seconds()*1000
##
##            if i<20:
##                ti = datetime.datetime.now()
##                Mol.CurrentFromMatrixAnalytic(0.01)
##                te = datetime.datetime.now()
##                delta = te-ti
##                timeCurrentAnal[i] = delta.total_seconds()*1000

            del Mol
            print "Ended sequence"

        fig, ((axGreen, axTrans), (axCurrentMat, axCurrent)) = plt.subplots(2,2)

        axGreen.plot(Natomss, timeGreenDiag, 'r')
        axGreen.plot(Natomss, timeGreenInv, 'b')
        axGreen.legend(('Diagonalized', 'Inverted'))
        axGreen.set_xlabel('Number of atoms')
        axGreen.set_ylabel('Execution time (s)')
        axGreen.set_title('Green')

        axTrans.plot(Natomss, timeTransSeldenthuis, 'r')
        axTrans.plot(Natomss, timeTransMatrix, 'b')
        axTrans.plot(Natomss, timeTransOrbitals, 'g')
        axTrans.legend(('Seldenthuis', 'From Matrix', 'From Orbitals'))
        axTrans.set_xlabel('Number of atoms')
        axTrans.set_ylabel('Execution time (s)')
        axTrans.set_title('Transmission')
            
        axCurrentMat.plot(Natomss, timeCurrentMatDiscr10, 'maroon')
        axCurrentMat.plot(Natomss, timeCurrentMatDiscr25, 'darkred')
        axCurrentMat.plot(Natomss, timeCurrentMatDiscr50, 'r')
        axCurrentMat.plot(Natomss, timeCurrentMatDiscr100, 'salmon')
        axCurrentMat.plot(Natomss, timeCurrentMatAnal, 'b')
        axCurrentMat.legend(('Discrete 10', 'Discrete 25', 'Discrete 50', 'Discrete 100', 'Analytic'))
        axCurrentMat.set_xlabel('Number of atoms')
        axCurrentMat.set_ylabel('Execution time (s)')
        axCurrentMat.set_title('Current matrix')
        
        axCurrent.plot(Natomss, timeCurrentSeldenthuis, 'g')
        axCurrent.plot(Natomss, timeCurrentDiscr10, 'maroon')
        axCurrent.plot(Natomss, timeCurrentDiscr25, 'darkred')
        axCurrent.plot(Natomss, timeCurrentDiscr50, 'r')
        axCurrent.plot(Natomss, timeCurrentDiscr100, 'salmon')
        axCurrent.plot(Natomss, timeCurrentAnal, 'b')
        axCurrent.legend(('Seldenthuis', 'Discrete 10', 'Discrete 25', 'Discrete 50', 'Discrete 100', 'Analytic'))
        axCurrent.set_xlabel('Number of atoms')
        axCurrent.set_ylabel('Execution time (s)')
        axCurrent.set_title('Current')
        plt.show()
        
    def TestGreen(self, N):
        Natomss   = np.zeros(N)
        timeGreenDiag = np.zeros(N)
        timeGreenInv  = np.zeros(N)
        for i in range(N):
            Mol = Graphene(7,i, "Armchair")
            Mol.CreateHam()
            Mol.CreateHamExt(0.001, 'Not self consistent')
            Mol.CalcGamma()
            Natomss[i] = Mol.N
            
            ti = datetime.datetime.now()
            for j in range(10):
                Mol.GreenDiag(j-5)
            te = datetime.datetime.now()
            delta = te-ti
            timeGreenDiag[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            for j in range(10):
                Mol.GreenInv(j-5)
            te = datetime.datetime.now()
            delta = te-ti
            timeGreenInv[i] = delta.total_seconds()*1000
           
            print i
            del Mol
        print "Done"
        fig, axGreen = plt.subplots(1,1)
        axGreen.plot(Natomss, timeGreenDiag, 'r')
        axGreen.plot(Natomss, timeGreenInv, 'b')
        axGreen.legend(('Diagonalized', 'Inverted'))
        axGreen.set_xlabel('Number of atoms')
        axGreen.set_ylabel('Execution time (s)')
        axGreen.set_title('Green')
        plt.show()
        
    def TestTransmission(self, N):
        Natomss   = np.zeros(N)
        timeTransSeldenthuis = np.zeros(N)
        timeTransMatrix      = np.zeros(N)
        timeTransOrbitals    = np.zeros(N)
        for i in range(N):
            Mol = Graphene(7,i, "Armchair")
            Mol.CreateHam()
            Mol.CreateHamExt(0.001, 'Not self consistent')
            Mol.CalcGamma()
            Natomss[i] = Mol.N
            
            ti = datetime.datetime.now()
            Mol.TransmissionFromSeldenthuis()
            te = datetime.datetime.now()
            delta = te-ti
            timeTransSeldenthuis[i] = delta.total_seconds()

            if i<71:
                ti = datetime.datetime.now()
                Mol.TransmissionFromMatrix()
                te = datetime.datetime.now()
                delta = te-ti
                timeTransMatrix[i] = delta.total_seconds()

            ti = datetime.datetime.now()
            Mol.TransmissionFromOrbitals()
            te = datetime.datetime.now()
            delta = te-ti
            timeTransOrbitals[i] = delta.total_seconds()        

            print i
            del Mol
            
        fig, axTrans = plt.subplots(1,1)
        axTrans.plot(Natomss, timeTransSeldenthuis, 'r')
        axTrans.plot(Natomss[0:70], timeTransMatrix[0:70], 'b')
        axTrans.plot(Natomss, timeTransOrbitals, 'g')
        #axTrans.legend(('Seldenthuis', 'From Matrix', 'From Orbitals'))
        axTrans.set_xlabel('Number of atoms')
        axTrans.set_ylabel('Execution time (s)')
        axTrans.set_title('Transmission')
        plt.show()

    def TestCurrentMatrix(self, N=100):
        Natomss   = np.zeros(N)
        timeCurrentMatDiscr10  = np.zeros(N)
        timeCurrentMatDiscr25  = np.zeros(N)
        timeCurrentMatDiscr50  = np.zeros(N)
        timeCurrentMatDiscr100 = np.zeros(N)
        timeCurrentMatAnal     = np.zeros(N)
        timeCurrentMatSolomon  = np.zeros(N)
        for i in range(N):
            Mol = Graphene(7, i, "Armchair")
            Mol.CreateHam()
            Mol.CreateHamExt(0.001, 'Not self consistent')
            Mol.CalcGamma()
            Natomss[i] = Mol.N
            
            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(0.1, 0.01, 10)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr10[i] = delta.total_seconds()

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(0.1, 0.01, 25)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr25[i] = delta.total_seconds()

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(0.1, 0.01, 50)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr50[i] = delta.total_seconds()

            ti = datetime.datetime.now()
            Mol.CurrentMatrixDiscrete(0.1, 0.01, 100)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatDiscr100[i] = delta.total_seconds()
            
            if i<71:
                ti = datetime.datetime.now()
                Mol.CurrentMatrixAnalytic(0.1, 0.01)
                te = datetime.datetime.now()
                delta = te-ti
                timeCurrentMatAnal[i] = delta.total_seconds()

            ti = datetime.datetime.now()
            Mol.LocalCurrents(0.1, 0.01)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentMatSolomon[i] = delta.total_seconds()

            print i
            del Mol
            
        fig, axCurrentMat = plt.subplots(1,1)
        axCurrentMat.semilogy(Natomss, timeCurrentMatDiscr10, 'maroon')
        axCurrentMat.semilogy(Natomss, timeCurrentMatDiscr25, 'darkred')
        axCurrentMat.semilogy(Natomss, timeCurrentMatDiscr50, 'r')
        axCurrentMat.semilogy(Natomss, timeCurrentMatDiscr100, 'salmon')
        axCurrentMat.semilogy(Natomss[0:70], timeCurrentMatAnal[0:70], 'b')
        axCurrentMat.semilogy(Natomss, timeCurrentMatSolomon, 'gold')
       # axCurrentMat.legend(('Discrete 10', 'Discrete 25', 'Discrete 50', 'Discrete 100', 'Analytic', 'Solomon'))
        axCurrentMat.set_xlabel('Number of atoms')
        axCurrentMat.set_ylabel('Execution time (s)')
        axCurrentMat.set_title('Current matrix')
        plt.show()
        
    def TestCurrent(self, N):
        Natoms   = np.zeros(N)
        timeCurrentSeldenthuis = np.zeros(N)
        timeCurrentDiscr10  = np.zeros(N)
        timeCurrentDiscr25  = np.zeros(N)
        timeCurrentDiscr50  = np.zeros(N)
        timeCurrentDiscr100 = np.zeros(N)
        timeCurrentAnal  = np.zeros(N)
        
        for i in range(N):
            Mol = Graphene(7, i, "Armchair")
            Mol.CreateHam()
            Mol.CreateHamExt(0.001, 'Not self consistent')
            Mol.CalcGamma()
            Natoms[i] = Mol.N

            ti = datetime.datetime.now()
            Mol.CurrentFromSeldenthuis(0.01)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentSeldenthuis[i] = delta.total_seconds()*1000
            
            ti = datetime.datetime.now()
            Mol.CurrentFromMatrixDiscrete(0.01, 10)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentDiscr10[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentFromMatrixDiscrete(0.01, 25)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentDiscr25[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentFromMatrixDiscrete(0.01, 50)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentDiscr50[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentFromMatrixDiscrete(0.01, 100)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentDiscr100[i] = delta.total_seconds()*1000

            ti = datetime.datetime.now()
            Mol.CurrentFromMatrixAnalytic(0.01)
            te = datetime.datetime.now()
            delta = te-ti
            timeCurrentAnal[i] = delta.total_seconds()*1000
       

            print i
            del Mol
        print "Done papa"
        fig, axCurrent = plt.subplots(1,1)
        axCurrent.semilogy(Natoms, timeCurrentSeldenthuis, 'g')
        axCurrent.semilogy(Natoms, timeCurrentDiscr10, 'maroon')
        axCurrent.semilogy(Natoms, timeCurrentDiscr25, 'darkred')
        axCurrent.semilogy(Natoms, timeCurrentDiscr50, 'r')
        axCurrent.semilogy(Natoms, timeCurrentDiscr100, 'salmon')
        axCurrent.semilogy(Natoms, timeCurrentAnal, 'b')
        axCurrent.legend(('Seldenthuis', 'Discrete 10', 'Discrete 25', 'Discrete 50', 'Discrete 100', 'Analytic'))
        axCurrent.set_xlabel('Number of atoms')
        axCurrent.set_ylabel('Execution time (ms)')
        axCurrent.set_title('Current')
        plt.show()
        
    def TestSC(self, N):
        Natoms   = np.zeros(N)
        Bias = 0.1
        Delta = [10**-4, 10**-7, 10**-10]
        Beta  = [0.2, 0.5, 1]
        DeltaHubbard = np.zeros((len(Delta),N))
        DeltaPPP     = np.zeros((len(Delta),N))
        BetaHubbard = np.zeros((len(Delta),N))
        BetaPPP     = np.zeros((len(Delta),N))

        for i in range(N):
            Mol = Graphene(7,i, "Armchair")
            Mol.CreateHam()
            Mol.CreateV()
            Mol.CreateHamExt(Bias, 'Not self consistent')
            Mol.CalcGamma()
            Natoms[i] = Mol.N
            
            for j in range(len(Delta)):
                DeltaHubbard[j,i] = Mol.CreateHamExt(Bias, 'Hubbard', Delta[j])
                DeltaPPP[j,i] = Mol.CreateHamExt(Bias, 'PPP', Delta[j])
                BetaHubbard[j,i] = Mol.CreateHamExt(Bias, 'Hubbard', Delta[2], Beta[j])
                BetaPPP[j,i] = Mol.CreateHamExt(Bias, 'PPP', Delta[2], Beta[j])
                

            print i
            del Mol

        fig, ((axDeltaHubbard, axDeltaPPP),(axBetaHubbard, axBetaPPP)) = plt.subplots(2,2)
        color = ['lawngreen', 'turquoise', 'orchid']
        
        color2 = ['moccasin', 'lightsteelblue', 'orchid']
        for j in range(len(Delta)):
            axDeltaHubbard.plot(Natoms, DeltaHubbard[j,:], color[j], linewidth = len(Delta)-j, alpha = 1-0.3*j)
            axDeltaPPP.plot(Natoms, DeltaPPP[j,:],     color[j], linewidth = len(Delta)-j)
            axBetaHubbard.plot(Natoms, BetaHubbard[j,:],     color2[j], linewidth = len(Delta)-j)
            axBetaPPP.plot(Natoms, BetaPPP[j,:],     color2[j], linewidth = len(Delta)-j)
        axDeltaHubbard.set_title('Hubbard different delta')
        axDeltaHubbard.set_xlabel('Number of atoms')
        axDeltaHubbard.set_ylabel('Number of iterations')
        axDeltaPPP.set_title('PPP different delta')
        axDeltaPPP.set_xlabel('Number of atoms')
        axDeltaPPP.set_ylabel('Number of iterations')
        axBetaHubbard.set_title('Hubbard different beta')
        axBetaHubbard.set_xlabel('Number of atoms')
        axBetaHubbard.set_ylabel('Number of iterations')
        axBetaPPP.set_title('PPP different beta')
        axBetaPPP.set_xlabel('Number of atoms')
        axBetaPPP.set_ylabel('Number of iterations')
        plt.show()


if __name__ == '__main__':
    Time = TimeTest()
    Time.TestCurrentMatrix(145)
