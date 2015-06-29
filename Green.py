import numpy as np
from numpy import sqrt
import numpy.linalg as la
import matplotlib.pyplot as plt
from copy import copy, deepcopy
import math as math
from AtomParam import gLeft, gRight, Hubbard
import timeit

def DisplayMatrix (Matrix):
    fig, ax = plt.subplots(figsize=(20,11))
    table = plt.table(cellText = np.round(Matrix, 8), cellLoc='center', bbox=[0, 0, 1, 1])
    table.set_fontsize(25)
    ax.axis('off')
    plt.show()

class Green(object):
    """ Represents a Molecule including hamiltonian and transmission and current calculations 
        Attributes
        HamExt                      Hamiltonian of the Molecule, with interaction with leads
        N                           Size of the Hamiltonian
        Gleft, Gright               Interaction strengths with the leads
        gam_l, gam_r                Gamma matrices
        eigval, eigvec              Respectively real and imaginary part of the eigenvalues of the Hamiltonian
        C                           Matrix with the eigenvectors
        Gates                       Array of gate voltages which have to be calculated
        T,P,I                       Transmission amplitude and phase and current

        Methods
        InitGreen()                 Sets the initial parameters for the interactino with the leads
        SetG(Gleft, Gright)         Set the interaction strengths with the leads
        SetOverlap()
        CreateHamExt(Bias, scType, delta, Beta)             Creates the Hamiltonian that includes the interaction with the leads in different consistency ways (no consistency, Hubbard, Coulomb)
        GetHamExt(D, scType)
        CalcEigValVec()     
        UpdateGates(Gmin,Gmax)      Sets the range of the gate voltages betweeen the lowest and heighest eigenenergie plus 5 times gamma
        CalcGamma()                 Calculates the gamma matrices in the basis of the HamExt
        CalcGe()
        CalcBetaMat()               Calculates the Beta matrices
        
        Green(E)                    Calculates the Green matrix
        GreenDiag(E)
        GreenInv (E)

        Density(Bias)               Calculates the Density matrix, with an upper and lower energy boundary
        DensityDiscrete(Bias)
        DensityAnalytic(Bias)

        CalcDOS()

        TransmissionMatrix(E)       Calculates the complete transmission matrix

        Transmission()              Calculates the transmission based on one of the methods below
        TransmissionFromSeldenthuis Calculates the transmission based on the expression from Seldenthuis
        TransmissionFromMatrix()    Calculates the transmission based on the trace of the transmission matrix
        TransmissionFromOrbitals()  Calculates the transmission based on the individual transmissions from the orbitals

        TransmissionOrbital(iOrbit) Calculates the transmission of an orbital based on one of the methods below
        TransmissionOrbitalFromSeldenthuis(iOrbit) Calculates the transmission of an orbital based on the method described by Seldenthuis
        
        Current(bias)               Calculates the current based on one of the methods below
        CurrentFromSeldenthuis(bias)Calculates the current based on the expression from Seldenthuis
        CurrentFromMatrixDiscrete(bias) Calculates the current from the trace of the discretely obtained current matrix
        CurrentFromMatrixAnalytic(bias) Calculates the current from the trace of the analytically obtained current matrix

        CurrentMatrix(bias, E)      Calculates the current matrix based on one of the methods below
        CurrentMatrixDiscrete(bias, E, Nstep) Calculates the current matrix based on a Riemman sum of the transmission matrices
        CurrentMatrixAnalytic(bias, E)        Calculates the current matrix based on an analytic expression

        LocalCurrents(bias, E)      Calculates the local currents based on an expression by Solomon
        """
    
    def InitGreen(self):
        self.Gleft  = gLeft
        self.Gright = gRight
        self.bOverlap = False

    def SetG (self, Gleft, Gright):
        self.Gleft  = Gleft
        self.Gright = Gright

    def SetOverlap(self, bOverlap=False):
        self.bOverlap=bOverlap

    def CreateHamExt (self, Bias, scType, delta=10**-10, Beta=1.0):
        self.HamExt = deepcopy(self.Ham)
        """g models the interaction with the contact leads at the end atoms of the molecule"""
        for i in self.lpL:
            self.HamExt[i,i] -= 0.5j*self.Gleft
        for i in self.lpR:
            self.HamExt[i,i] -= 0.5j*self.Gright

        if scType == 'Not self consistent':
            return self.HamExt
    
        HamExtOrg = deepcopy(self.HamExt)
        Nmax    = 100 #maximum number of iterations
        RhoAvg  = np.zeros((Nmax), dtype=np.dtype('c16'))
        
        Dp = self.Density(Bias)
        RhoAvg[0] = np.trace(Dp)/self.N
        print "Average density", RhoAvg[0], "at initial iteration"
        self.HamExt = HamExtOrg + self.GetHamExt(Dp, scType)
        self.CalcEigValVec()
        Dn = self.Density(Bias)        
        RhoAvg[1] = np.trace(Dn)/self.N
        print "Average density", RhoAvg[1], "at second iteration"
        if abs(RhoAvg[1]-RhoAvg[0])<=delta:
            print "Final delta", abs(RhoAvg[1]-RhoAvg[0]), "<", delta
            print 1, "iteration"
            return 1

        for i in range(Nmax-2):

            self.HamExt = HamExtOrg + self.GetHamExt(Beta*Dn+(1-Beta)*Dp, scType)
            self.CalcEigValVec()
            Dp = deepcopy(Dn)
            Dn = self.Density(Bias)
            RhoAvg[i+2]=np.trace(Dn)/self.N
            #print i
            print "Average density", RhoAvg[i+2], "at iteration", i+2
            nIter = i+2
            
            if abs(RhoAvg[i+2]-Beta*RhoAvg[i+1]-(1-Beta)*RhoAvg[i])<=delta:
                nIter = i+2
                print "Final delta", abs(RhoAvg[i+2]-Beta*RhoAvg[i+1]-(1-Beta)*RhoAvg[i]), "<", delta
                print i+2, "iterations"
                break
            
            
            
            
        fig, ax = plt.subplots(1,1)
        ax.plot(RhoAvg[0:nIter+1])
        ax.set_title(scType + ", Delta=" + str(delta) + ", N=" + str(self.N) + ", Beta=" + str(Beta))
        return nIter

    def GetHamExt (self, D, scType):
        HamExtAdd = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        if scType=='Hubbard':
            for i in range(self.N):
                HamExtAdd[i,i] += (D[i,i]*0.5)*self.V[i,i]
                
        if scType=='PPP':
            for i in range (self.N):
                for j in range(self.N):
                    HamExtAdd[i,i] += (D[j,j])*self.V[i,j]
        return HamExtAdd
                   
    def CalcEigValVec(self):
        """ Calculate eigenvalues, normalized eigenvectors in order to diagonalize Green matrix """
        if self.bOverlap:
            SHS = np.dot(self.HamExt, self.S)
            SHs = np.dot(self.S, SHS)
            self.eigval,self.eigvec = la.eig(SHS)                                             
            for i in range(self.N):
                self.eigvec[i] = np.dot(self.S, self.eigvec[i])
        else:
            self.eigval,self.eigvec = la.eig(self.HamExt)

        self.eigvec /= np.sqrt((self.eigvec ** 2).sum(0))
        la.norm(self.eigvec[:,0]), la.norm(self.eigvec[:,1])
        dummy = np.real(self.eigval)
        idx = dummy.argsort()
        self.eigval = self.eigval[idx]
        self.eigvec = self.eigvec[idx]

##        print "Eigenvalues:"
##        print self.eigval
##        print "Eigenvectors:"
##        print self.eigvec
        
        self.C = np.zeros((self.N, self.N), dtype = np.dtype('c16'))
        for i in range(self.N):
            for j in range(self.N):
                self.C[j,i] = self.eigvec[i][j]

##        print "C:"
##        print self.C
##
##        print "Gdiag:"
##        print self.GreenDiag(0).round(3)
##        print "Ginv:"
##        print self.GreenInv(0)
            
        self.e_arr = np.real(self.eigval)
        self.gam   = -np.imag(self.eigval)
        self.UpdateGates()
    
    def UpdateGates (self, Gmin=None, Gmax=None):
        if Gmin==None:
            Gmin = min(self.e_arr)-5*max(abs(self.gam))
        if Gmax==None:
            Gmax = max(self.e_arr)+5*max(abs(self.gam))
        self.Gates = np.arange(Gmin, Gmax, float(Gmax-Gmin)/1000)
        return self.Gates
        
    def CalcGe (self, e):
        ge = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        for i in range(self.N):
            ge[i,i] = 1/(e-self.eigval[i])
        return ge

    def CalcGamma(self):
        """Gamma Matrices, modeling of the interaction with both leads"""
        gmat_l = np.zeros((self.N, self.N),np.dtype('c16')) 
        gmat_r = np.zeros((self.N, self.N),np.dtype('c16'))
        for i in self.lpL:
            gmat_l[i,i] = 0.5j*self.Gleft
        for i in self.lpR:
            gmat_r[i,i] = 0.5j*self.Gright

        self.CalcEigValVec()

        """ Diagonalize Gamma matrices into basis of HamExt"""
        self.gam_l = np.dot(gmat_l, np.conj(self.eigvec))
        self.gam_l = np.dot(np.transpose(self.eigvec), self.gam_l)
        self.gam_r = np.dot(gmat_r, np.conj(self.eigvec))
        self.gam_r = np.dot(np.transpose(self.eigvec), self.gam_r)
##        print "Gam L:"
##        print self.gam_l
##        print "Gam R:"
##        print self.gam_r

        self.CalcBetaMat ()
          
    def CalcBetaMat (self):
        """ Calculates the beta's used in the calculation of the current. """
        #self.beta_mat = np.zeros((self.N,self.N),dtype=float)
        self.Beta = np.zeros((self.N,self.N),dtype=np.dtype('c16'))
        for i in range(self.N):
           for j in range(self.N):
              #self.beta_mat[i,j] = (self.e_arr[i]-self.e_arr[j])/(self.gam[i]+self.gam[j])
               self.Beta[i,j] = 1j/(self.gam[i]+self.gam[j]+1j*(self.e_arr[i]-self.e_arr[j]))
        return self.Beta
    
        """-----------------------------------Green----------------------------------------"""
    def Green (self, e):
        return self.GreenInv (e)   #Atom basis

    def GreenDiag(self, e):
        G = np.dot(self.CalcGe(e), np.transpose(self.C))
        G = np.dot(self.C, G)
        return G    #Atom basis

    def GreenInv (self, e):
        G = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        np.fill_diagonal(G, e)
        return la.inv(G-self.HamExt)    #Atom basis
    
        """----------------------------------Density---------------------------------------"""
    def Density(self, Bias):
        return self.DensityDiscrete(Bias)
    
    def DensityDiscrete(self, Bias, Nsteps=100):
        D = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        Emin = min(self.Gates)
        Emax = Bias/2
        de = float(Emin-Emax)/Nsteps
        for i in range(Nsteps):
            e = Emin+(i+0.5)*de
            D += -1/math.pi*np.imag(self.Green(e))*de
        self.LD = np.real(abs(D))
        return self.LD  #Atom basis

##    def DensityAnalytic(self, Bias): # Does not work
##        Emin = min(self.Gates)
##        Emax = Bias/2
##        Geint = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
##        for i in range(self.N):
##            tan_term_plus = np.arctan(Emax**2/self.e_arr[i]**2)
##            tan_term_min  = np.arctan(Emin**2/self.e_arr[i]**2)
##            Geint[i,i] = -self.gam[i]/self.e_arr[i]**2*(tan_term_plus-tan_term_min)
##        LD = -1/math.pi*np.dot(Geint, self.C)
##        LD = np.dot(np.transpose(self.C), LD)
##        return np.real(LD)
    
        """------------------------------------DOS-----------------------------------------"""
    def CalcDOS (self):
        self.DOS = np.zeros((len(self.Gates)),dtype=np.dtype('c16'))

        for i, ei in enumerate(self.Gates):
            self.DOS[i] = -1/math.pi*np.trace(np.imag(self.Green(ei)))
                
        """----------------------------TransmissionMatrix----------------------------------"""
    def TransmissionMatrix(self, E):
        ge = self.CalcGe(E)
            
        TransMat = np.dot(self.gam_r, np.conj(ge))
        TransMat = np.dot(ge,         TransMat)
        TransMat = np.dot(self.gam_l, TransMat)
        return TransMat     #Orbital basis

##    def TransmissionMatrix2(self, e):
##        gmat_l = np.zeros((self.N, self.N),np.dtype('c16')) 
##        gmat_r = np.zeros((self.N, self.N),np.dtype('c16'))
##        for i in self.lpL:
##            gmat_l[i,i] = 0.5j*self.Gleft
##        for i in self.lpR:
##            gmat_r[i,i] = 0.5j*self.Gright
##
##        G = self.Green(e)
##        TransMat = np.dot(np.conjugate(G), self.C)
##        TransMat = np.dot(gmat_r, TransMat)
##        TransMat = np.dot(G, TransMat)
##        TransMat = np.dot(gmat_l, TransMat)
##        TransMat = np.dot(np.transpose(self.C), TransMat)
##        return TransMat
##
##    def TransMat (self, e):
##        gmat_l = np.zeros((self.N, self.N),np.dtype('c16')) 
##        gmat_r = np.zeros((self.N, self.N),np.dtype('c16'))
##        for i in self.lpL:
##            gmat_l[i,i] = 0.5j*self.Gleft
##        for i in self.lpR:
##            gmat_r[i,i] = 0.5j*self.Gright
##        G = self.Green(e)
##        print "Gmat_L:"
##        print gmat_l
##        print "Gmat_R:"
##        print gmat_r
##        print "G:"
##        print G
##        T = np.dot(gmat_r, np.conjugate(G))
##        T = np.dot(G, T)
##        T = np.dot(gmat_l, T)
##        return T
        
        """--------------------------------Transmission------------------------------------"""
    def Transmission (self):
        trans = self.TransmissionFromSeldenthuis ()
        self.T = np.real(trans)                     # Only for From Seldenthuis, From Matrix
        #self.T = np.real(trans*np.conjugate(trans)) #Only for From Orbitals
        self.P = np.arctan(np.imag(trans)/np.real(trans))
        return self.T

    def TransmissionFromSeldenthuis (self):
        """Calculate transmission"""
        trans = np.zeros((len(self.Gates)),dtype=np.dtype('c16'))
  
        for i in range(self.N):
            for j in range(self.N):
                trans = trans - self.gam_l[i,j]*self.gam_r[i,j]/((self.Gates - self.e_arr[i] + 1j*self.gam[i])*(self.Gates - self.e_arr[j] - 1j*self.gam[j]))
        return trans

    def TransmissionFromMatrix (self):
        trans = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))

        for i in range(len(self.Gates)):
            trans[i] = np.trace(self.TransmissionMatrix(self.Gates[i]))
        return -trans

    def TransmissionFromOrbitals (self):
        trans = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        d0    = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        d1    = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        d2    = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        for i in range(self.N):
            d0, d1, d2 = self.TransmissionOrbital(i)
            trans+=d0
        return trans
    
        """--------------------------Transmission - Orbital--------------------------------"""
    def TransmissionOrbital (self, iOrbital):
        return self.TransmissionOrbitalFromSeldenthuis (iOrbital)
    
    def TransmissionOrbitalFromSeldenthuis (self, iOrbital):
        
        Ttot = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        Tamp = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))
        Tpha = np.zeros((len(self.Gates)), dtype=np.dtype('c16'))

        gmat_l = np.zeros((self.N, self.N),np.dtype('c16')) 
        gmat_r = np.zeros((self.N, self.N),np.dtype('c16'))
        for i in self.lpL:
            gmat_l[i,i] = 0.5j*self.Gleft
        for i in self.lpR:
            gmat_r[i,i] = 0.5j*self.Gright

        vL=0
        vR=0
        for alpha in range(self.N):
            vL += self.eigvec[alpha][iOrbital]*sqrt(gmat_l[alpha,alpha])
            vR += self.eigvec[alpha][iOrbital]*sqrt(gmat_r[alpha,alpha])
##        print vL
##        print vR
        Ttot = vL*vR/(self.Gates-self.eigval[iOrbital])
        Amplitude = np.real(Ttot*np.conjugate(Ttot))
        Phase     = np.arctan(np.imag(Ttot)/np.real(Ttot))
        return (Ttot, Amplitude, Phase)

##    def TransmissionOrbitalExt(self, iOrbital, bias):
##        """Calculate phase information"""        
##
##        mu_L = self.Gates+0.5*bias
##        mu_R = self.Gates-0.5*bias
##
##        len = mu_L.size
##        gmat_l = np.zeros((self.N, self.N),np.dtype('c16')) 
##        gmat_r = np.zeros((self.N, self.N),np.dtype('c16'))
##        for i in self.lpL:
##            gmat_l[i,i] = 0.5j*self.Gleft
##        for i in self.lpR:
##            gmat_r[i,i] = 0.5j*self.Gright
##            
##
##        
##        t   = np.zeros((2,len), dtype=np.dtype('c16'))
##        T   = np.zeros((2,len), dtype=np.dtype('c16'))
##        phi = np.zeros((2,len), dtype=np.dtype('c16'))
##        for iOrbital in range(2):
##            vL = 0
##            vR = 0
##            print iOrbital
##            print self.eigvec
##            for alpha in range(self.N):
##                print alpha
##                print self.eigvec[alpha][iOrbital]*sqrt(gmat_l[alpha,alpha])
##                vL += self.eigvec[alpha][iOrbital]*sqrt(gmat_l[alpha,alpha])
##                vR += self.eigvec[alpha][iOrbital]*sqrt(gmat_r[alpha,alpha])
##
##            if vL*vR>0:
##                phi0 = -math.pi/2
##            else:
##                phi0 = math.pi/2
##            print "VL = ", vL
##            print "VR = ", vR
##            t0 = abs(vL*vR)/np.imag(self.eigval[iOrbital])
##            print t0
##            for j, ei in enumerate(self.Gates):
##                T[iOrbital, j] = vL*vR/(ei-self.eigval[iOrbital])
##                t  [iOrbital,j] = t0*sqrt(np.imag(self.eigval[iOrbital])**2/((np.real(self.eigval[iOrbital])-ei)**2+np.imag(self.eigval[iOrbital])**2))
##                phi[iOrbital,j] = phi0 + math.atan((self.Gates[j]-np.real(self.eigval[iOrbital]))/np.imag(self.eigval[iOrbital]))
##        fig, (Amplitude, Phase) = plt.subplots(2,1)
##        Amplitude.semilogy(self.Gates,  -t[0,:], 'r')
##        Amplitude.semilogy(self.Gates,  -t[1,:], 'b')
##        Amplitude.set_ylabel('Transmission')
##
##        """Calculate total from individual transmissions"""
####        Ttot   = np.zeros((1,len), dtype=np.dtype('c16'))
####        Ttot = T[0, :]+T[1, :]
####        for i, Ttoti in enumerate(Ttot):
####            Ttot[i] = Ttoti*np.conjugate(Ttoti)
####        Amplitude.semilogy(self.Gates, Ttot, 'y.')
##
##        """Calculate total from transmission matrix"""
##        Ttot2   = np.zeros((len), dtype=np.dtype('c16'))
##        t2   = np.zeros((4,len), dtype=np.dtype('c16'))
##        for j, ei in enumerate(self.Gates):
##            TransMat = self.TransmissionMatrix(ei)
##            Ttot2[j] = TransMat[0,0]+TransMat[1,1]
####            TransMat = np.dot(TransMat, self.C)
####            TransMat = np.dot(np.transpose(self.C), TransMat)
##            t2[0,j] = sqrt(TransMat[0,0]*np.conjugate(TransMat[0,0]))
##            t2[1,j] = sqrt(TransMat[1,1]*np.conjugate(TransMat[1,1]))
##            t2[2,j] = sqrt(TransMat[1,0]*np.conjugate(TransMat[1,0]))
##            t2[3,j] = sqrt(TransMat[0,1]*np.conjugate(TransMat[0,1]))
##        #Amplitude.semilogy(self.Gates, -Ttot2, 'm*')
##
##            
##        Amplitude.semilogy(self.Gates, self.T, 'g')
####
####        print t2[0, :]
####        print t2[1, :]
##        Amplitude.semilogy(self.Gates,  t2[0,:], 'dimgray', linewidth = 5)
##        Amplitude.semilogy(self.Gates,  t2[1,:], 'lightgray', linewidth = 5)
##        Amplitude.semilogy(self.Gates,  t2[2,:], 'yellow')
##        Amplitude.semilogy(self.Gates,  t2[3,:], 'gold')
####        #Amplitude.semilogy(self.Gates,  +t2[1,:], 'b.')
##        Amplitude.legend(('Orbital 1', 'Orbital 2', 'Total', '00', '11', '10', '01'))
##
##        Phase.plot(self.Gates, phi[0,:],'r')
##        Phase.plot(self.Gates, phi[1,:],'b')
##        Phase.set_xlabel('E-Ef (eV)')
##        Phase.set_ylabel('Phase')
##        
##        
##        return [t, phi]

        """----------------------------------Current---------------------------------------"""
    def Current (self, bias):
        self.I = self.CurrentFromSeldenthuis (bias)
        return self.I
    
    def CurrentFromSeldenthuis (self, bias):
        mu_L = self.Gates+0.5*bias
        mu_R = self.Gates-0.5*bias

        len = mu_L.size
        
        """Calculate current"""
        I = np.zeros((len),dtype=np.dtype('c16'))
        log_term_plus = np.zeros((self.N,len),np.dtype('c16'))
        log_term_minus = np.zeros((self.N,len),np.dtype('c16'))
        for i in range(self.N):
            log_term_plus [i,:] = np.log((mu_R-self.e_arr[i] + 1j*self.gam[i])/ \
                                 (mu_L-self.e_arr[i] + 1j*self.gam[i]))
            log_term_minus[i,:] = np.log((mu_R-self.e_arr[i] - 1j*self.gam[i])/ \
                                 (mu_L-self.e_arr[i] - 1j*self.gam[i]))
        for i in range(self.N):
            for j in range(self.N):
                #I = I + self.gam_l[i,j]*self.gam_r[i,j]/(self.gam[i] + self.gam[j])*1j/(1+1j*self.beta_mat[i,j])*(log_term_plus[i,:] - log_term_minus[j,:])
                I += self.gam_l[i,j]*self.gam_r[i,j]*self.Beta[i,j]*(log_term_plus[i,:]-log_term_minus[j,:])
        return np.real(I)

    def CurrentFromMatrixDiscrete (self, bias, Nsteps=100):
        Itot = np.zeros((len(self.Gates)),dtype=np.dtype('c16'))
        for i, ei in enumerate(self.Gates):
            CM = self.CurrentMatrixDiscrete(ei, bias, Nsteps)
            Itot[i] = np.trace(CM*np.conjugate(CM))
        return Itot
        
    def CurrentFromMatrixAnalytic (self, bias):
        Itot = np.zeros((len(self.Gates)),dtype=np.dtype('c16'))
        CM = self.CurrentMatrixAnalytic(bias)
        for i in range(len(self.Gates)):
            Itot[i] = np.trace(CM[:,:,i])
        return Itot

##    def CurrentDiscrete (self, bias, Nstep=100): #Doesn't work yet, complicated to obtain energy dependent transmission
##        I = np.zeros((len),dtype=np.dtype('c16'))
##        mu_L = e+bias/2
##        mu_R = e-bias/2
##        de = float(mu_L-mu_R)/Nstep
##        for i in range(Nsteps):
##            e = mu_R+(i+0.5)*de
##            I += self.Transmission()*de
##        return I

        """-------------------------------CurrentMatrix------------------------------------"""
    def CurrentMatrix (self, bias, e):
        return self.CurrentMatrixAnalytic (bias, e)
    
    def CurrentMatrixDiscrete (self, bias, e, Nstep = 100):
        CurrentMat = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        mu_L = e+bias/2
        mu_R = e-bias/2
        de = float(mu_L-mu_R)/Nstep
        for i in range(Nstep):
            e = mu_R+(i+0.5)*de
            CurrentMat +=self.TransmissionMatrix(e)*de
        CurrentMat = np.dot(CurrentMat, np.transpose(self.C))
        CurrentMat = np.dot(self.C, CurrentMat)
        return CurrentMat
    
    def CurrentMatrixAnalytic (self, bias, e=None):
        if e==None:
            mu_L = self.Gates+0.5*bias
            mu_R = self.Gates-0.5*bias
            Ne = len(mu_L)
        else:
            mu_L = e+0.5*bias
            mu_R = e-0.5*bias
            Ne=1
            
        CurrentMat = np.zeros((self.N, self.N, Ne), dtype=np.dtype('c16'))
        log_term_plus = np.zeros((self.N, Ne),np.dtype('c16'))
        log_term_minus = np.zeros((self.N, Ne),np.dtype('c16'))
        for i in range(self.N):
            log_term_plus [i,:] = np.log((mu_R-self.e_arr[i] + 1j*self.gam[i])/ \
                                 (mu_L-self.e_arr[i] + 1j*self.gam[i]))
            log_term_minus[i,:] = np.log((mu_R-self.e_arr[i] - 1j*self.gam[i])/ \
                                 (mu_L-self.e_arr[i] - 1j*self.gam[i]))
            
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    CurrentMat[i,j,:] += self.gam_l[i,k]*self.gam_r[k,j]*self.Beta[k,j]*(log_term_plus[k]-log_term_minus[j])
        
      
        for i in range(Ne):
            CurrentMat[:,:,i] = np.dot(CurrentMat[:,:,i], self.C)
            CurrentMat[:,:,i] = np.dot(np.transpose(self.C), CurrentMat[:,:,i])
##        print "Trace of current matrix:"
##        print np.trace(CurrentMat)
##        print "Sum of off diagonal elements:"
##        print CurrentMat.sum()-np.trace(CurrentMat)
        return CurrentMat

    def LocalCurrents(self, bias, e):
        Ga = self.Green(e)
        Gr = np.conjugate(Ga)
        K= np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        Gam_L = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        Gam_R = np.zeros((self.N, self.N), dtype=np.dtype('c16'))
        for i in self.lpL:
            Gam_L [i,i] -= 0.5j*self.Gleft
        for i in self.lpR:
            Gam_R [i,i] -= 0.5j*self.Gright
##        print "Gam_L:"
##        print Gam_L
##        print "Gam_R:"
##        print Gam_R
##        print "G:"
##        print G
        GGGL = np.dot(Gam_L, Ga)
        GGGL = np.dot(Gr, GGGL)
        GGGR = np.dot(Gam_R, Ga)
        GGGR = np.dot(Gr, GGGR)

        for i in range(self.N):
            for j in range(self.N):
                K[i,j] = self.Ham[i,j]*GGGL[j,i]-self.Ham[j,i]*GGGL[i,j]

##        for i in range(self.N):
##            for j in range(self.N):
##                for k in range(self.N):
##                    for l in range(self.N):
##                        if i is not j:
##                            K[i,j] += 1j*(self.Ham[i,j]*Gr[j,k]*Gam_L[k,l]*Ga[l,i] - self.Ham[j,i]*Gr[i,l]*Gam_L[l,k]*Ga[k,j])
##        print K
        return K
        

def TestTransmissionMatrixMethods (Mol):
    TM1 = np.zeros((Mol.N, Mol.N, len(Mol.Gates)))
    TM2 = np.zeros((Mol.N, Mol.N, len(Mol.Gates)))

    for i, ei in enumerate(Mol.Gates):
        TM1[:,:,i] = Mol.TransmissionMatrix(ei)
        TM2[:,:,i] = Mol.TransmissionMatrix2(ei)

    axTM = np.ndarray((Mol.N, Mol.N), dtype=np.object)
    fig, axTM = plt.subplots(Mol.N, Mol.N)

    for i in range(Mol.N):
        for j in range(Mol.N):
            axTM[i,j].semilogy(Mol.Gates, TM1[i,j,:]*np.conjugate(TM1[i,j,:]), 'r', linewidth = 3, alpha=0.5)
            axTM[i,j].semilogy(Mol.Gates, TM2[i,j,:]*np.conjugate(TM2[i,j,:]), 'g')
    fig.suptitle("Transmission Matrix Methods")
            
def TestTransmissionMethods (Mol):
    Tseldenthuis = Mol.TransmissionFromSeldenthuis()
    Tindivorbit  = Mol.TransmissionFromOrbitals()
    Tmatrix      = Mol.TransmissionFromMatrix()

    fig, (Amplitude, Phase) = plt.subplots(2,1)
    fig.suptitle("Transmission Methods")

    Amplitude.semilogy(Mol.Gates, np.real(Tseldenthuis), 'r', linewidth = 5, alpha=0.3)
    Amplitude.semilogy(Mol.Gates, Tindivorbit*np.conjugate(Tindivorbit), 'g', linewidth = 3, alpha = 0.6)
    Amplitude.semilogy(Mol.Gates, np.real(Tmatrix), 'b')
    Amplitude.legend(('Seldenthuis', 'From Orbitals', 'From matrix'))
    Phase.plot(Mol.Gates, np.imag(Tseldenthuis), 'r', linewidth = 5, alpha = 0.3)
    Phase.plot(Mol.Gates, np.imag(Tindivorbit), 'g', linewidth = 3, alpha = 0.6)
    Phase.plot(Mol.Gates, np.imag(Tmatrix), 'b')

def TestTransmissionOrbitalMethods (Mol):
    TSeldenthuis = np.zeros((Mol.N, len(Mol.Gates)), dtype=np.dtype('c16'))
    PSeldenthuis = np.zeros((Mol.N, len(Mol.Gates)), dtype=np.dtype('c16'))
    for i in range(Mol.N):
        T, TSeldenthuis[i,:], PSeldenthuis[i,:]= Mol.TransmissionOrbitalFromSeldenthuis(i)

    AFromMatrix = np.zeros((Mol.N, len(Mol.Gates)), dtype=np.dtype('c16'))
    PFromMatrix = np.zeros((Mol.N, len(Mol.Gates)), dtype=np.dtype('c16'))
    for j, ei in enumerate(Mol.Gates):
        T = Mol.TransmissionMatrix(ei)
        for i in range(Mol.N):
            AFromMatrix[i,j] = np.real(T[i,i])
            PFromMatrix[i,j] = np.arctan(np.imag(T[i,i])/np.real(T[i,i]))     #rotated by pi/2???  
        


    fig, (Amplitude, Phase) = plt.subplots(2,1)
    fig.suptitle("Transmission of Orbitals methods")
    for i in range(Mol.N):
        Amplitude.semilogy (Mol.Gates, TSeldenthuis[i,:])

    for i in range(Mol.N):
        Amplitude.semilogy (Mol.Gates, abs(AFromMatrix[i,:]), linewidth=3, alpha=0.5)

    for i in range(Mol.N):
        Phase.plot (Mol.Gates, PSeldenthuis[i,:])
    
    for i in range(Mol.N):
        Phase.plot (Mol.Gates, PFromMatrix[i,:], linewidth=3, alpha=0.5)

    Phase.set_xlabel('E-Ef')
    Amplitude.set_title('Transmission')
    Amplitude.set_ylabel('Amplitude')
    Phase.set_ylabel('Phase')
    fig.suptitle("Orbital Transmission Methods")

def TestCurrentMatrixMethods (Mol):
    CM1 = np.zeros((Mol.N, Mol.N, len(Mol.Gates)))
    CM2 = np.zeros((Mol.N, Mol.N, len(Mol.Gates)))

    for i, ei in enumerate(Mol.Gates):
        CM1[:,:,i] = Mol.CurrentMatrixDiscrete(0.1, ei)
    CM2 = Mol.CurrentMatrixAnalytic(0.1)

    axTM = np.ndarray((Mol.N, Mol.N), dtype=np.object)
    fig, axCM = plt.subplots(Mol.N, Mol.N)

    for i in range(Mol.N):
        for j in range(Mol.N):
            axCM[i,j].semilogy(Mol.Gates, CM1[i,j,:]*np.conjugate(CM1[i,j,:]), 'r', linewidth = 3, alpha=0.5)
            axCM[i,j].semilogy(Mol.Gates, CM2[i,j,:]*np.conjugate(CM2[i,j,:]), 'g')
    axCM[0,0].legend(('Discrete','Analytic'))
    fig.suptitle("Current Matrix Methods")

def TestCurrentMethods(Mol):
    ISeldenthuis = Mol.CurrentFromSeldenthuis(0.1)
    ImatrixDiscr = Mol.CurrentFromMatrixDiscrete(0.1, 100)
    ImatrixAnal  = Mol.CurrentFromMatrixAnalytic(0.1)

    fig, (Current) = plt.subplots(1,1)
    fig.suptitle("Current Methods")

    Current.semilogy(Mol.Gates, ISeldenthuis, 'r', linewidth = 5, alpha=0.3)
    Current.plot(Mol.Gates, ImatrixDiscr, 'g', linewidth = 3, alpha = 0.6)
    Current.semilogy(Mol.Gates, ImatrixAnal, 'b')
    Current.legend(('Seldenthuis', 'Matrix Discrete', 'Matrix Analytic'))
    return


            
if __name__ == '__main__':
    from Molecule import *
    from AtomParam import cAtom

 ##   Mol1 = Graphene(3,3, 'Armchair')
    Mol1 = Molecule()
    Mol1.AddAtom('C', 0, 0, 0)
    Mol1.AddAtom('C', 1.42, 0, 0)
    Mol1.SetLead (0, "Left")
    Mol1.SetLead (0, "Right")
    


    Graph = Graphics()
    Graph.SetMolecule(Mol1)
    print Mol1.Beta
 #   Graph.RunSequence()

##    TestTransmissionMatrixMethods (Mol1)
##    TestTransmissionMethods(Mol1)
##    TestTransmissionOrbitalMethods(Mol1)
##    TestCurrentMatrixMethods(Mol1)
##    TestCurrentMethods(Mol1)
    plt.show()



    
    
