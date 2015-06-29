import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from matplotlib.widgets import CheckButtons, RadioButtons, Button, Slider
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math as math
from Molecule import *
from AtomParam import cAtom, gLeft, gRight
import matplotlib.cm as cm
from numpy.random import randn

def DisplayMatrix (Matrix):
    fig, ax = plt.subplots(figsize=(20,11))
    table = plt.table(cellText = np.round(Matrix, 8), cellLoc='center', bbox=[0, 0, 1, 1])
    table.set_fontsize(25)
    ax.axis('off')
    plt.show()

class GraphMolecule (object):
    """ Manages the displaying of the molecule and selected atoms, currents, densities
        Attributes
        axMol                   Handle to the window with the atoms
        axAtoms                 Handle to the individual atoms
        AtomPlot                Handle to the atom indicator spot
        AtomSel                 Index of the selected atom
        LDPlot                  Handle to the local densities
        LCPlot                  Handle to the local current arrows

        Methods
        InitMolecule()
        ChangeMolecule()        Creates handles to the density and selection data points
        DrawMolecule()          Draws the molecule itself
        DrawLocalDensity(bDraw) Draws the local density onto the atom sites
        DrawLocalCurrents(bDraw)Draws the local currents from the atoms to each other
        OnClickMolecule(event)  Unselect atom for double click
        OnPickMolecule(event)   Select clicked atom
        OnPress()               Select next/previous atom or change atoms 
        DrawSelAtom()           Draws the indicator for the selected atom or series of atoms
        """
    
    def InitMolecule (self):
        self.axMol = self.axBias.twinx()
        self.axMol.set_position (self.axBias.get_position())
        self.axAtoms = None
        self.AtomPlot = []
        self.AtomSel = None
        self.LDPlot = []
        self.LCPlot = None

        cmap = mpl.cm.YlGnBu
        norm = mpl.colors.Normalize(-10, vmax=0)
        self.cb1 = mpl.colorbar.ColorbarBase(self.axCM, cmap=cmap, norm=norm, orientation='vertical')
        self.cb1.set_label('Current (10^)')

        
    def ChangeMolecule(self):
        if len(self.LDPlot) is not 0:
            for LDPl in self.LDPlot:
                LDPl.remove()
        del self.LDPlot
        if len(self.AtomPlot) is not 0:
            for AtPl in self.AtomPlot:
                AtPl.remove()
        del self.AtomPlot


        self.LDPlot = []
        self.AtomPlot = []
        for i in range(self.Mol.N):
            self.LDPlot.append(None)
            self.LDPlot[i],     = self.axMol.plot(self.Mol.Atom[i][1],self.Mol.Atom[i][2], 'ko', ms=0, alpha=0.25, visible=True)
            self.AtomPlot.append(None)
            self.AtomPlot[i],   = self.axMol.plot(self.Mol.Atom[i][1],self.Mol.Atom[i][2], 'ro', ms=0, alpha=0.5, visible=True)

        if self.LCPlot is not None:
            for i in range(self.Mol.N):
                for j in range(self.Mol.N):
                    if i!=j:
                        self.LCPlot[i,j].remove()
        del self.LCPlot
        self.LCPlot = np.ndarray((self.Mol.N, self.Mol.N), dtype=np.object)
        for i in range(self.Mol.N):
            for j in range(self.Mol.N):
                if i!=j:
                    x0 = self.Mol.Atom[i][1]
                    y0 = self.Mol.Atom[i][2]
                    dx = self.Mol.Atom[j][1]-x0
                    dy = self.Mol.Atom[j][2]-y0
                    self.LCPlot[i,j] = self.axMol.arrow(x0+0.1*dx, y0+0.1*dy, dx*0.8, dy*0.8, head_width=0.20, head_length=+0.1, linewidth = 30, fc='b', ec='b', alpha=0.5, visible=False)

        
    def DrawMolecule (self):
        Atom,x,y,z = self.Mol.GetCoord()
        N = len(x)
        for i in range(N):
            for j in range(i+1,N):
                if self.Mol.Ham[i,j]!=0:
                    self.axMol.plot ([self.Mol.Atom[i][1], self.Mol.Atom[j][1]], [self.Mol.Atom[i][2], self.Mol.Atom[j][2]], linewidth=1, color='k', alpha=0.8)        

        self.axMol.plot(x[self.Mol.lpL]-0.15, y[self.Mol.lpL], '>y', markersize=7, alpha=0.5)
        self.axMol.plot(x[self.Mol.lpR]+0.15, y[self.Mol.lpR], '<y', markersize=7, alpha=0.5)
        
        if self.axAtoms!=None:
            self.axAtoms.remove()
        self.axAtoms = self.axMol.scatter(x,y, s=60, c=[cAtom[i] for i in Atom], lw=0, picker=True)

        self.axMol.set_xlabel('x (Angstrom)')
        self.axMol.set_ylabel('y (Angstrom)')

        self.axMol.axis('equal')
        self.axMol.relim()
        self.axMol.autoscale_view()
        self.fig.canvas.draw()

    def DrawLocalDensity(self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions1.lines[1][0].get_visible()
            
        if bDraw:
            D = self.Mol.Density(self.Bias)
            Dmax = abs(D).max()
            for i, iLDPlot in enumerate(self.LDPlot):
                iLDPlot.set_ms((abs(self.Mol.LD[i,i])/Dmax)**0.3*30)
                if self.Mol.LD[i,i]>0:
                    iLDPlot.set_color('greenyellow')
                else:
                    iLDPlot.set_color('yellow')
        else:
            for iLDPlot in self.LDPlot:
                iLDPlot.set_ms(0)
        self.fig.canvas.draw()
        
    def DrawLocalCurrents(self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions1.lines[2][0].get_visible()

        self.axCM.set_visible(bDraw)

        if bDraw and self.Gate is not None:
            #I = self.Mol.CurrentMatrix(self.Bias, self.Gate)
            I = self.Mol.LocalCurrents(self.Bias, self.Gate)
            
            #DisplayMatrix(I)
            I = np.real(I)
            maxCurrent = abs(I).max()
            print "Maximum current:", maxCurrent
            for i in range(self.Mol.N):
                for j in range(i+1,self.Mol.N):
                    Inet = (I[i,j]-I[j,i])/2
                    if abs(Inet)>maxCurrent*0.1:
                        color = cm.YlGnBu(1-math.log10(abs(Inet))/-10)
                        if Inet>0:
                            self.LCPlot[i,j].set_linewidth(Inet/maxCurrent*10)
                            self.LCPlot[i,j].set_color(color)
                            self.LCPlot[i,j].set_visible(True)
                            self.LCPlot[j,i].set_visible(False)
                            print "Current from", i, "to", j, "is", I[i,j]
                        else:
                            self.LCPlot[j,i].set_linewidth(-Inet/maxCurrent*10)
                            self.LCPlot[j,i].set_color(color)
                            self.LCPlot[j,i].set_visible(True)
                            self.LCPlot[i,j].set_visible(False)
                            print "Current from", j, "to", i, "is", -I[i,j]
                    else:
                        self.LCPlot[i,j].set_visible(False)
                        self.LCPlot[j,i].set_visible(False)
        else:
            for i in range(self.Mol.N):
                for j in range(self.Mol.N):
                    if i!=j:
                        self.LCPlot[i,j].set_visible(False)
            
        self.fig.canvas.draw()
        
    def OnClickMolecule (self, event):
        if event.dblclick:
            self.UpdateAtomSel(None)

    def OnPickMolecule (self, event):
        self.UpdateAtomSel(event.ind[0])

    def OnPress (self, event):
        if (self.AtomSel is None) or isinstance(self.AtomSel, np.ndarray): return
        if event.key in ('right', 'left'): 
            if event.key=='right': inc = 1
            else:  inc = -1
            self.UpdateAtomSel((self.AtomSel + inc)%len(self.Mol.Atom), True)
        elif cAtom.has_key(event.key.upper()):
            print 'Update atom', self.AtomSel, 'into', event.key.upper()
            self.Mol.UpdateAtom (self.AtomSel, (event.key.upper(), self.Mol.Atom[self.AtomSel][1], self.Mol.Atom[self.AtomSel][2], self.Mol.Atom[self.AtomSel][3]))
            self.UpdateMolecule()

    def DrawSelAtom (self):
        if isinstance(self.AtomSel, np.ndarray):
            for i, ms in enumerate(self.AtomSel):
                if np.real(ms)>=0:
                    self.AtomPlot[i].set_color('r')
                else:
                    self.AtomPlot[i].set_color('b')
                self.AtomPlot[i].set_ms(np.real(ms*np.conjugate(ms))*100)
        else:
            if self.AtomSel is None:
                for AtomPl in self.AtomPlot:
                    AtomPl.set_ms(0)
            else:
                for i, iAtomPlot in enumerate(self.AtomPlot):
                    if i==self.AtomSel:
                        iAtomPlot.set_ms(15)
                        iAtomPlot.set_color('r')
                    else:
                        iAtomPlot.set_ms(0)
        self.fig.canvas.draw()
        

class GraphBias(object):
    """ Manages the displaying and control of the bias
        Attributes
        axBias                  Handle to the bias window

        Methods
        InitBias()
        DrawBias()              Draws the bias
        OnClickBias(event)      Changes the bias based on a click of the right handed side of the mouse
        """
    
    def InitBias (self):
        return
    
    def DrawBias (self):
        self.axBias.cla()
        self.axBias.axhspan(-self.Bias/2, self.Bias, 0.00, 0.05, facecolor='y', alpha=0.5)
        self.axBias.axhspan(-self.Bias/2,         0, 0.95, 1.00, facecolor='y', alpha=0.5)
        self.axBias.set_ylim([-self.Bias/2 ,float(3)/2*self.Bias])
        self.axBias.set_ylabel('Bias [V]', color='y')
        for tl in self.axBias.get_yticklabels():
            tl.set_color('y')
        self.fig.canvas.draw()

    def OnClickBias (self, event):
        if event.dblclick:
            self.UpdateBias(self.Bias*10)
        else:
            [ymin, ymax] = self.axMol.get_ylim()
            self.UpdateBias(-self.Bias/2+2*self.Bias*((event.ydata-ymin)/(ymax-ymin)))
        

class GraphTrans(object):
    """ Manages the displaying of the transmission and current and the control on the gate
        Attributes
        axTrans                 Handle to the transmission window
        GatePlot                Handle to the plot of the gate
        CurrentPlot             Handle to the plot of the current
        TransPlot               Handle to the plot of the transmisssion
        PhasePlot               Handle to the plot of the phase
        DOSPlot                 Handle tot the plot of the density of states

        Methods
        InitTrans()
        DrawTransmission()      Draws the transmission and its phase
        DrawCurrent()           Draws the current
        DrawDOS()               Draws the density of states
        DrawGate()              Draws the selected gate voltage
        OnClickTrans(event)     Changes the gate based on a click with the mouse
        """
    
    def InitTrans (self):
        (self.axTrans, self.axPhase, self.axOrb)  = [self.axTrans, self.axTrans.twinx(), self.axTrans.twinx()]
        self.GatePlot,   = self.axTrans.plot([1, 1],[1, 1], 'r', alpha=0.5, linewidth = 4, visible=False)
        self.CurrentPlot,= self.axTrans.semilogy(1,1, 'b', visible=True)
        self.TransPlot,  = self.axTrans.semilogy(1,1, 'g', visible=True)
        self.PhasePlot,  = self.axPhase.plot(1,1, 'm', visible=True)
        self.DOSPlot,    = self.axTrans.semilogy(1,1, 'darkorange', visible=True)
        self.axTrans.set_xlabel('E-Ef (V)')

    def DrawTransmission (self, bDraw=None, bDrawPhase=None):
        if bDraw==None:
            bDraw = self.cbOptions2.lines[0][0].get_visible()
        if bDrawPhase==None:
            bDrawPhase = self.cbOptions2.lines[4][0].get_visible()
        self.TransPlot.set_visible(bDraw)
        self.PhasePlot.set_visible(bDrawPhase)

        if bDraw:
            self.TransPlot.set_data(self.Mol.Gates, self.Mol.T)
            self.PhasePlot.set_data(self.Mol.Gates, self.Mol.P)
        self.axTrans.relim()
        self.axTrans.autoscale_view()
        self.axTrans.set_xlim([min(self.Mol.Gates), max(self.Mol.Gates)])

##        T = np.zeros(1000, dtype=np.dtype('c16'))
##        for i, e in enumerate(self.Mol.Gates):
##            T[i] = np.trace(self.Mol.CurrentMatrix(e, self.Bias), i==50)
##        self.axTrans.semilogy(self.Mol.Gates, np.real(T), 'purple', linewidth=6, alpha=0.3)
        self.fig.canvas.draw()
        #plt.legend(['Current', 'Transmission'])

    def DrawCurrent (self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions2.lines[1][0].get_visible()
        self.CurrentPlot.set_visible(bDraw)
        
        if bDraw:            
            self.CurrentPlot.set_data(self.Mol.Gates, self.Mol.I/self.Bias)
        self.axTrans.relim()
        self.axTrans.autoscale_view()
        self.axTrans.set_xlim([min(self.Mol.Gates), max(self.Mol.Gates)])
        self.fig.canvas.draw()

    def DrawDOS (self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions2.lines[2][0].get_visible()
        self.DOSPlot.set_visible(bDraw)

        if bDraw:
            self.DOSPlot.set_data(self.Mol.Gates, np.real(self.Mol.DOS))
        self.axTrans.relim()
        self.axTrans.autoscale_view()
        self.axTrans.set_xlim([min(self.Mol.Gates), max(self.Mol.Gates)])
        self.fig.canvas.draw()
        
    def DrawGate (self):
        if self.Gate == None or isinstance(self.OrbitalSel, int):
            self.GatePlot.set_visible(False)
        else:
            self.GatePlot.set_visible(True)
            self.GatePlot.set_data([self.Gate, self.Gate], [self.axTrans.get_ylim()])
        self.fig.canvas.draw()
        
    def OnClickTrans (self, event):
        if event.dblclick:
            self.UpdateGate(None)
        else:
            self.UpdateGate(event.xdata)
            
        
class GraphOrbitals (object):
    """ Manages the displaying of the eigenstates(orbitals?) of the molecule at the eigenvalues/energies
        Attributes
        axPhase                 Handle to the phase window
        axOrb                   Handle to the orbitals window
        plOrb                   Handle to the orbital plots
        rcOrb                   Handle to the rectangles of the orbital plots
        plSelOrb                Handle to the selected orbital plots
        rcSelOrb                handle to the rectangles of the selected orbital plots
        SelPhasePlot            Handle to the phase plot of the selected orbital
        SelTransPlot            Handle to the transmission plot of the selected orbital

        Methods
        InitTrans()
        ChangeMoleculeOrbitals  Creates a new number of handles for the orbitals and selected orbitals when the molecule changes
        DrawOrbitals()          Draws the eigenstates/orbitals at the appropriate eigenvalues/energies
        OnPickOrbital(event)    Changes the selected eigenstate/orbital
        DrawSelOrbitals()       Draws the indicators for the selected orbitals or series of orbitals
        DrawSelTransmission()   Draws the transmission and its phase from the selected orbital
        """

    def InitOrbitals (self):
        self.axPhase.set_position (self.axTrans.get_position())
        self.axOrb.set_position (self.axTrans.get_position())
        self.axOrb.spines['right'].set_position(('axes', 1.07))
        self.axOrb.tick_params(axis='y', colors='lightsteelblue')
        self.axOrb.set_ylabel("Eigenenergies of transmission", color='lightsteelblue')
        #self.axOrb.set_ylim(0, 1.2)
        self.axPhase.tick_params(axis='y', colors='m')
        self.axPhase.set_ylabel("Phase", color='m')
        self.axPhase.set_ylim(-math.pi, math.pi)
        self.plOrb = []
        self.rcOrb = []
        self.plSelOrb = []
        self.rcSelOrb = []
        self.SelPhasePlot, = self.axPhase.plot(1,1, 'm--', visible=True)
        self.SelTransPlot, = self.axTrans.semilogy(1,1, 'g--', visible=True)
##        self.SelTransPlot2, = self.axTrans.semilogy(1,1, 'b*', visible=True) #temporary for testing

    def ChangeMoleculeOrbitals(self):
        for i in range(len(self.plOrb)):
            self.rcOrb[i].remove()
            self.rcSelOrb[i].remove()
            
        del self.plOrb, self.rcOrb, self.plSelOrb, self.rcSelOrb
        self.plOrb = []
        self.rcOrb = []
        self.plSelOrb = []
        self.rcSelOrb = []
        
        for i in range(self.Mol.N):
            self.plOrb.append(None)
            self.plOrb[i] = self.axOrb.bar(0, 0, color='lightsteelblue', width=0.1, linewidth=1, alpha=0.4, picker=0.0)
            for rect in self.plOrb[i]:
                self.rcOrb.append(None)
                self.rcOrb[i] = rect
            self.plSelOrb.append(None)
            self.plSelOrb[i] = self.axOrb.bar(0, 0, color='r', width=0.1, linewidth=0, alpha=0.4)
            for rect in self.plSelOrb[i]:
                self.rcSelOrb.append(None)
                self.rcSelOrb[i] = rect
        
    def DrawOrbitals (self, bDraw=None):
        ymax = 0
        if bDraw==None:
            bDraw = self.cbOptions2.lines[3][0].get_visible()
        if bDraw:
            width = self.rcOrb[1].get_width()
            for i in range(self.Mol.N):
                if i==0:
                    y=0
                elif abs(self.Mol.e_arr[i-1]-self.Mol.e_arr[i])<0.01:
                    y+=0
                else:
                    y=0
                T = self.Mol.TransmissionMatrix(self.Mol.e_arr[i])
                height = min(1, max(0.1, abs(np.real(T[i,i]))))
                self.rcOrb[i].set_xy((self.Mol.e_arr[i]-width/2, y))
                self.rcOrb[i].set_height(height)
                self.rcSelOrb[i].set_xy((self.Mol.e_arr[i]-width/2, y))
                self.rcSelOrb[i].set_height(0)
                y+=height
                ymax = max(ymax, y)
        else:
            for i in range(self.Mol.N):
                self.rcOrb[i].set_xy((0,0))
                self.rcOrb[i].set_height(0)
            self.OrbitalSel = None

        self.axOrb.set_ylim([0, ymax*1.1])
        self.fig.canvas.draw()

    def OnPickOrbital (self, event):
        for i, rect in enumerate(self.rcOrb):
            if rect==event.artist:
                self.UpdateOrbitalSel(i)

    def DrawSelOrbitals (self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions2.lines[3][0].get_visible()
        
        if isinstance(self.OrbitalSel, list) and bDraw:
            for i, height in enumerate(self.OrbitalSel):
                self.rcSelOrb[i].set_height(abs(height)*self.rcOrb[i].get_height())
        else:
            for i in range(self.Mol.N):
                self.rcSelOrb[i].set_height(0)
            if self.OrbitalSel!=None and bDraw:
                self.rcSelOrb[self.OrbitalSel].set_height(self.rcOrb[self.OrbitalSel].get_height())
                
        self.fig.canvas.draw()

    def DrawSelTransmission (self, bDraw=None):
        if bDraw==None:
            bDraw = self.cbOptions2.lines[4][0].get_visible()
        if isinstance(self.OrbitalSel, list) or self.OrbitalSel==None:
            self.SelPhasePlot.set_visible(False)
            self.SelTransPlot.set_visible(False)
        else:
            T, Trans, Phase = self.Mol.TransmissionOrbital(self.OrbitalSel)
            self.SelTransPlot.set_data(self.Mol.Gates, Trans)
            self.SelPhasePlot.set_data(self.Mol.Gates, Phase)
            self.SelPhasePlot.set_visible(bDraw)
            self.SelTransPlot.set_visible(self.cbOptions2.lines[0][0].get_visible())
        self.axTrans.relim()
        self.axTrans.autoscale_view()
        self.fig.canvas.draw()

class GraphIVCurve (object):
    """ Manages the displaying of the eigenstates(orbitals?) of the molecule at the eigenvalues/energies
        Attributes
        axIVCurve               Handle to the IV Curve window
        Biases                  Array of biases to be calculated
        IVPlot                  Handle to the plot of the IVCurve

        Methods
        InitIVCurve()
        DrawIVCurve()           Draws the current against the applied voltage
        """
    def InitIVCurve (self):
        self.Biases = np.linspace(-2, 2, 100)
        self.IVPlot,  = self.axIVCurve.plot(1,1, 'b', visible=True)
        self.axIVCurve.set_ylabel("Current (A)")
        self.axIVCurve.set_xlabel("Bias (V)")

    def DrawIVCurve (self):
        current = np.zeros(len(self.Biases),np.dtype('c16'))
        if self.Gate==None:
            self.IVPlot.set_visible(False)
            self.fig.canvas.draw()
            return
        
        self.IVPlot.set_visible(True)
        idx = int((self.Gate-self.Mol.Gates[0])/(self.Mol.Gates[len(self.Mol.Gates)-1]-self.Mol.Gates[0])*len(self.Mol.Gates))

        for i, bias in enumerate(self.Biases):
            currents = self.Mol.Current(bias)
            current[i] = np.real(currents[idx])
        self.IVPlot.set_data(self.Biases, np.real(current))
        self.axIVCurve.relim()
        self.axIVCurve.autoscale_view()
        self.fig.canvas.draw()

        
class Graphics (GraphMolecule, GraphBias, GraphTrans, GraphOrbitals, GraphIVCurve):
    """ Manages the graphical representation of the project
        Attributes
        fig                     Handle to the entire figure
        Bias                    The selected bias
        Gate                    The selected gate voltage
        OrbitalSel              The selected orbital or series of orbitals

        axSC                    Handle to the radio buttons window for the self consistency
        cbSC                    Handle to the radio buttons for the self consistency

        axOptions1              Handle to the checkbox window with options 1
        cbOptions1              Handle to the checkboxes with options 1
        
        axOptions2              Handle to the checkbox window with options 2
        cbOptions2              Handle to the checkboxes with options 2

        axGleft, axGright       Handle to the slider windows for selecting the lead interaction strength
        sGleft,  sGright        Handle to the sliders for selecting the lead interaction strength

        axSave                  Handle to the save button window
        bSave                   Handle to the save button

        Methods
        init()
        OnPick(event)           Manages the pick events
        OnClick(event)          Manages the on click events
        Save()                  Manages the input from the save button
        Options1(label)         Manages the input from the options 1 window
        Options2(label)         Manages the input from the options 2 window

        SetMolecule(Mol)        Sets the molecule class from which the graphics gets its information
        UpdateMolecule()        Updates everything that changes when the molecule has changed

        UpdateConsistency(label)Updates the selected consistency method
        UpdateG(val)            Updates the interaction strength with the leads
        UpdateBias(bias)        Updates the selected bias
        UpdateHamExt()          Updates everything that changes when the extended hamiltonian is changed
        
        UpdateGate(gate)        Updates the selected gate voltage
        UpdateAtomSel(iAtom)    Updates the selected atom
        UpdateOrbitalSel(iOrb)  Updates the selected orbital
        UpdateSelection()       Updates everything that changes after one of the selections has changed

    """

    def __init__(self):
        self.Consistency = 'Not self consistent'
        self.Bias = 0.1
        self.Gate = None
        self.OrbitalSel = None
        self.fig, (self.axBias, self.axTrans, self.axIVCurve) = plt.subplots(3,1)
        self.fig.patch.set_facecolor('ghostwhite')
        plt.subplots_adjust(left=0.3)
        pos = self.axBias.get_position()
        self.axBias.set_position ([pos.x0, 0.55, pos.width, 0.40])
        self.axTrans.set_position([pos.x0, 0.05, pos.width, 0.40])
        self.axIVCurve.set_position([0.05, 0.05, 0.2, 0.3])
        self.axCM = self.fig.add_axes([0.94, 0.55, 0.01, 0.35])

        self.fig.canvas.mpl_connect('button_press_event', self.OnClick)
        self.fig.canvas.mpl_connect('pick_event', self.OnPick)
        self.fig.canvas.mpl_connect('key_press_event', self.OnPress)

        self.InitMolecule()
        self.InitBias()
        self.InitTrans()
        self.InitOrbitals()
        self.InitIVCurve()
        
        self.axSC = plt.axes([0.05, 0.85, 0.15, 0.10], axisbg='white')
        self.cbSC = RadioButtons(self.axSC, ('Not self consistent', 'Hubbard', 'PPP'))
        self.cbSC.on_clicked(self.UpdateConsistency)

        self.axOptions1 = plt.axes([0.05, 0.7, 0.15, 0.10], axisbg='white')
        self.cbOptions1 = CheckButtons(self.axOptions1, ('Overlap', 'Show Local Density','Show Local Currents'), (False, False, True))
        self.cbOptions1.on_clicked(self.Options1)
        
        self.axOptions2 = plt.axes([0.05, 0.5, 0.15, 0.15], axisbg='white')
        self.cbOptions2 = CheckButtons(self.axOptions2, ('Show Transmission', 'Show Current', 'Show DOS', 'Show Orbitals', 'Show Phase'), (True, True, False, False, False))
        self.cbOptions2.on_clicked(self.Options2)
        c = ['seagreen', 'b', 'darkorange', 'lightsteelblue', 'm']    
        [rec.set_facecolor(c[i]) for i, rec in enumerate(self.cbOptions2.rectangles)]
        
        self.axGleft  = plt.axes([0.05, 0.43, 0.15, 0.02], axisbg='white')
        self.axGright = plt.axes([0.05, 0.40, 0.15, 0.02], axisbg='white')
        self.sGleft   = Slider(self.axGleft,  'Gl', 0.0, 1.0, valinit = gLeft)
        self.sGright  = Slider(self.axGright, 'Gr', 0.0, 1.0, valinit = gRight)
        self.sGleft.on_changed(self.UpdateG)
        self.sGright.on_changed(self.UpdateG)
        
        self.axSave = plt.axes([0.92, 0.95, 0.07, 0.04])
        self.bSave = Button(self.axSave, 'Save')
        self.bSave.on_clicked(self.Save)

    def OnPick(self, event):
        if isinstance(event.artist, Rectangle):
            self.OnPickOrbital(event)
        else:
            self.OnPickMolecule(event)
                   
    def OnClick (self, event):
        if event.inaxes==self.axMol:
            if event.button==1:
                self.OnClickMolecule (event)
            elif event.button==3:
                self.OnClickBias (event)
        elif event.inaxes==self.axOrb:
            if event.button==1:
                return
            if event.button==3:
                self.OnClickTrans (event)
        return

    def Save(self, event):
        self.Mol.Save()

    def Options1(self, label):
        if label == 'Overlap':
            self.Mol.SetOverlap(self.cbOptions1.lines[0][0].get_visible())
            self.UpdateMolecule()
        elif label == 'Show Local Density':
            self.DrawLocalDensity(self.cbOptions1.lines[1][0].get_visible())
        elif label == 'Show Local Currents':
            self.DrawLocalCurrents(self.cbOptions1.lines[2][0].get_visible())
        return

    def Options2(self, label):
        if label == 'Show Transmission':
            self.DrawTransmission()
        elif label == 'Show Current':
            self.DrawCurrent()
        elif label == 'Show DOS':
            self.DrawDOS()
        elif label == 'Show Orbitals':
            self.DrawOrbitals()
            self.DrawSelOrbitals()
        elif label == 'Show Phase':
            self.DrawTransmission()
            self.DrawSelTransmission()
        return

     
    def SetMolecule(self, Molecule):
        self.Mol = Molecule
        self.UpdateMolecule()
        
    def UpdateMolecule (self):
        self.Mol.CreateHam()
        self.Mol.CreateOverlap()
        self.Mol.CreateV()
        self.ChangeMolecule ()
        self.ChangeMoleculeOrbitals ()
        self.DrawMolecule ()
        return self.UpdateHamExt ()
        

    def UpdateConsistency(self, label):
        self.Consistency = label
        print "Consistency set to:", self.Consistency
        return self.UpdateHamExt ()     
       
    def UpdateG (self, val):
        print "Interaction strengths set to:", self.sGleft.val, self.sGright.val
        self.Mol.SetG (self.sGleft.val, self.sGright.val)
        return self.UpdateHamExt ()
    
    def UpdateBias(self, bias):
        self.Bias = bias
        print "Bias voltage set to: ", self.Bias, "V"
        return self.UpdateHamExt ()

    def UpdateHamExt (self):
        self.Mol.CreateHamExt (self.Bias, self.Consistency)
        self.Mol.CalcGamma()
        self.Mol.Density(self.Bias)
        self.Mol.Transmission()
        self.Mol.Current(self.Bias)
        self.Mol.CalcDOS()
        
        self.DrawLocalDensity ()
        self.DrawBias ()
        self.DrawTransmission ()
        self.DrawDOS ()
        self.DrawCurrent ()
        self.DrawOrbitals ()
        return self.UpdateSelection()

    def UpdateGate (self, gate):
        self.Gate = gate
        print "Gates voltage set to: ", self.Gate, "V"
        self.OrbitalSel = None
        self.AtomSel    = None
        return self.UpdateSelection()

    def UpdateAtomSel(self, iAtom):
        self.AtomSel = iAtom
        if self.AtomSel == None:
            print "No atom selected"
            self.axMol.set_title('No atom selected')
            self.OrbitalSel = None
        elif isinstance(iAtom, int):
            print "Selected atom", self.Mol.Atom[self.AtomSel][0], "at", self.AtomSel, " Local density = ", self.Mol.LD[self.AtomSel, self.AtomSel]
            self.axMol.set_title('Selected atom: %c at %d. Density = %f'%(self.Mol.Atom[self.AtomSel][0],self.AtomSel, self.Mol.LD[self.AtomSel, self.AtomSel]))
            Orbitals = []
            for i in range(self.Mol.N):
                Orbitals.append(np.real(self.Mol.eigvec[i][self.AtomSel]*np.conjugate(self.Mol.eigvec[i][self.AtomSel])))
            self.OrbitalSel = Orbitals
        return self.UpdateSelection()
        
    def UpdateOrbitalSel(self, iOrbital):
        self.OrbitalSel = iOrbital
        
        if isinstance(self.OrbitalSel, int):
            print "Orbital set to:", self.OrbitalSel, "with energy", self.Mol.e_arr[self.OrbitalSel], "eV"
            self.AtomSel    = self.Mol.eigvec[iOrbital]
            self.Gate       = self.Mol.e_arr[self.OrbitalSel]
        return self.UpdateSelection()

    def UpdateSelection (self):
        self.DrawLocalCurrents ()
        self.DrawGate()
        self.DrawIVCurve()
        self.DrawSelAtom()
        self.DrawSelOrbitals()
        self.DrawSelTransmission()
        return

    def RunSequence (self):
        self.fig.set_size_inches(18.5, 13.5) #default 18.5, 10.5
        
        N=100
        for i in range(N):
            e = self.Mol.Gates[0]+(self.Mol.Gates[-1]-self.Mol.Gates[0])/N*i
            self.UpdateGate(e)            
            self.fig.savefig('Output/PN/Armchair 7,21, seq ' + str(i) + ' Energy=' + str(math.ceil(e*100)/100) + 'eV.png', dpi=self.fig.dpi)
            print e
            
if __name__ == '__main__':

##    Mol1.Load("molecules/L_Bend_Graphene")
##    Mol1.SetLeads((-1, 0))
##    Mol1.SetLead(136, 'Left')
##    Mol1.SetLead(137, 'Left')
##    Mol1.SetLead(138, 'Left')
##    Mol1.SetLead(139, 'Left')

    
    Mol1 = Graphene(7,11, 'Armchair')
    

##    Mol1 = Graphene(7,11, 'Armchair')
##    atom = Mol1.Atom[14]
##    atom = ('N', atom[1], atom[2], atom[3])
##    Mol1.UpdateAtom (14, atom)
##    atom = Mol1.Atom[18]
##    atom = ('N', atom[1], atom[2], atom[3])
##    Mol1.UpdateAtom (18, atom)
##    atom = Mol1.Atom[22]
##    atom = ('N', atom[1], atom[2], atom[3])
##    Mol1.UpdateAtom (22, atom)
    
    Graph = Graphics()
    Graph.SetMolecule(Mol1)
   # Graph.RunSequence()
    
    
    plt.show()
