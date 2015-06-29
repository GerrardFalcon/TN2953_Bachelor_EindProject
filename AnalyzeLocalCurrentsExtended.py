import matplotlib.pyplot as plt

import numpy as np
import math as math
from Molecule import *

Mol1 = Graphene(7,11, 'Armchair')
atom = Mol1.Atom[17]
atom = ('N', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (17, atom)
Mol1.CreateHam()
Mol1.CreateHamExt(0.1, 'Not self consistent')
Mol1.CalcGamma()


wcy = np.zeros((len(Mol1.Gates)))
wcx = np.zeros((len(Mol1.Gates)))
(Atoms, x, y, z) = Mol1.GetCoord()
yavg = np.mean(y)
print "Y average:", str(yavg)
xavg = np.mean(x)
print "X average:", str(xavg)
for k, ei in enumerate(Mol1.Gates):
    I = Mol1.LocalCurrents(600, ei)
    Imax = abs(I).max()
    print k
    for i in range(Mol1.N):
        for j in range(i+1, Mol1.N):
            Inet = abs(I[i,j]-I[j,i])
            wcy[k] += Inet/Imax*abs((y[i]+y[j])/2-yavg)**2
            wcx[k] += Inet/Imax*abs((x[i]+x[j])/2-xavg)**2
    wcy[k] = wcy[k]/Mol1.N
    wcx[k] = wcx[k]/Mol1.N

fig, ( ax1, ax2) = plt.subplots(2,1)
ax1.plot(Mol1.Gates, wcy)
ax12 = ax1.twinx()
ax12.plot(Mol1.Gates, wcx, 'r')
ax1.set_xlim([min(Mol1.Gates), max(Mol1.Gates)])
ax12.set_xlim([min(Mol1.Gates), max(Mol1.Gates)])
##      
##wcreg = [5.0, 12.8, 19,5]
##creg = ['tomato', 'gold', 'teal']
##
##reg = np.zeros((len(Mol1.Gates)))
##for i, wci in enumerate(wc):
##    delta = np.zeros((len(wcreg)))
##    for j in range(len(wcreg)):
##        delta[j] = abs(wc[i]-wcreg[j])
##
##    print delta
##    reg[i] = int(np.argmin(delta))
##    print reg[i]
##
##y = np.zeros((len(Mol1.Gates)))
##ax2.scatter(Mol1.Gates,y, s=2, c=[creg[int(i)] for i in reg], lw=0, marker=",")
##
##ax2.set_xlim([min(Mol1.Gates), max(Mol1.Gates)])
    
plt.show()
