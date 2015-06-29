import matplotlib.pyplot as plt

import numpy as np
import math as math
from Molecule import *

Mol1 = Graphene(7,21, 'Armchair')
atom = Mol1.Atom[23]
atom = ('B', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (23, atom)
    
atom = Mol1.Atom[27]
atom = ('B', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (27, atom)
    
atom = Mol1.Atom[31]
atom = ('B', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (31, atom)
    
atom = Mol1.Atom[34]
atom = ('N', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (34, atom)
    
atom = Mol1.Atom[38]
atom = ('N', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (38, atom)
    
atom = Mol1.Atom[42]
atom = ('N', atom[1], atom[2], atom[3])
Mol1.UpdateAtom (42, atom)
Mol1.CreateHam()
Mol1.CreateHamExt(0.1, 'Not self consistent')
Mol1.CalcGamma()
Mol1.UpdateGates(-8.2, 8.2)


wc = np.zeros((len(Mol1.Gates)))
(Atoms, x, y, z) = Mol1.GetCoord()
yavg = np.mean(y)
print "Y average:", str(yavg)
for k, ei in enumerate(Mol1.Gates):
    I = Mol1.LocalCurrents(600, ei)
    Imax = abs(I).max()
    print k
    for i in range(Mol1.N):
        for j in range(i+1, Mol1.N):
            Inet = abs(I[i,j]-I[j,i])
            wc[k] += Inet/Imax*abs((y[i]+y[j])/2-yavg)**2
    wc[k] = wc[k]/Mol1.N

fig, (ax1, ax2) = plt.subplots(2,1)
ax1.plot(Mol1.Gates, wc)
ax1.set_xlim([min(Mol1.Gates), max(Mol1.Gates)])
wcreg = [3.2, 7.4, 9.5, 11]
creg = ['tomato', 'gold' ,'darksage', 'teal']

reg = np.zeros((len(Mol1.Gates)))
for i, wci in enumerate(wc):
    delta = np.zeros((len(wcreg)))
    for j in range(len(wcreg)):
        delta[j] = abs(wc[i]-wcreg[j])

    print delta
    reg[i] = int(np.argmin(delta))
    print reg[i]

y = np.zeros((len(Mol1.Gates)))
ax2.scatter(Mol1.Gates,y, s=2, c=[creg[int(i)] for i in reg], lw=0, marker=",")

ax2.set_xlim([min(Mol1.Gates), max(Mol1.Gates)])
    
plt.show()
