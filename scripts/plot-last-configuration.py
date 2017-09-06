# ****************************************************************************
#
# Plot the last configuration encountered in the simulation of a specific file
#
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import numpy as np
import sys

def add_arrow(ax, y,x,size,c,theta=0):
    ax.arrow(x-size/2*np.cos(theta),y-size/2*np.sin(theta),size*np.cos(theta),size*np.sin(theta),width=0.01*size,length_includes_head=True,color=c)

def plot_spins(spins, a=1., b=1.):
    plt.figure()
    ax=plt.gca()
    N=len(spins)
    L=int(np.sqrt(N))
    for x in range(L):
        for y in range(L):
            add_arrow(ax,x*a,y*b,0.5,'#000000',spins[y*L+x])
    ax.set_xlim(-1,(L-1)*b+1)
    ax.set_ylim(-1,(L-1)*a+1)
    plt.show()

EXMC=False
SSF =True
if ('--exmc' in sys.argv):
    EXMC=True
    SSF =False
if (len(sys.argv)<2 or '-h' in sys.argv or '--help' in sys.argv):
    print("Usage: alpspython "+str(sys.argv[0])+" parm.taskX.out.h5")
else:
    a=1.
    b=1.
    try:
        a=pyalps.loadProperties([sys.argv[-1]])[0]["a"]
    except KeyError:
        pass
    try:
        b=pyalps.loadProperties([sys.argv[-1]])[0]["b"]
    except KeyError:
        pass
    if(pyalps.loadProperties([sys.argv[-1]])[0]["LATTICE"]!="square lattice" and
       pyalps.loadProperties([sys.argv[-1]])[0]["LATTICE"]!="anisotropic square lattice"):
        print("Error: at the moment this script only works for square lattice configurations")
        exit()
    if SSF:
        data=pyalps.ResultsToXY(pyalps.loadMeasurements([sys.argv[-1]],['Last Configuration']),x='EXMC: Temperature',y='Last Configuration')[0]
        spins=data.y
    if EXMC:
        ar=pyalps.h5.archive(sys.argv[-1],'r')
        data= ar['simulation/results']
        T=data['EXMC: Temperature']['mean']['value']
        spins=data['Last Configuration']['mean']['value']
    plot_spins(spins, a,b)
