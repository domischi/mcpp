# ****************************************************************************
#
# Plot the last configuration encountered in the simulation of a specific file
#
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import numpy as np
import sys
import re


def add_arrow(ax, x,y,size,c,theta=0):
    ax.arrow(x-size/2*np.cos(theta),y-size/2*np.sin(theta),size*np.cos(theta),size*np.sin(theta),width=0.1*size,length_includes_head=True,color=c)

def plot_spins(T, spins, coordinates, is_deleted):
    plt.figure()
    ax=plt.gca()
    for i in range(len(coordinates)/2):#only 2D
        x=coordinates[2*i]
        y=coordinates[2*i+1]
        if not(bool(is_deleted[i])):
            add_arrow(ax,x,y,0.5,'#000000',spins[i])
    w=max(coordinates)-min(coordinates)
    ax.set_xlim(min(coordinates)-w/10, max(coordinates)+w/10)
    ax.set_ylim(min(coordinates)-w/10, max(coordinates)+w/10)
    plt.title('T='+str(T))
    if not SAVE:
        plt.show()
    else:
        plt.savefig('last-configuration.'+SAVENAME+'.png')
if (len(sys.argv)<2 or '-h' in sys.argv or '--help' in sys.argv):
    print("""Usage: alpspython """+str(sys.argv[0])+""" parm.task[X].out.h5
            --help shows this message
            --save Save the file as a png
            --ssf-clones switches for multiple-clone mode, the usage changes slightly to:
                alpspython """+str(sys.argv[0])+"""--ssf-clones parm.task[X].clone[Y].h5
            --exmc switches for exmc mode, the usage changes slightly to:
                alpspython """+str(sys.argv[0])+"""--exmc parm.task[X].clone[Y].h5 [sector]""")
else:
    FILENAME=''
    SAVE=False
    SSF =True
    EXMC=False
    SSF_CLONES=False
    if ('--ssf-clones' in sys.argv):
        SSF_CLONES=True
        SSF=False
    if ('--exmc' in sys.argv):
        EXMC=True
        SSF =False
    if SSF_CLONES:
        FILENAME=sys.argv[-1]
        ar=pyalps.h5.archive(FILENAME,'r')
        try:
            CLONE_NR=re.match('.*clone(.*)\..*',sys.argv[-1]).group(1)
        except:
            print('Regex failed...')
            exit()
        data= ar['simulation/realizations/0/clones/'+str(CLONE_NR)+'/results/']
        spins     =data['Last Configuration']['mean']['value']
        coords    =data['Coordinates']['mean']['value']
        is_deleted=data['Is Deleted']['mean']['value']
        T=ar['parameters/T']

    if SSF:
        FILENAME=sys.argv[-1]
        data=pyalps.ResultsToXY(pyalps.loadMeasurements([FILENAME],['Last Configuration']),x='EXMC: Temperature',y='Last Configuration')[0]
        spins=data.y
    if EXMC:
        FILENAME=sys.argv[-2]
        ar=pyalps.h5.archive(FILENAME,'r')
        data= ar['simulation/realizations/0/clones/0/results/sections/'+sys.argv[-1]]
        T         =data['EXMC: Temperature']['mean']['value']
        spins     =data['Last Configuration']['mean']['value']
        coords    =data['Coordinates']['mean']['value']
        is_deleted=data['Is Deleted']['mean']['value']
    if ('--save' in sys.argv):
        SAVE=True
        SAVENAME=FILENAME[:-3]
        if(EXMC):
            SAVENAME=SAVENAME+'.sector'+((sys.argv[-1]).rjust(4,str(0)))

    plot_spins(T,spins, coords, is_deleted)
