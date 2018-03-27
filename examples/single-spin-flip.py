# ****************************************************************************
# 
# This tutorial shows how to simply generate the input files, start the task in python (not recommended on a cluster) and then evaluate. It is designed to the dipolar XY system in 2D, where the order parameter is a staggered magnetic field (stripe like)
#
# To a big proportion this tutorial is written by Brigitte Surer <surerb@phys.ethz.ch> in the spinmc part of the ALPS library. As it is similar the code was reused. 
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
from numpy import linspace

#prepare the input parameters
parms = []
for l in [8,16]: 
    for t in linspace(0,1.5,16):
        parms.append(
            { 
                 'LATTICE'        : "square lattice", 
                 'T'              : t,
                 'D'              : 1.,
                 'THERMALIZATION' : 50000,
                 'SWEEPS'         : 100000,
                 'ALGORITHM'      : "xy",
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 3.,
                 'L'              : l
               }
           )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm',parms)
#pyalps.runApplication('mc++',input_file,Tmin=5)
# use the following instead if you have MPI
pyalps.runApplication('mc++',input_file,Tmin=5,MPI=7)

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),['M staggered', 'c_V', 'BinderCumulant staggered', 'susceptibility staggered'])
M = pyalps.collectXY(data,x='T',y='M staggered',foreach=['L'])
chi = pyalps.collectXY(data,x='T',y='susceptibility staggered',foreach=['L'])
c_v= pyalps.collectXY(data,x='T',y='c_V',foreach=['L'])
binder = pyalps.collectXY(data,x='T',y='BinderCumulant staggered',foreach=['L'])

for mag in M:
    mag.props['label']='L='+str(mag.props['L'])
for c in chi:
    c.props['label']='L='+str(c.props['L'])
for c in c_v:
    c.props['label']='L='+str(c.props['L'])
for b in binder:
    b.props['label']='L='+str(b.props['L'])


#make plots
plt.figure()
pyalps.plot.plot(M)
plt.xlabel('Temperature $T$')
plt.ylabel('staggered Magnetization $M$')
plt.title('2D XY model with dipolar interaction')

plt.figure()
pyalps.plot.plot(chi)
plt.xlabel('Temperature $T$')
plt.ylabel('staggered Susceptibility $\chi$')
plt.title('2D XY model with dipolar interaction')

plt.figure()
pyalps.plot.plot(c_v)
plt.xlabel('Temperature $T$')
plt.ylabel('Specific Heat $c_V$')
plt.title('2D XY model with dipolar interaction')

plt.figure()
pyalps.plot.plot(binder)
plt.xlabel('Temperature $T$')
plt.ylabel('Binder Cumulant')
plt.title('2D XY model with dipolar interaction')

plt.show()
