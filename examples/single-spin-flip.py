# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Brigitte Surer <surerb@phys.ethz.ch> 
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
from numpy import linspace
#prepare the input parameters
parms = []
for l in [4,8,16]: 
    for t in linspace(0,2.5,25):
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 10000,
              'SWEEPS'         : 4000000,
              'UPDATE'         : "ssf",
              'L'              : l
            }
        )
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm',parms)
pyalps.runApplication('mc++',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('mc++',input_file,Tmin=5,MPI=2)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm'))

##load the susceptibility and collect it as function of temperature T
#data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7a'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
#magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
#connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
#spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
#binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
#binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])
#
##make plots
#plt.figure()
#pyalps.plot.plot(magnetization_abs)
#plt.xlabel('Temperature $T$')
#plt.ylabel('Magnetization $|m|$')
#plt.title('2D Ising model')
#
#plt.figure()
#pyalps.plot.plot(connected_susc)
#plt.xlabel('Temperature $T$')
#plt.ylabel('Connected Susceptibility $\chi_c$')
#plt.title('2D Ising model')
#
#plt.figure()
#pyalps.plot.plot(spec_heat)
#plt.xlabel('Temperature $T$')
#plt.ylabel('Specific Heat $c_v$')
#plt.title('2D Ising model')
#
#plt.figure()
#pyalps.plot.plot(binder_u4)
#plt.xlabel('Temperature $T$')
#plt.ylabel('Binder Cumulant U4 $g$')
#plt.title('2D Ising model')
#
#plt.figure()
#pyalps.plot.plot(binder_u2)
#plt.xlabel('Temperature $T$')
#plt.ylabel('Binder Cumulant U2 $g$')
#plt.title('2D Ising model')
#plt.show()

