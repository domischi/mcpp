#**********************************************************************#
# Allows for Finite Size Scaling by using a Manipulate similar feature #
#**********************************************************************#

import os
import pyalps
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import pyalps.plot
from numpy import linspace, log

for l in range(12,36,4): #12,16,...,40 
    for t in linspace(0,1.4,14):
        parms.append(
            { 
                 'LATTICE'        : "square lattice", 
                 'T'              : t,
                 'D'              : 1.,
                 'THERMALIZATION' : 50000,
                 'SWEEPS'         : 100000,
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 1.8,
                 'L'              : l,
                 'Each_Measurement': 15
               }
           )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm',parms)
#pyalps.runApplication('mc++',input_file,Tmin=5)
# use the following instead if you have MPI
pyalps.runApplication('mc++',input_file,Tmin=5,MPI=1)

def f_alpha(Le, Lo, b=2, d=2):
    return 2-d*log(b)/log(Le)
def f_beta(Le, Lo, b=2, d=2):
    return log(b)/log(Le)*(log(Lo)/log(b)-2)
def f_gamma(Le, Lo, b=2, d=2):
    return log(b)/log(Le)*(2*log(Lo)/log(b)-d)
def f_nu(Le, Lo, b=2, d=2):
    return log(b)/log(Le)

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),['M staggered', 'c_V', 'BinderCumulant staggered', 'susceptibility staggered'])

Tc0=0.82
Tc_min=0.7
Tc_max=0.9
le0=1.78
le_min=1.5
le_max=2.0
lo0=1.78
lo_min=1.5
lo_max=2.0
alpha0=f_alpha(le0,lo0)
beta0 =f_beta(le0,lo0)
gamma0=f_gamma(le0,lo0)
nu0   =f_nu(le0,lo0)

left_start=0.05
right_end=1-left_start
bottom_start=0.3
top_end=1
boarder_width=0.1
width=0.5-boarder_width
height=0.35-boarder_width

fig = plt.subplots()
plt.subplots_adjust(left=left_start, bottom=bottom_start, right = width-boarder_width, top=bottom_start+height-boarder_width)
axM= plt.axes([left_start, bottom_start,width,height])
axChi = plt.axes([left_start, bottom_start+height+boarder_width,width,height])
axBinder = plt.axes([left_start+width+boarder_width, bottom_start,width,height])
axCv = plt.axes([left_start+width+boarder_width,  bottom_start+height+boarder_width, width,height]) 
plt.sca(axChi)
chi = pyalps.collectXY(data,x='T',y='susceptibility staggered',foreach=['L'])
for d in chi:
    d.props['label']='L='+str(d.props['L'])
    d.x -= Tc0
    d.x = d.x/Tc0
    l = d.props['L']
    d.x = d.x * pow(float(l),1/nu0)
    d.y = d.y * pow(float(l),-gamma0/nu0)
pyalps.plot.plot(chi)
plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
plt.ylabel(r'$\chi L^{\frac{-\gamma}{\nu}}$', fontsize='x-large')
plt.title(r'$\chi$', fontsize='x-large')
binder = pyalps.collectXY(data,x='T',y='BinderCumulant staggered',foreach=['L'])
plt.sca(axBinder)
for d in binder:
    d.props['label']='L='+str(d.props['L'])
    d.x -= Tc0
    d.x = d.x/Tc0
    l = d.props['L']
    d.x = d.x * pow(float(l),1/nu0)
plt.cla()
pyalps.plot.plot(binder)
plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
plt.ylabel(r'$U^4$', fontsize='x-large')
plt.title('Binder', fontsize='x-large')
M = pyalps.collectXY(data,x='T',y='M staggered',foreach=['L'])
plt.sca(axM)
for d in M:
    d.props['label']='L='+str(d.props['L'])
    d.x -= Tc0
    d.x = d.x/Tc0
    l = d.props['L']
    d.x = d.x * pow(float(l),1/nu0)
    d.y = d.y * pow(float(l),beta0/nu0)
plt.cla()
pyalps.plot.plot(M)
plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
plt.ylabel(r'$ML^{\frac{\beta}{\nu}}$', fontsize='x-large')
plt.title(r'$M$', fontsize='x-large')
M = pyalps.collectXY(data,x='T',y='c_V',foreach=['L'])
plt.sca(axCv)
for d in M:
    d.props['label']='L='+str(d.props['L'])
    d.x -= Tc0
    d.x = d.x/Tc0
    l = d.props['L']
    d.x = d.x * pow(float(l),1/nu0)
    d.y = d.y * pow(float(l),-alpha0/nu0)
plt.cla()
pyalps.plot.plot(M)
plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
plt.ylabel(r'$C_VL^{\frac{-\alpha}{\nu}}$', fontsize='x-large')
plt.title(r'$C_V$', fontsize='x-large')

axcolor = 'lightgoldenrodyellow'
axle = plt.axes([0.1, 0.10, 0.7, 0.03], axisbg=axcolor)
axlo = plt.axes([0.1, 0.05, 0.7, 0.03], axisbg=axcolor)
axtc = plt.axes([0.1, 0.15, 0.7, 0.03], axisbg=axcolor)

sle = Slider(axle, r'$\lambda^e$', le_min, le_max, valinit=le0)
slo = Slider(axlo,r'$\lambda^o$', lo_min,lo_max, valinit=lo0)
stc = Slider(axtc, r'$T_c$', Tc_min, Tc_max, valinit=Tc0)

def update(val):
    lo = slo.val
    le = sle.val
    Tc = stc.val
    alpha = f_alpha(le,lo)
    beta  = f_beta(le,lo) 
    gamma = f_gamma(le,lo)
    nu    = f_nu(le,lo) 
    chi = pyalps.collectXY(data,x='T',y='susceptibility staggered',foreach=['L'])
    plt.sca(axChi)
    for d in chi:
        d.props['label']='L='+str(d.props['L'])
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),1/nu)
        d.y = d.y * pow(float(l),-gamma/nu)
    plt.cla()
    pyalps.plot.plot(chi)
    plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
    plt.ylabel(r'$\chi L^{\frac{-\gamma}{\nu}}$', fontsize='x-large')
    plt.title(r'$\chi$', fontsize='x-large')
    cv = pyalps.collectXY(data,x='T',y='c_V',foreach=['L'])
    plt.sca(axCv)
    for d in cv:
        d.props['label']='L='+str(d.props['L'])
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),1/nu)
        d.y = d.y * pow(float(l),-alpha/nu)
    plt.cla()
    pyalps.plot.plot(cv)
    plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
    plt.ylabel(r'$C_V L^{\frac{-\alpha}{\nu}}$', fontsize='x-large')
    plt.title(r'$cv$', fontsize='x-large')
    binder = pyalps.collectXY(data,x='T',y='BinderCumulant staggered',foreach=['L'])
    plt.sca(axBinder)
    for d in binder:
        d.props['label']='L='+str(d.props['L'])
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),1/nu)
    plt.cla()
    pyalps.plot.plot(binder)
    plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
    plt.ylabel(r'$U^4$', fontsize='x-large')
    plt.title('Binder', fontsize='x-large')
    M = pyalps.collectXY(data,x='T',y='M staggered',foreach=['L'])
    plt.sca(axM)
    for d in M:
        d.props['label']='L='+str(d.props['L'])
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),1/nu)
        d.y = d.y * pow(float(l),beta/nu)
    plt.cla()
    pyalps.plot.plot(M)
    plt.xlabel(r'$(T-T_c)L^{\frac{1}{\nu}}$', fontsize='x-large')
    plt.ylabel(r'$ML^{\frac{\beta}{\nu}}$', fontsize='x-large')
    plt.title(r'$M$', fontsize='x-large')
    fig.canvas.draw_idle()

sle.on_changed(update)
stc.on_changed(update)
slo.on_changed(update)
plt.show()
