# ****************************************************************************
# 
# This tutorial shows how to simply generate the input files, start the task in python (not recommended on a cluster) and then evaluate. It is designed to the dipolar XY system in 2D, where the order parameter is a staggered magnetic field (stripe like)
#
# To a big proportion this tutorial is written by Brigitte Surer <surerb@phys.ethz.ch> in the spinmc part of the ALPS library. As it is similar the code was reused. 
# 
# ****************************************************************************

import os
import pyalps
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import pyalps.plot
from numpy import linspace

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),['M staggered', 'c_V', 'BinderCumulant staggered', 'susceptibility staggered'])
#M = pyalps.collectXY(data,x='T',y='M staggered',foreach=['L'])
#chi = pyalps.collectXY(data,x='T',y='susceptibility staggered',foreach=['L'])
#c_v= pyalps.collectXY(data,x='T',y='c_V',foreach=['L'])
#binder = pyalps.collectXY(data,x='T',y='BinderCumulant staggered',foreach=['L'])

Tc0=1.2
Tc_min=1.
Tc_max=1.5
nu0=1.6
nu_min=0.1
nu_max=3
gamma0= 7./4
gamma_min=0.
gamma_max=3
beta0=0.18
beta_min=0.1
beta_max=0.3


fig, axM = plt.subplots()
plt.subplots_adjust(left=0.375, bottom=0.3, right = 0.625)
axChi = plt.axes([0.05, 0.3, 0.2 ,0.6])
axBinder = plt.axes([0.7, 0.3, 0.2, 0.6])
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
    #d.y = d.y * pow(float(l),-gamma/nu)
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

axcolor = 'lightgoldenrodyellow'
axbeta = plt.axes([0.1, 0.05, 0.7, 0.03], axisbg=axcolor)
axgamma = plt.axes([0.1, 0.1, 0.7, 0.03], axisbg=axcolor)
axnu = plt.axes([0.1, 0.15, 0.7, 0.03], axisbg=axcolor)
axtc = plt.axes([0.1, 0.2, 0.7, 0.03], axisbg=axcolor)

sgamma = Slider(axgamma, r'$\gamma$', gamma_min, gamma_max, valinit=gamma0)
snu = Slider(axnu, r'$\nu$', nu_min, nu_max, valinit=nu0)
stc = Slider(axtc, r'$T_c$', Tc_min, Tc_max, valinit=Tc0)
sbeta = Slider(axbeta,r'$\beta$', beta_min,beta_max, valinit=beta0)

def update(val):
    gamma = sgamma.val
    nu = snu.val
    beta = sbeta.val
    Tc = stc.val
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
    plt.ylabel(r'$\chi^{\frac{-\gamma}{\nu}}$', fontsize='x-large')
    plt.title(r'$\chi$', fontsize='x-large')
    binder = pyalps.collectXY(data,x='T',y='BinderCumulant staggered',foreach=['L'])
    plt.sca(axBinder)
    for d in binder:
        d.props['label']='L='+str(d.props['L'])
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),1/nu)
        #d.y = d.y * pow(float(l),-gamma/nu)
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
    plt.ylabel(r'$ML^{\frac{-\gamma}{\nu}}$', fontsize='x-large')
    plt.title(r'$M$', fontsize='x-large')
    fig.canvas.draw_idle()

sgamma.on_changed(update)
snu.on_changed(update)
stc.on_changed(update)
sbeta.on_changed(update)
plt.show()
