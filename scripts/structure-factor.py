# ****************************************************************************
# 
# This tutorial shows how to simply generate the input files, start the task in python (not recommended on a cluster) and then evaluate. It is designed to the dipolar XY system in 2D, where the order parameter is a staggered magnetic field (stripe like)
#
# To a big proportion this tutorial is written by Brigitte Surer <surerb@phys.ethz.ch> in the spinmc part of the ALPS library. As it is similar the code was reused. 
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
import pyalps.plot
from numpy import linspace, sqrt
from os import path, mkdir, listdir
import numpy as np
import re
#prepare the input parameters
parms = []
#for l in [16,32,64,128]: 
#for l in [32]: 
#    for t in [0,0.5]:
#    #for t in linspace(0.5,0.6,2):
#        parms.append(
#            { 
#                 'LATTICE'        : "square lattice", 
#                 'T'              : t,
#                 'D'              : 1.,
#                 'THERMALIZATION' : 110,
#                 'SWEEPS'         : 210,
#                 'UPDATE'         : "ssf",
#                 'cutoff_distance': 1.8,
#                 'L'              : l,
#                 'structure_factor': True,
#                 'Targeted Acceptance Ratio': 0.4,
#                 'Each_Measurement': 15
#            }
#           )
#
##write the input file and run the simulation
#input_file = pyalps.writeInputFiles('parm',parms)
##pyalps.runApplication('mc++',input_file,Tmin=5)
## use the following instead if you have MPI
#pyalps.runApplication('mc++',input_file,Tmin=5,MPI=1)
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),['|Structure Factor|^2'])

def Get_0k_range(L):
	return range(0,L*L,L)
def Get_k0_range(L):
	return range(0,L,1)
def Get_kk_range(L):
	return range(0,L*L,L+1)
def Get_0k_data(data,L):
	return data.flatten()[Get_0k_range(L)]
def Get_k0_data(data,L):
	return data.flatten()[Get_k0_range(L)]
def Get_kk_data(data,L):
	return data.flatten()[Get_kk_range(L)]

data = pyalps.collectXY(data,x='T',y='|Structure Factor|^2')[0]
T_arr=np.unique(data.x)
N_T=len(T_arr)
S2 = list()
S2E= list()
for i in range(0,len(data.y)):
    S2.append(data.y[i].mean)
    S2E.append(data.y[i].error)

L=int(sqrt(len(S2)/N_T))
S2 =np.array(S2 ).reshape((N_T,L,L))
S2E=np.array(S2E).reshape((N_T,L,L))

plt.figure()
left_start=0.05
right_end=1-left_start
bottom_start=0.3
top_end=1
boarder_width=0.1
width=0.5-boarder_width
height_plots=(1-bottom_start)/3.-boarder_width

# HM= heat map, TS= temperature slider
axHM = plt.axes([left_start                    , bottom_start,                      width, 1-bottom_start])
ax0k = plt.axes([left_start+width+boarder_width, bottom_start+(height_plots+boarder_width)*2, width, height_plots])
axk0 = plt.axes([left_start+width+boarder_width, bottom_start,                      width, height_plots])
axkk = plt.axes([left_start+width+boarder_width, bottom_start+height_plots+boarder_width, width, height_plots])

#GUI
height_slider=0.05
bottom_slider= bottom_start - 0.1
width_slider=right_end-left_start-0.05
axTS = plt.axes([left_start, bottom_slider, width_slider, height_slider])
SlT = Slider(axTS, 'T', T_arr[0], T_arr[-1], valinit=T_arr[0])

left_save  = 0.75
bottom_save= 0.05 
height_save= 0.1
width_save = 0.2
axSv = plt.axes([left_save, bottom_save, width_save, height_save])
BuSv = Button(axSv,'Save to File')

left_In   = left_start 
bottom_In = bottom_save
height_In = height_save
width_In  = 0.3
axIn = plt.axes([left_In, bottom_In, width_In, height_In])
BuIn = CheckButtons(axIn, ('Interpolate',), (False,))
Interpolate=False
def helper_interpolate(val):
	global Interpolate
	Interpolate=not Interpolate
BuIn.on_clicked(helper_interpolate)

def GetClosestTIndex(T):
	return np.abs(T_arr-T).argmin()
def GetClosestLowerIndex(T):
	idx=GetClosestTIndex(T)
	if T<T_arr[idx]:
		return idx-1
	return idx
def GetMixingRatio(i,T):
	if(i!=len(T_arr)-1):
		T_lower=T_arr[i]
		T_upper=T_arr[i+1]
		return 1-(T-T_lower)/(T_upper-T_lower)
	else:
		return 1

def plot_data(S2_data, S2_err):
	normalization_constraint=max(S2_data.flatten())
	S2_data=S2_data/normalization_constraint
	S2_err =S2_err /normalization_constraint
	plt.sca(axHM)
	plt.cla()
	plt.xlabel(r'$k_x$')
	plt.ylabel(r'$k_y$')
	plt.imshow(S2_data, interpolation='nearest', origin='lower')
	plt.sca(ax0k)
	plt.cla()
	plt.xlabel(r'$k$')
	plt.ylabel(r'$|S(0,k)|^2$')
	plt.errorbar(linspace(0,1,L), Get_0k_data(S2_data,L),yerr=Get_0k_data(S2_err,L))
	plt.sca(axk0)
	plt.cla()
	plt.xlabel(r'$k$')
	plt.ylabel(r'$|S(k,0)|^2$')
	plt.errorbar(linspace(0,1,L), Get_k0_data(S2_data,L),yerr=Get_k0_data(S2_err,L))
	plt.sca(axkk)
	plt.cla()
	plt.xlabel(r'$k$')
	plt.ylabel(r'$|S(k,k)|^2$')
	plt.errorbar(linspace(0,1,L), Get_kk_data(S2_data,L),yerr=Get_kk_data(S2_err,L))
	

def update(dummy):
	T=SlT.val
	if not Interpolate:
		i=GetClosestTIndex(T)
		if(T!=T_arr[i]):
			SlT.set_val(T_arr[i])
		data =S2[i,:,:]
		error=S2E[i,:,:]
	else:
		i=GetClosestLowerIndex(T)
		p=GetMixingRatio(i,T)
		if i+1!=len(T_arr):
			data =p*S2[i,:,:] +(1-p)*S2[i+1,:,:]
			error=p*S2E[i,:,:]+(1-p)*S2E[i+1,:,:]
		else:
			data =S2[i,:,:] 
			error=S2E[i,:,:]
		
	plot_data(data, error)
	
NEW_SESSION=True
SESSION_ID=0
def GetSaveName():
	global NEW_SESSION
	global SESSION_ID
	PathToWorkingFolder='.'
	PathToResultsFolder=PathToWorkingFolder+'/snapshots'
	prefix='snapshot'
	tmp_session_id=0
	tmp_file_id=0
	FILE_ID   =1
	if not path.exists(PathToResultsFolder):
		mkdir(PathToResultsFolder)
	for f in listdir(PathToResultsFolder):
		tmp=re.search('snapshot\.(\d+).(\d+).pdf', f)
		if tmp:
			tmp_session_id=max(int(tmp.group(1)), tmp_session_id) 
			tmp_file_id   =max(int(tmp.group(2)), tmp_file_id)
	if NEW_SESSION:
		NEW_SESSION=False
		SESSION_ID=tmp_session_id+1
	else: 		
		FILE_ID=tmp_file_id+1
	return PathToWorkingFolder+'/snapshots/snapshot.'+str(SESSION_ID)+'.'+str(FILE_ID)+'.pdf'

def save(dummy):
	plt.savefig(GetSaveName())
	
SlT.on_changed(update)
BuSv.on_clicked(save)
update(0)

plt.show()
