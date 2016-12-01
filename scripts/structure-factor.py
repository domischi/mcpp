#************************************************************************# 
# This simple script allows for a graphical temperature Sweep of the     #
# structure factor modulus squared (proportional to a scattering signal) #
#************************************************************************#

# I know that the script pollutes the standard output, however up to now 
# I didn't find a solution, as the problem is with a dangeling pointer.
# This issue is resolved with IPython, however as alps does link against
# the standard python this is not a possibility. For now just ignore the 
# ouput of the save function

import pyalps
from numpy import linspace, sqrt
from os import path, mkdir, listdir
import numpy as np
import re

T=True
F=False

ANALYZE       =T
GENERATE_INPUT=F
RUN_SIMULATION=F

if ANALYZE:
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, CheckButtons
    import pyalps.plot
if GENERATE_INPUT:
    #prepare the input parameters
    parms = []
    l=32
    for t in linspace(0.0,1.2,25):
        parms.append(
            { 
                 'LATTICE'        : "square lattice", 
                 'Initialization' : 'Vortex',
                 'T'              : t,
                 'D'              : 1.,
                 'measure last configuration' : True,
                 'THERMALIZATION' : 50000,
                 'SWEEPS'         : 15000,
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 3.0,
                 'L'              : l,
                 'structure_factor': True,
                 'Targeted Acceptance Ratio': 0.4,
                 'Each_Measurement': 15
            }
           )
    #write the input file and run the simulation
    input_file = pyalps.writeInputFiles('parm',parms)
if RUN_SIMULATION:
    pyalps.runApplication('mc++',input_file,Tmin=5,MPI=1)
if ANALYZE:
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

    fig=plt.figure()
    left_start=0.05
    right_end=1-left_start
    bottom_start=0.3
    top_end=1
    boarder_width=0.1
    width=0.5-boarder_width
    height_plots=(1-bottom_start)/3.-boarder_width

# HM= heat map, TS= temperature slider
    axHM = plt.axes([left_start                    , bottom_start,                      width, .98-bottom_start])
    axCB = plt.axes([left_start+width+boarder_width*0.2, bottom_start,                      boarder_width*0.2, .98-bottom_start])
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

    def AdjustCheckButton(Bu):
        CheckButtonsRectangleY=0.1
        CheckButtonsRectangleH=0.8
        CheckButtonsCrossL=0.05
        CheckButtonsCrossR=0.17
        CheckButtonsCrossD=0.1
        CheckButtonsCrossU=0.9
        Bu.rectangles[0].set_y     (CheckButtonsRectangleY)
        Bu.rectangles[0].set_height(CheckButtonsRectangleH)
        Bu.lines[0][0].set_xdata([CheckButtonsCrossL,CheckButtonsCrossR])
        Bu.lines[0][0].set_ydata([CheckButtonsCrossU,CheckButtonsCrossD])
        Bu.lines[0][1].set_xdata([CheckButtonsCrossL,CheckButtonsCrossR])
        Bu.lines[0][1].set_ydata([CheckButtonsCrossD,CheckButtonsCrossU])

    left_In   = left_start 
    bottom_In = bottom_save
    height_In = height_save*0.4
    width_In  = 0.3
    axIn = plt.axes([left_In, bottom_In, width_In, height_In])
    BuIn = CheckButtons(axIn, ('Interpolate',), (False,))
    AdjustCheckButton(BuIn)
    Interpolate=False
    def helper_interpolate(val):
        global Interpolate
        Interpolate=not Interpolate
    BuIn.on_clicked(helper_interpolate)

    left_SLy  = left_start 
    bottom_SLy= bottom_save+1.2*height_In
    height_SLy= height_In
    width_SLy = 0.3
    axSLy= plt.axes([left_SLy, bottom_SLy, width_SLy, height_SLy])
    BuSLy= CheckButtons(axSLy, ('Log y',), (False,))
    AdjustCheckButton(BuSLy)
    SLyPlot=False
    def helper_slyplot(val):
        global SLyPlot
        SLyPlot=not SLyPlot
        update(0)
    BuSLy.on_clicked(helper_slyplot)

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

    def plot_heatmap(ax, S2_data, S2_err):
        plt.sca(ax)
        plt.cla()
        plt.xlabel(r'$k_x$')
        plt.ylabel(r'$k_y$')
	dx=1./len(S2_data)
	dy=1./len(S2_data)
	y, x = np.mgrid[slice(0, 1 + dy, dy),
			slice(0, 1 + dx, dx)]
        if not SLyPlot:
            im = ax.pcolormesh(x, y,S2_data)
        else:
            im = ax.pcolormesh(x, y,np.log( S2_data))
        fig.colorbar(im, cax=axCB)
    def plot_0k(ax,S2_data,S2_err):
        plt.sca(ax)
        plt.cla()
        plt.xlabel(r'$k$')
        plt.ylabel(r'$|S(0,k)|^2$')
        if not SLyPlot:
            plt.errorbar(linspace(0,1,L), Get_0k_data(S2_data,L),yerr=Get_0k_data(S2_err,L))
        else:
            plt.semilogy(linspace(0,1,L), Get_0k_data(S2_data,L))
    def plot_k0(ax,S2_data,S2_err):
        plt.sca(ax)
        plt.cla()
        plt.xlabel(r'$k$')
        plt.ylabel(r'$|S(k,0)|^2$')
        if not SLyPlot:
            plt.errorbar(linspace(0,1,L), Get_k0_data(S2_data,L),yerr=Get_k0_data(S2_err,L))
        else:
            plt.semilogy(linspace(0,1,L), Get_k0_data(S2_data,L))
    def plot_kk(ax,S2_data,S2_err):
        plt.sca(ax)
        plt.cla()
        plt.xlabel(r'$k$')
        plt.ylabel(r'$|S(k,k)|^2$')
        if not SLyPlot:
            plt.errorbar(linspace(0,1,L), Get_kk_data(S2_data,L),yerr=Get_kk_data(S2_err,L))
        else:
            plt.semilogy(linspace(0,1,L), Get_kk_data(S2_data,L))

    def plot_data(S2_data, S2_err):
        plot_heatmap(axHM,S2_data, S2_err)
        plot_0k(ax0k,S2_data, S2_err)
        plot_kk(axkk,S2_data, S2_err)
        plot_k0(axk0,S2_data, S2_err)
        fig.canvas.draw() 
       
    def get_S2_val():
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
        normalization_constraint=max(data.flatten())
        data=data/normalization_constraint
        error =error /normalization_constraint
        return data, error

    def update(dummy):
        data, error=get_S2_val()
        plot_data(data, error)
        
    NEW_SESSION=True
    SESSION_ID=0
    
    def GetIDString():
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
        return str(SESSION_ID)+'.'+str(FILE_ID)
    def GetSnapshotName(SessionID):
        PathToWorkingFolder='.'
        PathToResultsFolder=PathToWorkingFolder+'/snapshots'
        return PathToWorkingFolder+'/snapshots/snapshot.'+SessionID+'.pdf'
    def GetHeatmapName(SessionID):
        PathToWorkingFolder='.'
        PathToResultsFolder=PathToWorkingFolder+'/snapshots'
        return PathToWorkingFolder+'/snapshots/heatmap.'+SessionID+'.pdf'
    def Get0kName(SessionID):
        PathToWorkingFolder='.'
        PathToResultsFolder=PathToWorkingFolder+'/snapshots'
        return PathToWorkingFolder+'/snapshots/Plot.0k.'+SessionID+'.pdf'
    def GetkkName(SessionID):
        PathToWorkingFolder='.'
        PathToResultsFolder=PathToWorkingFolder+'/snapshots'
        return PathToWorkingFolder+'/snapshots/Plot.kk.'+SessionID+'.pdf'
    def Getk0Name(SessionID):
        PathToWorkingFolder='.'
        PathToResultsFolder=PathToWorkingFolder+'/snapshots'
        return PathToWorkingFolder+'/snapshots/Plot.k0.'+SessionID+'.pdf'

    def save_helper(SessionID, NameFunction, data, err, PlotFunction):
        f=plt.figure()
        ax=f.gca()
        PlotFunction(ax,data,err)
        f.savefig(NameFunction(SessionID))
        try:
            plt.close(f)
        except:
            pass
    def save_heatmp(SessionID, NameFunction, data, err, PlotFunction):
        f=plt.figure()
        ax=f.gca()
        plt.xlabel(r'$k_x$')
        plt.ylabel(r'$k_y$')
	dx=1./len(data)
	dy=1./len(data)
	y, x = np.mgrid[slice(0, 1 + dy, dy),
			slice(0, 1 + dx, dx)]
        if not SLyPlot:
            im = ax.pcolormesh(x, y,data)
        else:
            im = ax.pcolormesh(x, y,np.log(data))
        fig.colorbar(im, ax=ax)
        f.savefig(NameFunction(SessionID))
        try:
            plt.close(f)
        except:
            pass

    def save(dummy):
        SessionID=GetIDString()
        plt.savefig(GetSnapshotName(SessionID))
        data,err=get_S2_val()
        save_heatmp(SessionID, GetHeatmapName, data, err, plot_heatmap)
        save_helper(SessionID, Getk0Name,      data, err, plot_k0)
        save_helper(SessionID, Get0kName,      data, err, plot_0k)
        save_helper(SessionID, GetkkName,      data, err, plot_kk)
    SlT.on_changed(update)
    BuSv.on_clicked(save)
    update(0)

    plt.show()
