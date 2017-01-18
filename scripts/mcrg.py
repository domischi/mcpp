# ****************************************************************************
# 
# This does generate input / analyze data for the mc++ code with respect to its MCRG component
# 
# ****************************************************************************

GENERATE_INPUT=False
ANALYZE       =True
mcrg_iteration_depth=3
if GENERATE_INPUT:
    import pyalps
    import pyalps.alea as alpsalea
    import numpy as np
    from pyalps.math import exp
    from numpy import linspace, sqrt, dot, zeros, log, sort, array, outer, savetxt, transpose,mean, std, delete, vectorize
    from numpy.linalg import eigvals, inv, cond, norm, eig
    from numpy.random import normal
    from functools import reduce
    T=linspace(2.2,2.3,5)
    n_clones=1
    #prepare the input parameters
    parms = []
    l=135
    for interactions in ['dXY']:
        for reduction in ['FerroBlockspin']:
            for t in T:
                parms.append(
                    { 
                         'LATTICE'        : "square lattice", 
                         'Initialization'  : "Random", 
                         'Ising'          : True,
                         'NUM_CLONES'     : n_clones,
                         'T'              : t,
                         'D'              : 1.,
                         'THERMALIZATION' : 5000,
                         'SWEEPS'         : 3000,
                         'UPDATE'         : "ssf",
                         'L'              : l,
                         'mcrg_iteration_depth': mcrg_iteration_depth,
                         'MCRG Interactions': interactions,
                         'MCRG Reduction Technique': reduction,
                         'Each_Measurement': 25 
                       }
                   )

    #write the input file and run the simulation
    input_file = pyalps.writeInputFiles('parm',parms)
if ANALYZE:
    import pyalps
    import matplotlib.pyplot as plt
    import pyalps.plot
    import pyalps.alea as alpsalea
    from pyalps.math import exp
    import numpy as np
    from numpy import linspace, sqrt, dot, zeros, log, sort, array, outer, savetxt, transpose,mean, std, delete, vectorize
    from numpy.linalg import eigvals, inv, cond, norm, eig
    from numpy.random import normal

    mcrg_iteration_depth=2

#If a value is just slightly imaginary, use this function
    def RawSanitize(val, tol=1e-3):
        if(val==0 or abs(val.imag/val.real)<tol):
            return val.real
        else:
            print('Attention required: '+str(val))
            return np.nan
            #return val
    Sanitize=vectorize(RawSanitize)

#Load the mean values of the arrays and the Jackknife data to do the error analysis
    def LoadJackknife(filename,observable_name,respath='/simulation/results/'):
        ar=pyalps.h5.archive(filename,'r')
        return ar[respath+observable_name+'/jacknife']['data']
    def LoadMean(filename,observable_name,respath='/simulation/results/'):
        ar=pyalps.h5.archive(filename,'r')
        return ar[respath+observable_name+'/mean/value']

#Calculates the critical exponents out of the arrays <S_a^n S_b^n>, <S_a^n-1 S_b^n>, <S_a^n>, <S_a^-1n>
    def EV(sasb_n, sasb_nm1, sa_n, sa_nm1,debug1=False,debug2=False):
        try:
            dim = len(sa_n)
            sasb_n  =  sasb_n.reshape(dim , dim)
            sasb_nm1=sasb_nm1.reshape(dim , dim)
            dsdk_n  =sasb_n  -outer(sa_n,sa_n)
            dsdk_nm1=sasb_nm1-outer(sa_nm1,sa_n)
            T=dot(dsdk_nm1,inv(dsdk_n))
            if debug2:
                print(cond(T))
            #evs=sort(abs(eigvals(T)))
            evs=sort(eigvals(T))
            if debug1:
                val,vec=eig(T)
                print(val)
                print('-----------------')
                for i in range(len(val)):
                    print(val[i],vec[:,i])
                print('-----------------')
            lambda_ =max(abs(evs))
            lambda_test =evs[-1]
            #assert(lambda_==lambda_test)
            if(len(evs)>1):
                lambda_2=evs[-2]
            else:
                lambda_2=0
            lambda_=RawSanitize(lambda_)
            return lambda_, lambda_2
        except np.linalg.linalg.LinAlgError:
            return 0,0

#Use the Jackknife data, do the Eigenvalue analysis over them, and estimate the error out of this
    def JackknifeEVPrime(j_sasb_n, j_sasb_nm1, j_sa_n, j_sa_nm1, lambda_):
        s=0.
        jackknife_length=len(j_sasb_n)
        assert(len(j_sasb_n)  ==jackknife_length) 
        assert(len(j_sasb_nm1)==jackknife_length)
        assert(len(j_sa_n)    ==jackknife_length)
        assert(len(j_sa_nm1)  ==jackknife_length)
        for i in range(jackknife_length):#leave away measurement at index i
            sasb_n    = mean( delete(j_sasb_n   ,i , axis=0), axis=0)
            sasb_nm1  = mean( delete(j_sasb_nm1 ,i , axis=0), axis=0)
            sa_n      = mean( delete(j_sa_n     ,i , axis=0), axis=0)
            sa_nm1    = mean( delete(j_sa_nm1   ,i , axis=0), axis=0)
            lambda_p=EV(sasb_n, sasb_nm1, sa_n, sa_nm1)
            
            s+=(lambda_-lambda_p)**2
        return sqrt((jackknife_length-1.)/(jackknife_length)*s)

    dictionary_reduction_types={}
    num_max_reduction_types=4
    def num_reduction_types():
        return len(dictionary_reduction_types)
    def add_to_reduction_dict(s):
        if s in dictionary_reduction_types:
            return
        else:
            dictionary_reduction_types[s]=len(dictionary_reduction_types)
    def index_reduction_type(s):
        return dictionary_reduction_types[s]
    def reduction_type_by_index(i):
        return [key for key, value in dictionary_reduction_types.items() if value == i][0]

    dictionary_interactions={}
    num_max_interactions=4
    def num_interactions():
        return len(dictionary_interactions)
    def add_to_interaction_dict(s):
        if s in dictionary_interactions:
            return
        else:
            dictionary_interactions[s]=len(dictionary_interactions)
    def index_interaction_set(s):
        return dictionary_interactions[s]
    def interaction_set_by_index(i):
        return [key for key, value in dictionary_interactions.items() if value == i][0]

    filenames=pyalps.getResultFiles(prefix='parm')
    num_T=len(filenames) #separate files for blockspin and decimation, however same file for different iterations and odd even

    LVsT          =np.zeros([2,mcrg_iteration_depth,num_max_reduction_types, num_max_interactions, num_T]) # (odd|even), iteration, reduction technique, interaction set
    PlotT         =np.zeros([2,mcrg_iteration_depth,num_max_reduction_types, num_max_interactions, num_T]) # (odd|even), iteration, reduction technique, interaction set
    counter_matrix=np.zeros([2,mcrg_iteration_depth,num_max_reduction_types, num_max_interactions])        # (odd|even), iteration, reduction technique, interaction set
    for f in filenames:
        filename=f
        if filename[-4:]=='.xml':
            filename=filename[:-4]+'.h5'
        T=pyalps.loadProperties([filename])[0]['T']
        reduction_type =pyalps.loadProperties([filename])[0]['MCRG Reduction Technique']
        interaction_set=pyalps.loadProperties([filename])[0]['MCRG Interactions']
        add_to_reduction_dict(reduction_type) 
        add_to_interaction_dict(interaction_set) 
        print(index_reduction_type(reduction_type),reduction_type, index_interaction_set(interaction_set), interaction_set ,f)
        for type_of_interaction in ['e', 'o']:
            for it in range(1,mcrg_iteration_depth+1): 
                data = (pyalps.loadMeasurements([filename], ['MCRG S_alpha'+str(it-1),'MCRG S_alpha'+str(it),'MCRG S_alpha'+str(it-1)+' S_beta'+str(it),'MCRG S_alpha'+str(it)+' S_beta'+str(it)]))

                mean_sasb_n=  LoadMean(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it)  +' S_beta'+str(it))
                mean_sasb_nm1=LoadMean(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it-1)+' S_beta'+str(it))
                mean_sa_n=    LoadMean(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it)  )
                mean_sa_nm1=  LoadMean(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it-1))
               
                #jknf_sasb_n=  LoadJackknife(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it)  +' S_beta'+str(it))
                #jknf_sasb_nm1=LoadJackknife(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it-1)+' S_beta'+str(it))
                #jknf_sa_n=    LoadJackknife(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it)  )
                #jknf_sa_nm1=  LoadJackknife(filename,'MCRG' + type_of_interaction+ ' S_alpha'+str(it-1))
                #print(reduction_type) 
                #print(EV(mean_sasb_n, mean_sasb_nm1, mean_sa_n, mean_sa_nm1, debug1=False, debug2=False))
                lambda_ = EV(mean_sasb_n, mean_sasb_nm1, mean_sa_n, mean_sa_nm1, debug1=False, debug2=False)[0]
                lambda_ = RawSanitize(lambda_)
                counter=int(counter_matrix[int('e'==type_of_interaction),it-1, index_reduction_type(reduction_type), index_interaction_set(interaction_set)])
                LVsT [int('e'==type_of_interaction),it-1, index_reduction_type(reduction_type), index_interaction_set(interaction_set), counter] = lambda_
                PlotT[int('e'==type_of_interaction),it-1, index_reduction_type(reduction_type), index_interaction_set(interaction_set), counter] = T
                counter_matrix[int('e'==type_of_interaction),it-1, index_reduction_type(reduction_type), index_interaction_set(interaction_set)]+= 1



    for actual_index in range(0,mcrg_iteration_depth):
        fig, (ax1,ax2) = plt.subplots(1,2)
        plt.title('Iteration '+str(actual_index+1))
        for reduction_type in range(num_reduction_types()):
            for interaction_set in range(num_interactions()):
                plt.sca(ax1)
                plt.plot(PlotT[0,actual_index,reduction_type,interaction_set], LVsT[0,actual_index,reduction_type,interaction_set],'*', label=reduction_type_by_index(reduction_type)+' '+interaction_set_by_index(interaction_set))
                plt.sca(ax2)
                plt.plot(PlotT[1,actual_index,reduction_type,interaction_set], LVsT[1,actual_index,reduction_type,interaction_set],'*', label=reduction_type_by_index(reduction_type)+' '+interaction_set_by_index(interaction_set))
        plt.sca(ax1)
        plt.plot([0,3],[7.84]*2,'b-')
        plt.plot([2.2]*2,[0,10],'b-')
        plt.legend(loc='best')
        plt.xlabel(r'$T$')
        plt.ylabel(r'$\lambda^o$')
        plt.sca(ax2)
        plt.plot([0,3],[3]*2,'b-')
        plt.plot([2.2]*2,[0,10],'b-')
        plt.legend(loc='best')
        plt.xlabel(r'$T$')
        plt.ylabel(r'$\lambda^e$')
    plt.show()
