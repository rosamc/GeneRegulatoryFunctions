#Example of script to perform a boundary search on the cluster. This should be submitted by a bash script such as the run_jobarray.py
#Increase the number of iterations and niters_save (commented values are good to begin with).


import sys,os
import numpy as np
sys.path.append('../../utilsGRF/')
import BoundaryFinder as BF
import time
import itertools
import json

sys.path.append('../TFbindingP')
import N5P
def func(pars):
    #Boundary finder is expecting either the correct result or None, so this function allows this. At this moment I don´t know how to return None with pybind.
    out=[None,None]
    result=N5P.interfaceps_a_GRF_N5P_x(pars,np.array([1]))
    #print(result)
    if result[0]>0:
        out=[result[0],result[1]]
    return out






step=0.005
stp_ar=np.arange(0.2,1.6+step,step)
pos_ar=np.arange(0.3,1.3+step,step)


minp=1e-3
maxp=1e3

wpidx0=31
wpidx1=61
npars=62



constraintdict={x:{'min':1,'max':maxp} for x in range(wpidx0,wpidx1+1)}
constraintdict[0]={'fixed':1,'min':1,'max':1}
settings={'pars_limit':[minp,maxp],
          'compute_x_y_f':func,
          'wpindices':[wpidx0,wpidx1],
          'nsites':5,
          'npars':npars,
          'constraints':constraintdict,
           'row_ar':stp_ar,
          'col_ar':pos_ar,
          'seed':1,
         'mat':None,
         'mat_pars':None}

 
jid=int(sys.argv[1])-1  #this is an integer from 1 to as many combinations of metaparameters are tried (16 in this case) so that they can be run in parallel.

extremesu=[[-2,2],[-1.5,1.5],[-1,1],[-0.5,0.5]]
prob_par=[0.2,0.5]
prob_replace=[0.2,0.6]

combination=list(itertools.product(extremesu,prob_par,prob_replace))[jid]


extremesu,prob_par,prob_replace=combination

niters_conv=500
extr_uniform=[extremesu[0],extremesu[1]]

niters=10 #5000
L=15

name_save='poln5'
outfolder=name_save+'_out_%d'%jid

if not os.path.isdir(outfolder):
    os.mkdir(outfolder)


args={'niters':niters,
      'niters_conv':niters_conv,
      'niters_conv_points':1000,
      'niters_save':2, #100,
      'folder_save':outfolder,
       'name_save':name_save, 
      'prob_par':prob_par,
      'prob_replace':prob_replace,
      'extr_uniform':extr_uniform,
      'L_project':L,
      'plotting':False,
      'verbose':True,
       'dofirstmutate':True,'dopullcentroids':True,'dopulltangents':True}


#save parameters used for running the algorithm
def function_tostring(x):
    if isinstance(x, np.ndarray):
        return ','.join(map(str,x))
    else:
        return x.__name__
outfolder_final='final_results'
outfnames=[os.path.join(outfolder,name_save+'_%d.sett'%jid),os.path.join(outfolder_final,name_save+'_%d.sett'%jid)]
for fname in outfnames:
    outf=open(fname,'w')
    #outf.write(time.ctime()+'\n')
    #with open(outf, 'w') as file:
    json.dump(dict({'time':time.ctime()},**settings),outf,default=function_tostring) # use `json.loads` to do the reverse
    outf.close()

    outf=open(fname.replace('.sett','.args'),'w')
    #outf.write(time.ctime()+'\n')
    #with open(outf, 'w') as file:
    json.dump(dict({'time':time.ctime()},**args),outf) # use `json.loads` to do the reverse
    outf.close()

    #json.dump(args, outf)#    file.write(json.dumps(args))
    
#run
ti=time.time() 
BE=BF.BoundaryExplorer(**settings)
if settings['mat'] is None:
    BE.get_initial_points(10)
print('First 10 points.')
sys.stdout.flush()
BE.extend_boundary(**args)
name='%s_%d_last'%(name_save,jid)
np.save(os.path.join(outfolder_final,'mat_'+name+'.npy'),BE.mat)
np.save(os.path.join(outfolder_final,'mat_pars_'+name+'.npy'),BE.mat_pars)
te=time.time()
sys.stdout.flush()
print('time difference',te-ti)
print(BE.converged)
