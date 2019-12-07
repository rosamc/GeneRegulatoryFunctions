import sys, os, re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import BoundaryFinder as BF
import json


def get_common_boundary(mats,matpars,col_ar=None,row_ar=None):
    """Gets the boundary common to a set of matrices. Arguments:  
    - mats: a set of matrices with 1 on the grid positions where a parameter set was found
    - matpars: the corresponding matrices with parameter values
    - col_ar: bins of the columns.
    - row_ar: bins of the rows. 
    It returns:
    - A matrix with 1 on the positions that correspond to the common boundary.
    - A pandas dataframe with the column, row and parameter values of these boundary points. 
    The same cell position can appear more than once. """
    if len(mats)>1:
        mat_sum=np.sum(mats,axis=0)
        mat_sum[mat_sum>0]=1
    else:
        mat_sum=mats[0]
    print(mat_sum.shape)
    B=BF.BoundaryExplorer(col_ar=col_ar,row_ar=row_ar,mat=mat_sum)
    B.find_boundary_points()
    rows=[]
    cols=[]
    parameters=[]
    for row in range(len(row_ar)):
        for col in range(len(col_ar)):
            if B.indices_boundary_all[row,col]>0:
                for m in range(len(mats)):
                    if mats[m][row,col]>0:
                        rows.append(row_ar[row])
                        cols.append(col_ar[col])
                        parameters.append(matpars[m][row,col])
                    
    return [B.indices_boundary_all,pd.DataFrame({'row':rows,'col':cols,'parameters':parameters})]

def read_settings(filename):
    """Old function to read file with settings and args dictionary. Should not be needed if using the json dumping system. """
    
    #else:
    #last_iter=last_iter_m
    #t_elapsed='NA'

    row_ar=[]
    col_ar=[]
    rows=False
    cols=False

    for l in open(filename,'r').readlines():
        l=l.strip()
        #print(l)
        if l.startswith('row_ar'):
            rows=True
            row_ar.extend(list(map(float,l.split('[')[1].split())))
        elif rows is True:
            if ']' in l:
                rows=False
                row_ar.extend(list(map(float,l.split(']')[0].split())))
            else:
                row_ar.extend(list(map(float,l.split())))
        if l.startswith('col_ar'):
            cols=True
            col_ar.extend(list(map(float,l.split('[')[1].split())))
        elif cols is True:
            if ']' in l:
                cols=False
                col_ar.extend(list(map(float,l.split(']')[0].split())))
            else:
                col_ar.extend(list(map(float,l.split())))
        if 'prob_par' in l:
            prob_par=l.strip(',').split(':')[1]
            
        if 'prob_replace' in l:
            prob_replace=l.strip(',').split(':')[1]
           
        if 'niters_conv' in l and not 'points' in l:
            niters_conv=l.strip(',').split(':')[1]
           
        if 'niters_conv' in l and 'points' in l:
            niters_conv_pt=l.strip(',').split(':')[1]
            
        if 'extr_uniform' in l:
            extr=l.strip(',').split(':')[1]
    return [np.array(row_ar), np.array(col_ar),prob_par,prob_replace,niters_conv,niters_conv_pt,extr]

def plot_boundaries_search(njobs=1,final=True, printtocheck=True, fldr='',basename='', 
                           joinmats=True,jid_num=None, reference=None, xlabel='position', ylabel='steepness',jsonf=True):
    """Plots the boundaries generated in a parallel search. 
    njobs: number of parallel jobs run.
    final: True/ False depending on whether the jobs finished (True) or were cut due to time limit on the cluster (False).
    printtocheck: if True, prints the boundary points to a file in order to check with mathematica.
    flder: directory where the results to analyse are.
    basename: name given to the matrices/settings file when executing the search. 
    joinmats: if True, will return the common boundary and the dataframe with the corresponding points (via get_common_boundary function)
    jid_num: jid of the parallel search in O2.
    reference: in case there is a reference boundary that wants to be overlayed, pass here as a 2D array where each row is a point of (col, row) 
    json: set to False only for backward compatibility, when settings and args dictionaries were not saved as json.
    """
    basename_mat='mat_%s'%basename
    basename_mat_pars='mat_pars_%s'%basename
    pat_mat=re.compile('%s_([0-9]*).npy'%basename_mat)
    if printtocheck:
        folder_tocheck=os.path.join(fldr,'tocheck') #here write output files to check with mathematica
        if not os.path.isdir(folder_tocheck):
            os.mkdir(folder_tocheck)
        print('folder to check',folder_tocheck)
    if joinmats:
        matslist=[]
        matsparslist=[]
    
    outf=os.path.join(fldr,'final_results')
    for i in range(njobs):
        if final is False:
            #In case the search did not finish due to maximum time allowed on cluster reached
            outfolder=os.path.join(fldr,basename+'_out_%d'%(i))
            mats=[f for f in os.listdir(outfolder) if basename_mat in f]
            iters=[int(pat_mat.findall(f)[0]) for f in mats]
            argsort=np.argsort(iters)
            last_iter_m=iters[argsort[-1]]
            last_iter=last_iter_m

            mat=np.load(os.path.join(outfolder,basename_mat+'_%d.npy'%(last_iter_m)))
            mat_pars=np.load(os.path.join(outfolder,basename_mat_pars+'_%d.npy'%(last_iter_m)))
            timediff='-'
            converged='-'
            cont=True
        else:       
            try:
                mat=np.load(os.path.join(outf,'%s_%d_last.npy'%(basename_mat,i)))
                mat_pars=np.load(os.path.join(outf,'%s_%d_last.npy'%(basename_mat_pars,i)))
                if jid_num is not None:
                    stdoutfh=open(os.path.join(fldr,'%s_%d.out'%(jid_num,i+1)),'r')
                    stdout=stdoutfh.readlines()[-2:]
                    stdoutfh.close()
                    timediff=stdout[0].split(':')[1].strip()
                    converged=stoud[1].strip()
                else:
                    print('no jid_num, timediff and converged unknown.')
                    timediff=''
                    converged=''
                cont=True
                
            except FileNotFoundError as e:
                print('%d Not found '%i)
                print(e)
                cont=False
        if i%4==0:
            if i>0:
                plt.tight_layout()
                plt.show()
            fig,axes=plt.subplots(1,4,figsize=(20,8))
        if cont:

            fnamesett=os.path.join(outf,'%s_%d.sett'%(basename,i))
            if jsonf:
                settings=json.load(open(fnamesett))
                fnameargs=os.path.join(outf,'%s_%d.args'%(basename,i))
                args=json.load(open(fnameargs))
                col_ar=list(map(float,settings['col_ar'].split(',')))
                row_ar=list(map(float,settings['row_ar'].split(',')))
                prob_par=args['prob_par']
                prob_replace=args['prob_replace']
                niters_conv=args['niters_conv']
                niters_conv_points=args['niters_conv_points']
                extr=args['extr_uniform']
            else:
                row_ar, col_ar,prob_par,prob_replace,niters_conv,niters_conv_points,extr=read_settings(fnamesett)


            #print(row_ar)
            

            ax=axes[i%4]
            ax.set_title('%d\n,p.par=%s,p.repl=%s\n niters_conv=%s,pt=%s,extr=%s\n %s\n converged=%s '%(i,prob_par,prob_replace,niters_conv,niters_conv_points,extr,timediff,converged))
            ax.imshow(mat,origin='lower',extent=[col_ar[0],col_ar[-1],row_ar[0],row_ar[-1]],cmap=plt.cm.Greys)
            if reference is not None:
                ax.scatter(reference[:,0],reference[:,1],color='r',s=4,alpha=0.5)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            
            if joinmats or printtocheck:
                B=BF.BoundaryExplorer(col_ar=col_ar,row_ar=row_ar,mat=mat)
                B.find_boundary_points()
                
            if joinmats:
                matslist.append(np.asarray(B.indices_boundary_all,dtype=int))
                mat_pars_c=mat_pars.copy()
                mat_pars_c[~B.indices_boundary_all]=0
                
                matsparslist.append(mat_pars_c)
            

            if printtocheck:
                #print x,y,parameters to check with mathematica
                outfile=open(os.path.join(folder_tocheck,'%s_%d.in'%(basename_mat,i)),'w')
                
                for row in range(len(row_ar)):
                    for col in range(len(col_ar)):
                        if B.indices_boundary_all[row,col]>0:
                            ax.scatter(col_ar[col],row_ar[row],color='g',s=4,alpha=0.5)
                            pars=list(map(str,mat_pars[row,col]))
                            outfile.write(str(col_ar[col])+','+str(row_ar[row])+','+','.join(pars)+'\n')
                outfile.close()
    if printtocheck:
        print('folder to check with mathematica is', folder_tocheck)
    if joinmats:
        return get_common_boundary(matslist,matsparslist,row_ar=row_ar,col_ar=col_ar)
    else:
        return


