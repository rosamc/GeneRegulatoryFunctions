"""
    Auxiliary code to plot boundaries.
    Copyright (C) <2021>  <Rosa Martinez Corral>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    """
import sys, os, re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import BoundaryFinder as BF
import json
import glob


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

def plot_boundaries_search(njobs=None,final=True, printtocheck=True, fldr='',basename='', 
                           joinmats=True,jid_num=None, reference=None, xlabel='position', ylabel='steepness',
                           jsonf=True,septime=":",getallpoints=False,unfinishedfolder=None,difparslimit=False):
    """Plots the boundaries generated in a parallel search. 
    njobs: number of parallel jobs run.
    final: True/ False depending on whether the jobs finished (True) or were cut due to time limit on the cluster (False).
    printtocheck: if True, prints the boundary points to a file in order to check with mathematica.
    flder: directory where the results to analyse are.
    basename: name given to the matrices/settings file when executing the search. 
    joinmats: if True, will return the common boundary and the dataframe with the corresponding points (via get_common_boundary function)
    jid_num: (string) jid of the parallel search in O2.
    reference: in case there is a reference boundary that wants to be overlayed, pass here as a 2D array where each row is a point of (col, row) 
    jsonf: set to False only for backward compatibility, when settings and args dictionaries were not saved as json.
    getallpoints: get all the points, not just the boundary, as a dataframe.
    unfinishedfolder: to be used with final=True in case not everything converged. Then this is the folderpath to a folder with the output folders were the intermediate results are saved.
    difparslimit: set to True if not all jobs were run for the same parameter limits, in this case it will group boundaries by parslimit
    """
    if final is True and unfinishedfolder is not None:
        print("final=True is incompatible with unfinishedfolder=True")
        raise(ValueError)
    basename_mat='mat_%s'%basename
    basename_mat_pars='mat_pars_%s'%basename
    pat_mat=re.compile('%s_([0-9]*).npy'%basename_mat)
    if printtocheck:
        folder_tocheck=os.path.join(fldr,'tocheck_%s'%basename) #here write output files to check with mathematica
        if not os.path.isdir(folder_tocheck):
            os.mkdir(folder_tocheck)
        print('folder to check',folder_tocheck)
    
    
    
    if final is False:
        folders=glob.glob(fldr+"/%s*"%basename) 
    else:
        outf=os.path.join(fldr,'final_results')
        folders=glob.glob(outf+"/mat_%s*.npy"%basename)
        folders=[x.replace("_last.npy",'') for x in folders]

    jids=np.array([int(x.split("_")[-1]) for x in folders])
    argsort_jids=np.argsort(jids)
    jids=jids[argsort_jids]
    folders=np.array(folders)[argsort_jids]
    #print(jids)
    #print(folders)
    if njobs is not None:
        if njobs!=len(jids): #sanity check
            print("njobs does not coincide with the number of jobs found. Exiting...")
            raise ValueError


    toprocess_unf=[]
    if unfinishedfolder is not None:
        #check if there is something here:
        folders_unfinished=glob.glob(unfinishedfolder+"/%s*"%basename)
        if len(folders_unfinished)>0:
            toprocess_unf=[os.path.join(unfinishedfolder,f) for f in folders_unfinished]

    
    #make lists of jids that correspond to same parameter limits
    group_jids=[]
    if difparslimit:
        #if final is False:
        #    raise("difparslimit is only prepared to deal with completed jobs")
        if jsonf is False:
            print("difparslimit is only prepared to deal with json settings")
            raise(ValueError)

        else:
            corresponding_parslimit=[]
            if final is True:
                
                for jid in jids:
                    fnamesett=os.path.join(outf,'%s_%d.sett'%(basename,jid))  
                    settings=json.load(open(fnamesett))
                    parslimit=",".join(list(map(str,settings["pars_limit"])))
                    corresponding_parslimit.append(parslimit)
            else:
                for folder in folders:
                    fnamesett=[x for x in os.listdir(folder) if x.endswith(".sett")][0]
                    settings=json.load(open(os.path.join(folder,fnamesett)))
                    parslimit=",".join(list(map(str,settings["pars_limit"])))
                    corresponding_parslimit.append(parslimit)
        
            unique_parslimit=np.unique(corresponding_parslimit)
            for parslimit in unique_parslimit:
                idxs=[x for x in range(len(corresponding_parslimit)) if corresponding_parslimit[x]==parslimit]
                #print(parslimit,idxs)
                jobs_p=[jids[i] for i in idxs]
                group_jids.append(jobs_p)

    else:
        group_jids.append(jids)

        #IDENTIFY THE ONES THAT SHARE SAME pars_limit and get the correspoding jids
        


    return_list=[]
    for jids2 in group_jids:
        print("Processing", jids2)
        if joinmats:
            matslist=[]
            matsparslist=[]
        if getallpoints:
            allpointslist=[]
            allpointscolnames=['row','col','parameters']

        if unfinishedfolder is not None:
            njobs=len(jids2)+len(toprocess_unf)
            itemslist=range(njobs)
        else:
            itemslist=range(len(jids2))

        for i_ in itemslist:
            if i_<len(jids2):
                i=jids2[i_]
                unfinished=False
            else:
                unfinished=True

            if final is False or unfinished:
                #In case the search did not finish due to maximum time allowed on cluster reached
                if final is False:
                    folder=folders[np.where(jids==i)[0][0]] #get folder corresponding to that jid
                    outfolder=os.path.join(fldr,folder)
                    outf=outfolder
                else:

                    outfolder=toprocess_unf[i_%len(jids)]
                    i=int(outfolder.split("_")[-1])

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
                        timediff=stdout[0].split(septime)[1].strip()
                        converged=stdout[1].strip()
                        if not converged in ["True","False"]: #if it was killed due to time limit, then it is only the iteration number. Discard.
                        	converged="-"
                        	timediff="-"
                    else:
                        print('no jid_num, timediff and converged unknown.')
                        timediff=''
                        converged=''
                    cont=True
                    
                except FileNotFoundError as e:
                    print('%d Not found '%i)
                    print(e)
                    cont=False
            if i_%4==0:
                if i_>0:
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
                

                ax=axes[i_%4]
                title='%d\n,p.par=%s,p.repl=%s\n niters_conv=%s,pt=%s,extr=%s\n %s, converged=%s '%(i,prob_par,prob_replace,niters_conv,niters_conv_points,extr,timediff,converged)
                #print(title)
                ax.set_title(title)
                ax.imshow(mat,origin='lower',extent=[col_ar[0],col_ar[-1],row_ar[0],row_ar[-1]],cmap=plt.cm.Greys)
                if reference is not None:
                    ax.scatter(reference[:,0],reference[:,1],color='r',s=4,alpha=0.5)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)

                if getallpoints:
                    for row in range(len(row_ar)):
                        for col in range(len(col_ar)):
                            if mat[row,col]>0:
                                allpointslist.append([row_ar[row],col_ar[col],mat_pars[row,col]])


                
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
                    checkfile='%s_%d.in'%(basename_mat,i)
                    outfile=open(os.path.join(folder_tocheck,checkfile),'w')
                    print("writing file to check in %s"%(checkfile))
                    
                    for row in range(len(row_ar)):
                        for col in range(len(col_ar)):
                            if B.indices_boundary_all[row,col]>0:
                                ax.scatter(col_ar[col],row_ar[row],color='g',s=4,alpha=0.5)
                                pars=list(map(str,mat_pars[row,col]))
                                outfile.write(str(col_ar[col])+','+str(row_ar[row])+','+','.join(pars)+'\n')
                    outfile.close()
        plt.tight_layout()
        plt.show()
        if printtocheck:
            print('folder to check with mathematica is', folder_tocheck)

        if getallpoints:
            allpdf=pd.DataFrame(np.array(allpointslist),columns=allpointscolnames)

        returnar=[]
        if joinmats:
            returnar.append(get_common_boundary(matslist,matsparslist,row_ar=row_ar,col_ar=col_ar))
        if printtocheck:
            returnar.append(folder_tocheck)
        if getallpoints:
            returnar.append(allpdf)
        return_list.append(returnar)

    if len(return_list)==1:#when all parameter have same pars_limit
        return return_list[0]
    else:#if different pars_limit
        return return_list
    
def approximate_ps_python(pars,func=None,verbose=True,min_cutoff=0.05, max_cutoff=0.95):
    """Func should compute the value of the GRF at a value x, taking as arguments the parameter vector and the x value.
    min_cutoff, max_cutoff: thresholds used to define x domain. It is expected that function values go from less than the min_cutoff to more than the max_cutoff."""
    
    #first find the appropriate x domain
    xvals=np.logspace(-20,20,200)
    y=np.array([func(pars,x) for x in xvals])
    if min(y)>min_cutoff or max(y)<max_cutoff:
        return [None,None]
    else:
        idx1=np.where(y>0.05)[0][0]
        idx2=np.where(y>0.95)[0][0]
        x1=0.01*xvals[idx1]
        x2=100*xvals[idx2]
        y1=func(pars,x1)
        y2=func(pars,x2)
        
        while y1>0.005:
            x1=x1*0.01
            y1=func(pars,x1)
        while y2<0.999:
            x2=x2*100
            y2=func(pars,x2)
        
        #get 1000 log-spaced points and find x_05 to normalise
        xvals=np.logspace(np.log10(x1),np.log10(x2),1000)
        y=[func(pars,x) for x in xvals]
        y=np.array(y)
        max_=max(y)
        argmax=np.argmax(y)
        min_=min(y[0:argmax])
        half=min_+((max_-min_)/2)
        tolerance=max_/100
        idx=np.where(np.abs(y-half)<tolerance)[0][0]
        x05=xvals[idx]
       
        
        #print("x05",x05,x05printed)
        xvalsn=xvals*x05
        #xvalsnprinted=xvals*x05printed
        yn=[func(pars,x) for x in xvalsn]
        ynmin=min(yn)
        ynmax=max(yn)

        #make sure it covers enough range of the function. this may not work for certain nonmonotonic functions
        dif=np.diff(y)
        if not np.all(dif>0):
            if verbose:
                print("The function seems not monotonic. Careful when selecting the ranges of the normalised function!")

        #assuming that the minimum is towards the beginning, extend the small values
        recompute=False
        nit=0
        if (ynmin-min_)>0.001 or (yn[0]-y[0])>0.001:
            recompute=True
            found=False
            x1n=xvalsn[0]
            while (not found) and nit<10: 
                nit+=1
                x1n=x1n*0.01
                y1n=func(pars,x1n)
                if (y1n-min_)<0.001 or (y1n-y[0])<0.001: 
                    found=True
            if nit==10:
                if verbose:
                    print("reached nit for min", xvals[0],y[0],min_, xvalsn[0],yn[0],ynmin, y1n)
                    print(y1n-min_,y1n-y[0])
                    print(",".join(map(str,pars)))
                
        else:
            x1n=xvalsn[0]
        
        nit=0
        if (max_-ynmax)>0.001 or (y[-1]-yn[-1])>0.001:
            recompute=True
            found=False
            x2n=xvalsn[-1]
            while (not found) and nit<10:
                nit+=1
                x2n=x2n*100
                y2n=func(pars,x2n)
                if (max_-y2n)<0.001 or (y[-1]-y2n)<0.001:
                    found=True
            if nit==10:
                if verbose:
                    print("reached nit for max", xvals[-1],yn[-1], ynmax, xvalsn[-1], y2n)
                    print(y2n-min_,y2n-y[-1])
                    print(",".join(map(str,pars)))
        else:
            x2n=xvalsn[-1]
        if recompute:
            xvalsn=np.logspace(np.log10(x1n),np.log10(x2n),1000)
            xvals=xvalsn/x05
            #y=[func(pars,x) for x in xvals]
        
        yn=[func(pars,x) for x in xvalsn]
        if (min(yn)>0.1) or (max(yn)<0.9):
            if verbose:
                print("min and max should be checked")
                print("min:", min(yn), "max:", max(yn))
                print(",".join(map(str,pars)))
                print()


        dif=np.diff(yn)
        dx=np.diff(xvals)
        derivative=dif/dx
        maxdif=np.argmax(derivative)
        pos=xvals[maxdif]
        stp=derivative[maxdif]
        return [(pos,stp), [xvals,xvalsn,yn]]

