#Script to write .cpp file and print compilation command for O2. Modify with apprpriate paths and names 
import numpy as np
import sys, os
sys.path.append('/home/rm335/repos/sharedposstpNov19/utilsGRF/')
import writescripts
import subprocess

obj=writescripts.PrepareFilesEqbindingmodels(varGRF='x',concvars=['x','P'],N=5,strategy='pol')
obj.get_rhos_and_GRFs()

folder='./'
if not os.path.isdir(folder):
    os.mkdir(folder)
fname=os.path.join(folder,'N5Peq.cpp') #N4an will be the name of the pybind module as well
typestring='long double' 
funcname='GRF_N5P'


obj.write_pybind_interface(fname=fname, funcname=funcname, typestring='long double')
#obj.write_checkfile_mathematica_singlevar(fname=fname.replace('.cpp','.wl'),additionallinespars="P=1;\n")
#obj.write_checksingleparset_mathematica_singlevar(fname=fname.replace('.cpp','.nb'),additionallinespars="P=1;\n")

#compile. Notice that this requires eigen. Modify path accordingly. On a mac, use -undefined dynamic_lookup
filename=fname
path_1,file=os.path.split(fname)
objectnamelinux=os.path.join(path_1,file.replace('.cpp',''))
compilestringlinux="c++ -O2 -DNDEBUG -Wall -shared -std=c++11  -fPIC -I /home/rm335/libs/eigen-eigen-323c052e1731 -I /home/rm335/repos/sharedposstpNov19/utilsGRF  `python3 -m pybind11 --includes` %s -o %s`python3-config --extension-suffix`"%(filename,objectnamelinux)
compilestring=compilestringlinux
print(compilestring)

fout=open(fname.replace('.cpp','_compile.sh'),'w')
fout.write("module load gcc/6.2.0\n")
fout.write("module load python/3.6.0\n")
fout.write("source /home/rm335/myenvpy36/bin/activate\n")
fout.write(compilestring)
fout.write("\n")
fout.close()
print("\n\ncpp file is ready to be compiled. Run bash %s_compile.sh"%fname.replace('.cpp',''))
