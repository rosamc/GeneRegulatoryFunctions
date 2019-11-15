import sys, os, sympy, re, glob
from sympy.parsing.sympy_parser import parse_expr
from sympy.printing.cxxcode import cxxcode
import numpy as np
import subprocess
import itertools

class PrepareFiles():
    """Parent class in the hierarchy, with common functionality to write to make calculations for graphs at equilibrium and away.
    It is intended to write files to compute the value of a GRF at a given input value,
    and position-steepness of the GRF, using pybind to call c++ code. In addition, it writes mathematica files to double check.
    Initialisation arguments:
    - varGRF:variable with respect to which the GRF is computed
    - concvars: all variables that can change concentration (including varGRF.) For example, rN the Pol model, this would be ['x','P']
    - parlist: list of parameter values. Leave empty for the equilibrium models.

    """ 


    def __init__(self, varGRF='x',concvars=[''],parlist=[]):
        self.varGRF=varGRF #variable with respect to which the GRF is computed
        self.concvars=concvars #all variables that can change concentraiton
        self.parlist=parlist #rates
        self.numsympy=None
        self.densympy=None
        self.coeffs_num=None
        self.coeffs_den=None
        self.coeffs_num_string=None
        self.coeffs_den_string=None
        self.numstring=None
        self.denstring=None
       
    def simpify_num_den(self,numstr,denstr):
        """Given the strings with the mathematical expressions for numerator and denominator of the GRF (numstr, denstr), computes their corresponding sympy expressions. 
        This can be useful to do further symbolic calculation, and is used by compute_coeffs function.
        It is likely to fail for large graphs. So it is probably best to use simpify_rhos and parse_GRF functions
        """
        
        self.numstring=numstr
        self.denstring=denstr

        for label in self.parlist:
            #print(label)
            exec("%s=sympy.symbols('%s')"%(label,label))
        for var in self.concvars:
            #print('symbol', label)
            exec("%s=sympy.symbols('%s')"%(var,var))

        print('Parsing into sympy expressions. This can take a while...')
        sympyexpr=parse_expr(numstr)
        self.numsympy=sympyexpr
        sympyexpr=parse_expr(denstr)
            #print('symbol', label)
        self.densympy=sympyexpr
   
        var=sympy.symbols(self.varGRF)
        self.coeffs_num=sympy.Poly(self.numsympy,var).all_coeffs()[::-1]
        self.coeffs_den=sympy.Poly(self.densympy,var).all_coeffs()[::-1]
        self.coeffs_num_string=','.join(map(cxxcode,self.coeffs_num))
        self.coeffs_den_string=','.join(map(cxxcode,self.coeffs_den))
        self.coeffs_num_fromrhos=None
        self.coeffs_den_fromrhos=None
    
    def simpify_rhos(self):
        """Given the expressions for the rhos in self.rhos, obtains the corresponding coefficients with respect to varGRF and puts them in coeffs_rhos."""
        for label in self.parlist:
            #print(label)
            exec("%s=sympy.symbols('%s')"%(label,label))
        for var in self.concvars:
            exec("%s=sympy.symbols('%s')"%(var,var))

        coeffs_rhos=[]
        #convert each rho into a polynomial and get corresponding coefficients
        for rho in self.all_rhos:
            expr=parse_expr(rho.split("=")[1])
            poly=sympy.Poly(expr,sympy.symbols(self.varGRF))
            coeffs=poly.all_coeffs()[::-1]
            coeffs_rhos.append(coeffs)
        self.coeffs_rhos=coeffs_rhos
        
    def parse_GRF(self,numexpr,denexpr):
        """The expressions should be of the form 1*(rho1)+0.5*(rho2)+0.75*(rho3)+ktan*(rho4). (Notice that even if multiplied by 1, this has to be specified.)
        Notice that the indexing of the rhos starts at 1, to be consistent with the MTT.py. 
        Is also used for the equilibrium graphs, in this case the rhos are the products along paths."""
        
        self.numstring=numexpr
        self.denstring=denexpr
        pat=re.compile("([\(\)0-9aA-zZ\./]*)\*(\([a-z0-9\+]+\))")
        patrho=re.compile("rho([0-9]+)")
        degrees_rhos=[len(x) for x in self.coeffs_rhos]
        for enum, expr in enumerate([numexpr,denexpr]):
            all_coeffs_str=[[] for i in range(max(degrees_rhos))]
            terms=pat.findall(expr)
            for term in terms:
                multiplier,expr=term
                idxs=patrho.findall(expr)
                for idx_ in idxs:
                    idx=int(idx_)-1
                    for i in range(degrees_rhos[idx]):
                        all_coeffs_str[i].append(multiplier+"*coeffs_%s[%d]"%(idx_,i))
            print(enum)
            if enum==0:
                self.coeffs_num_fromrhos=[]
                for coeff in all_coeffs_str:
                    self.coeffs_num_fromrhos.append("+".join(coeff))
            else:
                self.coeffs_den_fromrhos=[]
                for coeff in all_coeffs_str:
                    self.coeffs_den_fromrhos.append("+".join(coeff))


    def __write_header(self, fh):
        """Writes includes of cpp file."""
        fh.write(    
"""
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>
#include <unsupported/Eigen/Polynomials>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "posstpfunc_cpp_longdouble.h"
using namespace std;
using namespace Eigen;
namespace py=pybind11;\n
""")
        
    def __write_GRF_coeffs(self,fh,funcname_varGRF, typestring, additionallinespars=None):
        """write function that returns the coefficients of the numerator and denominator of the GRF (with respect to varGRF)."""
        
        if len(self.concvars)>1:
            fh.write("void %s(py::array_t<double> parsar, vector<%s> &num, vector<%s> &den, py::array_t<double>othervars){\n"%(funcname_varGRF,typestring,typestring))
        else:
            fh.write("void %s(py::array_t<double> parsar, vector<%s> &num, vector<%s> &den){\n"%(funcname_varGRF,typestring,typestring))
        fh.write("    typedef %s T;\n"%typestring)
        fh.write("""
    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;\n""")
        
        for pnum,par in enumerate(self.parlist):

            fh.write("    T %s=pars[%d];\n"%(par,pnum))

        if len(self.concvars)>1:
            fh.write("""
    auto varsarbuf=othervars.request();
    double *varsar=(double *) varsarbuf.ptr;\n""")
            i=0
            for var in self.concvars:
                if var != self.varGRF:
                    fh.write("    T %s=varsar[%d];\n"%(var,i))
                    i+=1
        if additionallinespars is not None:
            fh.write(additionallinespars)
        
        if self.coeffs_num_fromrhos is None:

            fh.write("    num={%s};\n"%self.coeffs_num_string)
            fh.write("    den={%s};\n"%self.coeffs_den_string)
        else:
            for i in range(len(self.coeffs_rhos)):
                string_=cxxcode(self.coeffs_rhos[i])
                string_=string_.replace('[','').replace(']','')
                fh.write("    vector<T> coeffs_%d={%s};\n"%(i+1,string_))
            
            for i in range(len(self.coeffs_num_fromrhos)):
                fh.write("    T numdeg%d=%s;\n"%(i,self.coeffs_num_fromrhos[i]))
            for i in range(len(self.coeffs_den_fromrhos)):
                fh.write("    T dendeg%d=%s;\n"%(i,self.coeffs_den_fromrhos[i]))
            
            
            fh.write("    num={")
            for i in range(len(self.coeffs_num_fromrhos)):
                if i<len(self.coeffs_num_fromrhos)-1:
                    end=','
                else:
                    end='};\n'
                fh.write("numdeg%d%s"%(i,end))
                
            fh.write("    den={")
            for i in range(len(self.coeffs_den_fromrhos)):
                if i<len(self.coeffs_den_fromrhos)-1:
                    end=','
                else:
                    end='};\n'
                fh.write("dendeg%d%s"%(i,end))
            
        fh.write("}\n\n")
        
    def __write_interface_posstp_simple(self,fh,funcname_varGRF,typestring):
        """write function which is called from python that computes position and steepness using the cpp code. 
        Position and steepness are calculated computing the critical points of the derivative of the GRF.
        Derivatives are computed by successively taking products of the coefficients of the polynomials with the exponents.
        In this case roots of polynomials to find x05 and critical points of derivative are computed using eigenvalues of companion matrix."""

        if len(self.concvars)>1:
            fh.write("py::array_t<double> interfaceps_s_%s(py::array_t<double> parsar, py::array_t<double> othervars ) {\n"%funcname_varGRF)
        else:
            fh.write("py::array_t<double> interfaceps_s_%s(py::array_t<double> parsar ) {\n"%funcname_varGRF)

        fh.write("    typedef %s T;\n"%typestring)
        fh.write("""
    vector<T> num;
    vector<T> den;
    vector<double>result;
""")

        if len(self.concvars)>1:
            fh.write("    %s(parsar,num,den,othervars);\n"%funcname_varGRF)
        else:
            fh.write("    %s(parsar,num,den);\n"%funcname_varGRF)

        fh.write("""
    result=compute_pos_stp(num,den,"simple");
    py::array_t<double> resultpy = py::array_t<double>(2);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    ptrresultpy[0]=result[0];
    ptrresultpy[1]=result[1];

    return  resultpy;
    }\n
""")
        
    def __write_interface_posstp_aberth(self,fh,funcname_varGRF,typestring):
        """write function which is called from python that computes position and steepness using the cpp code. 
        In this case roots of polynomial are computed using aberth method implemented by Chris. Notice that this is much slower than using the "simple" option, but more accurate."""

        if len(self.concvars)>1:
            fh.write("py::array_t<double> interfaceps_a_%s(py::array_t<double> parsar, py::array_t<double> othervars ) {\n"%funcname_varGRF)
        else:
            fh.write("py::array_t<double> interfaceps_a_%s(py::array_t<double> parsar ) {\n"%funcname_varGRF)

        fh.write("    typedef %s T;\n"%typestring)
        fh.write("""
    vector<T> num;
    vector<T> den;
    vector<double>result;
""")

        if len(self.concvars)>1:
            fh.write("    %s(parsar,num,den,othervars);\n"%funcname_varGRF)
        else:
            fh.write("    %s(parsar,num,den);\n"%funcname_varGRF)

        fh.write("""
    result=compute_pos_stp(num,den,"aberth");
    py::array_t<double> resultpy = py::array_t<double>(2);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    ptrresultpy[0]=result[0];
    ptrresultpy[1]=result[1];

    return  resultpy;
    }\n
""")
        
    def __write_interface_GRFatinput(self,fh,funcname_varGRF,typestring):
        """Write function which is called from python that computes the value of the GRF at a given input value."""

        if len(self.concvars)>1:
            fh.write("double interface_%s(py::array_t<double> parsar, py::array_t<double> othervars, double varGRFval ) {\n\n"%funcname_varGRF)
        else:
            fh.write("double interface_%s(py::array_t<double> parsar, double varGRFval ) {\n\n"%funcname_varGRF)

        fh.write("    typedef %s T;\n"%typestring)
        fh.write("""
    vector<T> num;
    vector<T> den;
    
""")

        if len(self.concvars)>1:
            fh.write("    %s(parsar,num,den,othervars);\n"%funcname_varGRF)
        else:
            fh.write("    %s(parsar,num,den);\n"%funcname_varGRF)

        fh.write("""
    T numsum=0;
    T densum=0;
    std::vector<T>::size_type i;
    for (i=0;i<num.size();i++){
    numsum+=num[i]*pow(varGRFval,(int)i);
    }
    for (i=0;i<den.size();i++){
    densum+=den[i]*pow(varGRFval,(int)i);
    }
    double result=numsum/densum;
    return result;
    }\n
""")
                     
    def __write_interface_monotonic(self,fh,funcname_varGRF,typestring):
        """Write function which is called from python that checks if the GRF has critical points (returns 0) or not (returns 1).
        This is experimental and should be tested further as of November 2019."""

        if len(self.concvars)>1:
            fh.write("py::array_t<double> interfacemonotonic_%s(py::array_t<double> parsar, py::array_t<double> othervars ) {\n"%funcname_varGRF)
        else:
            fh.write("py::array_t<double> interfacemonotonic_%s(py::array_t<double> parsar ) {\n"%funcname_varGRF)

        fh.write("    typedef %s T;\n"%typestring)
        fh.write("""
    vector<T> num;
    vector<T> den;
    vector<double> result;
""")

        if len(self.concvars)>1:
            fh.write("    %s(parsar,num,den,othervars);\n"%funcname_varGRF)
        else:
            fh.write("    %s(parsar,num,den);\n"%funcname_varGRF)

        fh.write("""
    result=compute_monotonic(num,den); //return {-1} if derivative is 0, {-2} if no roots for the derivative of the GRF, and the roots otherwise
    int n=result.size();
    py::array_t<double> resultpy = py::array_t<double>(n);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    for (int i=0;i<n;i++){
        if ((result[i]<1e15)&&(result[i]>1e-15)){
            ptrresultpy[i]=result[i];
        }
    }
    return resultpy;
    
    }\n
""")

        
        
    def write_pybind_interface(self,fname=None, funcname=None, typestring='long double',additionallinespars=None):
        """Writes .cpp file with GRF functions and computation of position and steepness. pybind is used to interface with them.

        """

        var=sympy.symbols(self.varGRF)
        funcname_varGRF=funcname+'_'+self.varGRF
        
        
        #print(coeffs_num)
        #print(coeffs_den)
        f=open(fname,'w')
        self.__write_header(f)
        self.__write_GRF_coeffs(f,funcname_varGRF,typestring,additionallinespars) #coefficients of num and den of GRF with respect to input
        self.__write_interface_posstp_simple(f,funcname_varGRF,typestring) #function to call from python to compute position-steepness. roots are found from eigenvalues of companion matrix
        self.__write_interface_posstp_aberth(f,funcname_varGRF,typestring) #function to call from python to compute position-steepness. roots are found using aberth method. 
        self.__write_interface_monotonic(f,funcname_varGRF,typestring) #function to call from python to assess monotonicity
        self.__write_interface_GRFatinput(f,funcname_varGRF,typestring) #function to call from python to compute GRF at given input value
                     
        f.write("PYBIND11_MODULE(%s, m) {\n"%fname.split('/')[-1].replace('.cpp',''))
        f.write("    m.def(\"interfaceps_s_%s\", &interfaceps_s_%s, \"A function which returns pos stp, roots with eigenvalues of companion matrix.\");\n"%(funcname_varGRF, funcname_varGRF))
        f.write("    m.def(\"interfaceps_a_%s\", &interfaceps_a_%s, \"A function which returns pos stp, roots with aberth method.\");\n"%(funcname_varGRF, funcname_varGRF))
        f.write("    m.def(\"interfacemonotonic_%s\", &interfacemonotonic_%s, \"A function which assessess whether GRF has a local maximum.\");\n"%(funcname_varGRF,funcname_varGRF))
        f.write("    m.def(\"interface_%s\", &interface_%s, \" A function that returns GRF at a given input value.\");\n"%(funcname_varGRF, funcname_varGRF))
        f.write("}\n")
        f.close()
        
    def write_checkfile_mathematica_singlevar(self,fname='./testmf.txt', additionallinespars=None):
        print('writing mathematica file')
        if len(self.concvars)>1:
            if additionallinespars is None:
                print("There is more than 1 input variable. Pass a string to additionallinespars argument with the values corresponding to those")
              
        f=open(fname,'w')
        f.write("""
Print[\"First line\"];
infolder=$ScriptCommandLine[[2]];
SetDirectory[infolder];
tolpos=0.005;
tolstp=0.005;
Print[\"starting at\"];
Print[infolder];
infiles=FileNames[\"mat*.in\"];\n
        """)
        parstring_underscore=''
        for pnum,par in enumerate(self.parlist):
            par=par.replace('_','U')
            parstring_underscore+=par+'_,'
        print(parstring_underscore)

        parsvars=parstring_underscore

        for var in self.concvars:
            if var != self.varGRF:
                parsvars+=var+'_,'
        parsvars=parsvars.strip(',') #remove last comma
        allunderscore=parsvars+','+self.varGRF+'_'
        allnounderscore=allunderscore.replace('_','')
        
        if self.coeffs_num_fromrhos is None:
            f.write('num[%s]:=%s;\n'%(allunderscore,self.numstring.replace('_','U')))
            f.write('den[%s]:=%s;\n'%(allunderscore,self.denstring.replace('_','U')))
        else:
            for i in range(len(self.coeffs_rhos)):
                rho_coeffs=self.coeffs_rhos[i] #rho_i expressed as the coefficients of the polynomial
                poly_string=[]
                for d in range(len(rho_coeffs)): #for each coeffiscient
                    poly_string.append('('+str(rho_coeffs[d]).replace('_','U')+')'+"*"+self.varGRF+"^%d"%d)
                
                f.write("rho%d[%s]:=%s;\n"%(i+1,allunderscore,"+".join(poly_string)))
            
            pat=re.compile("(rho[0-9]+)")
            numstring_pars=re.sub(pat,"\g<1>[%s]"%(allnounderscore),self.numstring)
            denstring_pars=re.sub(pat,"\g<1>[%s]"%(allnounderscore),self.denstring)
            f.write("num[%s]:=%s;\n"%(allunderscore,numstring_pars))
            f.write("den[%s]:=%s;\n"%(allunderscore,denstring_pars))
        
        
        f.write('GRF[%s]:=num[%s]/den[%s];\n'%(allunderscore,allnounderscore,allnounderscore))
        parstring_nounderscore=parstring_underscore.replace('_','')
        parstring_nounderscore=parstring_nounderscore.strip(',')
        npars=len(parstring_nounderscore.split(','))
        f.write('{%s}=Table[0,{i,1,%d}];\n'%(parstring_nounderscore,npars))
        if additionallinespars is not None:
            f.write(additionallinespars)
        parstring_nounderscore_commas=["\"%s\""%x for x in parstring_nounderscore.split(',')]
        f.write("parsliststring={%s};"%(",".join(parstring_nounderscore_commas)))
        f.write("""
Print[\"Defined GRF\"];
For[j=1,j<=Length[infiles],j++, 
infname=infiles[[j]];
WriteString[\"stdout\",infname,\"\\n\"];
outfname=StringJoin[StringSplit[infname,\".\"][[1]],\"_checked.out\"];
WriteString[\"stdout\",outfname,\"\\n\"];
data=Import[infname,\"CSV\"];
outf=OpenWrite[outfname];
WriteString[outf,StringRiffle[{\"pos\",\"rho\"},","],";", StringRiffle[parsliststring,","],"\\n"];
For[i =1,i<= Length[data], i++,
{pos0, stp0}=data[[i]][[1;;2]];
PLIST=data[[i]][[3;;]];
{%s}=PLIST[[1;;Length[PLIST]]];\n"""%parstring_nounderscore)
        f.write("f[%s_]:=GRF[%s];\n"%(self.varGRF,allnounderscore))
        f.write("halfX = Solve[f[%s]==1/2&&%s>0,%s];\n"%(self.varGRF,self.varGRF,self.varGRF))
        f.write("""
g[yv2]:=f[halfX[[1]][[1]][[2]]*yv2];
maxY = Solve[D[g[yv2],{yv2,2}]==0&&yv2>0&&yv2<10^(5),yv2];
maxD = D[g[yv2],yv2]/.maxY;
rho =Max[maxD]; (* TAKE MAX GLOBAL DERIVATIVE *)
pos = maxY[[Position[maxD,Max[maxD]][[1]][[1]]]][[1]][[2]];(* TAKE POS CORRESPONDING TO MAX GLOBAL DERIVATIVE *)
If[ pos-pos0<tolpos&&rho-stp0<tolstp,
WriteString["stdout",i,","];
WriteString[outf,StringRiffle[{pos,rho},","],";", StringRiffle[PLIST,","],"\\n"];,

WriteString["stdout","\\n not correct"];
WriteString["stdout",pos,,rho,,pos0,, stp0];

]
]
Close[outf];
]
        """)
        
    def write_checksingleparset_mathematica_singlevar(self,fname='./testmf.txt', additionallinespars=None):
        """Writes a mathematica file to check the calculation for a given parameter set. Adapted from code by Felix Wong."""
        print('writing mathematica file')
        if len(self.concvars)>1:
            if additionallinespars is None:
                print("There is more than 1 input variable. Pass a string to additionallinespars argument with the values corresponding to those")
              
        f=open(fname,'w')
        parstring_underscore=''
        for pnum,par in enumerate(self.parlist):
            par=par.replace('_','U')
            parstring_underscore+=par+'_,'
        print(parstring_underscore)

        parsvars=parstring_underscore

        for var in self.concvars:
            if var != self.varGRF:
                parsvars+=var+'_,'
        parsvars=parsvars.strip(',') #remove last comma
        allunderscore=parsvars+','+self.varGRF+'_'
        allnounderscore=allunderscore.replace('_','')
        
        if self.coeffs_num_fromrhos is None:
            f.write('num[%s]:=%s;\n'%(allunderscore,self.numstring.replace('_','U')))
            f.write('den[%s]:=%s;\n'%(allunderscore,self.denstring.replace('_','U')))
        else:
            for i in range(len(self.coeffs_rhos)):
                rho_coeffs=self.coeffs_rhos[i] #rho_i expressed as the coefficients of the polynomial
                poly_string=[]
                for d in range(len(rho_coeffs)): #for each coeffiscient
                    poly_string.append('('+str(rho_coeffs[d]).replace('_','U')+')'+"*"+self.varGRF+"^%d"%d)
                
                f.write("rho%d[%s]:=%s;\n"%(i+1,allunderscore,"+".join(poly_string)))
            
            pat=re.compile("(rho[0-9]+)")
            numstring_pars=re.sub(pat,"\g<1>[%s]"%(allnounderscore),self.numstring)
            denstring_pars=re.sub(pat,"\g<1>[%s]"%(allnounderscore),self.denstring)
            f.write("num[%s]:=%s;\n"%(allunderscore,numstring_pars))
            f.write("den[%s]:=%s;\n"%(allunderscore,denstring_pars))
        
        
        f.write('GRF[%s]:=num[%s]/den[%s];\n'%(allunderscore,allnounderscore,allnounderscore))
        parstring_nounderscore=parstring_underscore.replace('_','')
        parstring_nounderscore=parstring_nounderscore.strip(',')
        npars=len(parstring_nounderscore.split(','))
        f.write('{%s}=Table[0,{i,1,%d}];\n'%(parstring_nounderscore,npars))
        if additionallinespars is not None:
            f.write(additionallinespars)
        f.write("""
Print[\"Defined GRF\"];

PLIST={};
{%s}=PLIST[[1;;Length[PLIST]]];\n"""%parstring_nounderscore)
        f.write("f[%s_]:=GRF[%s];\n"%(self.varGRF,allnounderscore))
        f.write("halfX = Solve[f[%s]==1/2&&%s>0,%s];\n"%(self.varGRF,self.varGRF,self.varGRF))
        f.write("""
g[yv2]:=f[halfX[[1]][[1]][[2]]*yv2];
maxY = Solve[D[g[yv2],{yv2,2}]==0&&yv2>0&&yv2<10^(5),yv2];
maxD = D[g[yv2],yv2]/.maxY;
rho =Max[maxD]; (* TAKE MAX GLOBAL DERIVATIVE *)
pos = maxY[[Position[maxD,Max[maxD]][[1]][[1]]]][[1]][[2]];(* TAKE POS CORRESPONDING TO MAX GLOBAL DERIVATIVE *)
{pos,rho}
Plot[g[yv2] /. {yv2 -> y0}, {y0, 0, 2}]

        """)
    

class PrepareFilesNoneq(PrepareFiles):
    """Child class to be used for graphs away from equilibrium. Initialisation arguments:
    - edges: should be a list of edges in the form [inputnode,'label',targetnode].
    - MTT folder: folder where MTT.py script is found.
    - graphbasename: name for the grpah files. Notice that it erases existing files with this basename."""

    def __init__(self,varGRF='x',concvars=[''],parlist=[],edgelist=[],MTTfolder='',graphbasename='graph'):
        super().__init__(varGRF=varGRF,concvars=concvars,parlist=parlist)
        self.edges=edgelist
        self.MTTfolder=MTTfolder
        self.graphbasename=graphbasename
        self.graphname=os.path.join(self.MTTfolder,self.graphbasename+'.txt')
        files=glob.glob(os.path.join(self.MTTfolder,self.graphbasename)+"*")
        for f in files:
            print('removing ',f)
            os.remove(f)
            
    def write_execute_parse(self):
        #write input file for MTT.py
        
        outf=open(self.graphname,'w')
        for x in self.edges:
            x_str=list(map(str,x))
            outf.write('('+','.join(x_str)+') ')
        outf.close()
        
        #execute_MTT(self):
        
        args=["python",os.path.join(self.MTTfolder,'MTT.py'),self.graphname]
        print("executing MTT")
        print(subprocess.check_call(args))
        
    
        #parse_rhos_from_graphfiles(self):
        #this was previously called rhos_from_edges
        edges=np.array(self.edges)
        #find files that correspond to the spanning trees
        files=glob.glob(os.path.join(self.MTTfolder,self.graphbasename+'*'))
        pat=re.compile('[0-9]+.txt')
        n=0
        for f in files:
            if pat.findall(f):
                n+=1
        all_rhos=[] 
        for i in range(1,n+1):
            fname=os.path.join(self.MTTfolder,self.graphbasename+'-%d.txt'%i)

            fi=open(fname,'r')
            rho_expr='rho_%d='%i

            for lnum,l in enumerate(fi.readlines()):
                idxs=np.array(list(map(int,l.strip().split())),dtype=bool)
                subset=edges[idxs]
                expr=[x[1] for x in subset] #labels that multiply together in a term of the rho
                if lnum<1:
                    rho_expr=rho_expr+'*'.join(expr)
                else:
                    rho_expr=rho_expr+'+'+'*'.join(expr)

            all_rhos.append(rho_expr)
        
        self.all_rhos=[]
        for rho in all_rhos:
            for var in self.concvars:
                var_1='-'+var
                var_2='*'+var
                if var_1 in rho:
                    #print('replacing',var_1,var_2)
                    rho=rho.replace(var_1,var_2)
                    #print(rho)
            self.all_rhos.append(rho)

class PrepareFilesEqbindingmodels(PrepareFiles):
    """Child class to be used for hypercube graphs with only TF or TF and Pol. 
    N is the number of TF binding sites (excluding polymerase in case of the Pol model).
    Initialise only with varGRF and concvars."""
    def __init__(self,varGRF='x',concvars=['x','P'], N=3, strategy='pol'):
        
        super().__init__(varGRF=varGRF,concvars=concvars)
        self.N=N
        if 'P' in concvars and strategy != "pol":
            print("For Pol model strategy should be pol. Exiting...")
            sys.exit()
        self.strategy=strategy #an, av, oom, pol
    
        
    def __get_rho(self, sites, sitep=None):
        #product of labels from reference to node, such that labels are the bare constants and cooperativities with wij, i<j
        #sites is a string with the bound sites
        #sitep is the Pol site (should be the site with largest index)
        sites_=sites.split(',')
        if len(sites_)==1:
            if sitep is None:
                rho='K%s * x**1'%(sites_[0])
            else:
                if sites_[0]==sitep:
                    rho='P'
                else:
                    rho='K%s * x**1'%(sites_[0]) #the power to the 1 is used below for parsing
        else:
            rho=''
            if (sitep is None) or (sitep not in sites_):
                niter=len(sites_)-1
                for i in range(niter):
                    i1=sites_[i]
                    #if len(sites_[i+1:])>1:
                    i2=''.join(sites_[i+1:])
                    rho+='K%s * w%s%s * '%(i1, i1,i2)
            else:
                niter=len(sites_)-2
                for i in range(niter):
                    i1=sites_[i]
                    #if len(sites_[i+1:])>1:
                    i2=''.join(sites_[i+1:-1])
                    rho+='K%s * w%s%s * '%(i1, i1,i2)
                rho+='K%s * '%(sites_[niter])

            if sites_[-1]==sitep:
                ns=len(sites_)-1
                rho+='P * wp%s * x**%d'%(''.join(sites_[:-1]),ns)
            else:
                ns=len(sites_)
                rho+='K%s * x**%d'%(sites_[-1],ns)

        return rho
    
    
    def get_rhos_and_GRFs(self, verbose=True):
        """Generate expressions for rhos and GRFs."""
        
        if 'P' in self.concvars:
            sitep=str(self.N+1)
            Niter=self.N+2
        else:
            sitep=None
            Niter=self.N+1
        
        
        symbols_str=set() #to get parlist
        self.all_rhos=['rho1=1'] #reference node
        nrho=2
        for i in range(Niter):
            combis=itertools.combinations(np.arange(1,Niter),i)
            for x in combis:
                if len(x)>0:
                    combistr=','.join(map(str,x))
                    rho=self.__get_rho(combistr,sitep=sitep)
        
                    self.all_rhos.append('rho%d=%s'%(nrho,rho))
                    #print(nrho,rho)
                    nrho+=1
                    for word in rho.split(' '):
                        symbols_str.add(word.replace('(','').replace(')','').replace('*','').replace('+',''))
    
        #get parlist, sorted such that Ks are first, omegas are second and omegasP last
        

        symbols_str_sorted=[[],[],[]]
        for i in range(3):
            for snum,s in enumerate(symbols_str):
                if len(s)>0 and not 'P' in s and not 'x' in s:
                    if i==0:
                        if 'K' in s:
                            symbols_str_sorted[i].append(s)
                    if i==1:
                        if 'w' in s and not 'wp' in s:
                            symbols_str_sorted[i].append(s)
                    elif i==2:
                        if 'wp' in s:
                            symbols_str_sorted[i].append(s)
        symbols_str_sorted_=[]
        for x in symbols_str_sorted:
            #from a stackoverflow question: python sort is stable, so you can sort in multiple passes (in reverse order):
            x.sort() 
            x.sort(key=lambda x:len(x))
            symbols_str_sorted_.extend(x)
        print(symbols_str_sorted_)
        self.parlist=symbols_str_sorted_
        print('parlist is', self.parlist)
        
        #get GRF
        
        if self.strategy=='an':
            numstring='1*(rho%d)'%(len(self.all_rhos))
        if self.strategy=='av':
            xppat=re.compile('x\*\*([0-9]*)')
            numstring=''
            for i in range(len(self.all_rhos)):
                rho=self.all_rhos[i]
                power_=xppat.findall(rho)
                if len(power_)>0:
                    for power in power_:
                        numstring+='(%s.0/%d)*(rho%d)+'%(power,self.N,i+1)
            numstring=numstring.strip('+')
        if self.strategy=='oom':
            print("strategy oom has not been implemented. Exiting...")
            sys.exit()
        if self.strategy=='pol':
            numstring='1*('
            for i in range(len(self.all_rhos)):
                rho=self.all_rhos[i]
                if 'P' in rho and 'x' in rho:
                    numstring+='rho%d+'%(i+1)
            numstring=numstring.strip('+')
            numstring+=')'
        denstring='1*('
        for i in range(len(self.all_rhos)):
            denstring+='rho%d+'%(i+1)
        denstring=denstring.strip('+')
        denstring+=')'
        
        self.simpify_rhos()
        if verbose:
            print('GRF numerator:%s'%numstring)
            print('GRF denominator:%s'%denstring)
        self.parse_GRF(numstring,denstring)
        