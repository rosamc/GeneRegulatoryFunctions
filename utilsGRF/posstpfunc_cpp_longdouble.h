#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "unipolynoboost.hpp"

using namespace std;
using namespace Eigen;

/* Function to compute position and steepness for a GRF, and function to assess if it is monotonic or not.
Rosa Martinez Corral. 06/11/2019
*/

void get_positive_real_roots_aberth(vector<long double> coeffsx05, vector<long double> &pos_real){
    
    reverse(coeffsx05.begin(),coeffsx05.end()); //here first coefficient is larger degree
    utils::UniPoly<long double> coeffsx05poly = utils::UniPoly<long double>("x",coeffsx05); //this is implemented by Chris in unipolynoboost

    

    // Get the half-maximal concentration of the GRF
    std::pair<std::vector<mp_complex<long double> >, utils::SolverStats<long double> > solution_data;
    const long double imag_tol = 1e-20;
    solution_data = coeffsx05poly.roots("aberth", imag_tol);
    auto roots_halfmax = solution_data.first;
    long double x05;
    mp_complex<long double> root;
    //cout << "with Chris's method";

    std::vector<long double>::size_type i;
    for (i=0;i<roots_halfmax.size();i++){
        root=roots_halfmax[i];
        //cout<<root<<"\n";

            if (abs(root.imag()) < imag_tol && root.real() > 0.0){
                pos_real.push_back(root.real());
            } 
    }
}



void get_positive_real_roots(vector<long double> &coeffs, vector<long double> &pos_real){
    //Returns roots of polynomial given by coeffs that are positive and real. Coeffs should be in increasing order. 
    //Computes eigenvalues of companion matrix. It is not always accurate.

    //vector<double>coeffs;
    std::vector<long double>::size_type Nc=coeffs.size();
    std::vector<long double>::size_type i,j, Nrows, lastidx;
    
    Nrows=Nc-1;
    lastidx=Nrows-1;

    Matrix< long double, Dynamic, Dynamic> A(Nrows,Nrows);
    A.setZero(Nrows,Nrows);
    //MatrixXd A = MatrixXd::Zero(Nrows,Nrows); //This would be the same. Matrix with dynamic allocation of rows and columns, filled with 0
    //But do not use the code below. It produces inaccurate results!!
    //for (i=1;i<Nrows;i++){
    //    for (j=0;j<Nrows;j++){
    //    A(i,j)=0; // 
    //    }
    //}
    
    for (i=1;i<Nrows;i++){
        A(i,i-1)=1; // fill diagonal (-1 diagonal)
    }
    for (i=0;i<Nrows;i++){
        A(i,lastidx)=-1 * coeffs[i]/coeffs[Nrows]; //fill coefficients
    }

    //VectorXcd rootsout;
    Matrix<complex<long double>,Dynamic,1> rootsout;
    EigenSolver<Matrix<long double, Dynamic, Dynamic>> es(A, false); //no need to compute eigenvectors
    rootsout=es.eigenvalues();

    //vector<double> pos_real;
    //py::print("The roots out are");
    //for (i=0;i<rootsout.size();i++){
    //    py::print(rootsout[i]);
    //}

    if (rootsout.size()>0){
    for (i=0;i<rootsout.size();i++){
        //py::print("Checking root",rootsout[i]);
        if (rootsout[i].imag() ==0){
            if (rootsout[i].real()>0){
                pos_real.push_back(rootsout[i].real());
            }     
        }
        }
    }
    
    //py::print("The roots found are");
    //for (i=0;i<pos_real.size();i++){
    //    py::print(pos_real[i]);
    //}

}


void product_coeffs_vectors(vector<long double> &c1, vector<long double>  &c2, vector<long double> &result){
    //vector<double> result(degree+degree);
    std::vector<long double>::size_type idx, i, j, n1, n2;
    n1=c1.size();
    n2=c2.size();
    for ( i=0;i<n1;i++){
        for ( j=0;j<n2;j++){
            idx=i+j;      
            result[idx]+=c1[i]*c2[j];
        }
    }
} 


void remove_zeros_endvector(vector<long double> &v){ //, vector<double> &vout){
    //py::print("removing zeros");
    //cout << "removing zeros";

    while (abs(v.back())<1e-20){
        v.pop_back();
    }
    
    /*std::vector<double>::size_type i;
    std::vector<double>::size_type j;
    std::vector<double>::size_type Nc;
    //int i,j,Nc;
    Nc=vin.size();
    //py::print(Nc-1,vin[Nc-1]);

    //std::vector<double>::reverse_iterator rit = vin.rbegin();
    if (abs(vin[Nc-1])<1e-15){
        for (i = Nc-1;i>=0;i--){
            if (abs(vin[i])>1e-15){
                for (j=0;j<=i;j++){
                    vout.push_back(vin[j]);
                }
                //py::print("new Nc",Nc);
                break;
            }
        }
    }else{
        for (j=0;j<Nc;j++){
                    vout.push_back(vin[j]);
                }
        //vout=vin;
    }
    //py::print("length after removing",vout.size());
    */
   

}

void get_fraction_derivative_coeffs(vector<long double> &c1, vector<long double> &c2, vector<long double> &derivativenum, vector<long double> &derivativeden){
    //Computes derivative of the fraction given by c1/c2 in terms of the coefficients. Returns the coefficients of the numerator and the denominator of the derivative.
    //c1 contains coefficients of numerator, c2 contains the coefficients of the denominator. Both c1 and c2 should have 0 in the case of 0 coefficients before the coefficient with last degree. E.g. 0,0,x^2,0,x^3 but not 0,0,x^2,0,x^3,0."""
    
    std::vector<long double>::size_type nnum, nden, degreenum, degreeden, degreesum, i, j;
    nnum=c1.size();
    nden=c2.size();
    degreenum=nnum-1;
    degreeden=nden-1;
    degreesum=degreenum+degreeden;

    //printf("size nnum is %lu, nden is %lu\n",nnum,nden);
    //printf("size der vector is %lu, dder vector is %lu\n",derivativenum.size(),derivativeden.size());
    

    //compute coefficients of derivative of the numerator and denominator polynomials
    vector<long double> dc1(degreenum);
    vector<long double> dc2(degreeden);
    
    //start at 1 so that derivative of term of degree 0 disappears
    for (i=1;i<nnum;i++){
        dc1[i-1]=i*c1[i];
        //printf("dc1 %u, %Le , %Le",i, c1[i], dc1[i-1]);     
    }
    for (i=1;i<nden;i++){
        dc2[i-1]=i*c2[i];
    }

    
    vector<long double> dc1_dot_c2(degreesum); //derivative of numerator is degreenum-1, and denominator is degreeden. Need to add one for the degree 0. 
    product_coeffs_vectors(dc1,c2,dc1_dot_c2); //derivative of numerator * denominator
    vector<long double> dc2_dot_c1(degreesum);
    product_coeffs_vectors(dc2,c1,dc2_dot_c1); //derivative of denominator * numerator
    //vector<long double> densq(degreesum+1,0); //denominator^2;
    product_coeffs_vectors(c2,c2,derivativeden);

    /*
    py::print("dc1dot");

    for(i=0;i<dc1_dot_c2.size();i++){
        py::print(dc1_dot_c2[i]);
        }
    for(i=0;i<dc2_dot_c1.size();i++){
        py::print(dc2_dot_c1[i]);
        }
    */
    //vector<double> dernum;
    for (i=0;i<degreesum;i++){     
        derivativenum[i]=(dc1_dot_c2[i]-dc2_dot_c1[i]);
    }

}

vector<double> compute_pos_stp(vector<long double> &num, vector<long double> &den, string rootmethod){


    std::vector<long double>::size_type i, j, nnum, nden;
    long double i1;
    long double i2;

    nnum=num.size();
    nden=den.size(); 
    

    vector<long double> coeffsx05(max(nnum,nden)); 
    for (i=0;i<max(nnum,nden);i++){
        if (i<nnum){
            i1=num[i];
        }else{
            i1=0;
        }
        if (i<nden){
            i2=0.5*den[i];
        }else{
            i2=0;
        }

        coeffsx05[i]=i1-i2;
    }

    //cout << "Printing GRF num\n";
    //for (i=0;i<max(nnum,nden);i++){
    //    cout << num[i];
    //    cout << "\n";
    //}


     //Find roots of coeffsx05, and keep the ones that are positive and real. Throw an error if this happens more than once. 
    
    long double x05;
    vector<long double> x05v;
    if (rootmethod=="simple"){
    get_positive_real_roots(coeffsx05, x05v);
    }
    if (rootmethod=="aberth"){
    get_positive_real_roots_aberth(coeffsx05, x05v);
    }
    double maxder=-1;
    double xmaxder=-1;


    if (x05v.size()!=1){
        x05=-1;

    }else{
        x05=x05v[0];
    }
    //py::print("x05 is ", x05);
    //printf("x05: %Le\n",x05);
    

    double x05double = (double) x05;

    
    if (x05>0.0){

        //Normalise coefficients of num and den so that hopefully numbers become more reasonable
        
        
        for (i=0;i<nnum;i++){
            num[i]=num[i]*pow(x05,(int) i);
            //printf("term %d, power x05 %Le, num[i] %Le\n",i, pow(x05,i), num[i]);
        }

        for (i=0;i<nden;i++){
            den[i]=den[i]*pow(x05,(int) i);
        }
        
        long double xp, num_sum, den_sum;
       
        std::vector<long double>::size_type Niter;

        //derivative of GRF. It is a fraction: derivativenum/derivativeden
        //GRF numerator has degree nnum-1 and denominator has degree nden-1. So den*d(num)/dx -num*d(den)/dx is of degree nden-1+nnum-1-1. 
        // So the array has to be of size nnum+nden-3+1 to account for 0.
        //GRF denominator has degree nden-1. So squared of this is degree nden-1+nden-1. 

        vector<long double> derivativenum(nnum+nden-2); //(nsize)
        vector<long double> derivativeden(nden+nden-1);

        get_fraction_derivative_coeffs(num,den,derivativenum,derivativeden);
        remove_zeros_endvector(derivativenum); //,derivativenum);
        remove_zeros_endvector(derivativeden); //_,derivativeden);

        nnum=derivativenum.size();
        nden=derivativeden.size();

        vector<long double> derivative2num(nnum+nden-2);
        vector<long double> derivative2den(nden+nden-1);
        //cout << "going for second";

        

        //second derivative
        get_fraction_derivative_coeffs(derivativenum,derivativeden,derivative2num,derivative2den);
        
        //py::print("second");
        //cout << "done second\n";
        remove_zeros_endvector(derivative2num);
        remove_zeros_endvector(derivative2den);
        
        
        vector<long double> critpoints;
        if (rootmethod=="simple"){
        get_positive_real_roots(derivative2num,critpoints); //critical points are derivative2=0 so numerator of derivative2=0;
        }
        if (rootmethod=="aberth"){
        get_positive_real_roots_aberth(derivative2num,critpoints); //critical points are derivative2=0 so numerator of derivative2=0;
        }
        //py::print("criticalpoints");
        //cout << "critpoints\n ";
        //for(i=0;i<critpoints.size();i++){
        //    printf("%Le,", critpoints[i]);
        //}
        //cout << "\n";
        

        
        //if (critpoints.size()==0){
             //py::print("no critical points found");
        //    cout << "no critical points found\n";

        if (critpoints.size()>0){
            
            long double thirdderx;
            nnum=derivative2num.size();
            nden=derivative2den.size();


            vector<long double> derivative3num(nnum+nden-2);
            vector<long double> derivative3den(nden+nden-1);

            get_fraction_derivative_coeffs(derivative2num,derivative2den,derivative3num,derivative3den);

            remove_zeros_endvector(derivative3num);
            remove_zeros_endvector(derivative3den);
            //py::print("derivative3", derivative3num.size());
            

            vector <long double> maxderv;
            vector <long double> xmaxderv;

            
            for (j=0;j<critpoints.size();j++){
                num_sum=0.0;
                den_sum=0.0;
                nnum=derivative3num.size();
                nden=derivative3den.size();
                for (i=0;i<nnum;i++){
                    xp=pow(critpoints[j],i);
                    num_sum+=derivative3num[i]*xp;
                }
                for (i=0;i<nden;i++){
                    xp=pow(critpoints[j],i);
                    den_sum+=derivative3den[i]*xp;
                    //py::print("Adding to num",derivative3num[i],xp,num_sum);
                    //py::print("Adding to den",derivative3den[i],xp,den_sum);
                }
                thirdderx=num_sum/den_sum;
                // py::print("For point",critpoints[j],num_sum,den_sum,thirdderx);
                if (thirdderx<0){ //local maximum

                    num_sum=0;
                    den_sum=0;
                    nnum=derivativenum.size();
                    nden=derivativeden.size();
                    for (i=0;i<nnum;i++){
                        xp=pow(critpoints[j],(int)i);
                        num_sum+=derivativenum[i]*xp;
                    }
                    for (i=0;i<nden;i++){
                        xp=pow(critpoints[j],(int)i);
                        den_sum+=derivativeden[i]*xp;
                        //py::print("Adding to num",derivative3num[i],xp,num_sum);
                        //py::print("Adding to den",derivative3den[i],xp,den_sum);
                    }


                    maxderv.push_back(num_sum/den_sum);
                    xmaxderv.push_back(critpoints[j]);
                }
            }

            Niter=maxderv.size();
            if (Niter>0){
                i=distance(maxderv.begin(),max_element(maxderv.begin(),maxderv.end()));
                maxder=maxderv[i]; //*x05;
                xmaxder=xmaxderv[i]; ///x05;
                //py::print("maxder",maxder);
                //py::print("xmaxder",xmaxder);


            }//else{
                //py::print("Could not find max derivative despite finding critical points.");
                //cout << "Could not find max derivative despite finding critical points.";
            //}
        }


    }
    
    //printf("at the end: %g, %g\n", xmaxder, maxder);
    vector<double> result = {xmaxder, maxder};
    return result;
}

int compute_monotonic(vector<long double> &num, vector<long double> &den){


    std::vector<long double>::size_type i, j, nnum, nden;
    long double i1;
    long double i2;

    nnum=num.size();
    nden=den.size(); 


    //derivative of GRF. It is a fraction: derivativenum/derivativeden
    //GRF numerator has degree nnum-1 and denominator has degree nden-1. So den*d(num)/dx -num*d(den)/dx is of degree nden-1+nnum-1-1. 
    // So the array has to be of size nnum+nden-3+1 to account for 0.
    //GRF denominator has degree nden-1. So squared of this is degree nden-1+nden-1. 

    vector<long double> derivativenum(nnum+nden-2); //(nsize)
    vector<long double> derivativeden(nden+nden-1);

    get_fraction_derivative_coeffs(num,den,derivativenum,derivativeden);
    remove_zeros_endvector(derivativenum); //,derivativenum);
    remove_zeros_endvector(derivativeden); //_,derivativeden);

    //critical points of derivative correspond to roots of derivativenum
    vector<long double> critpoints;
    //get_positive_real_roots(derivative2num,critpoints); //critical points are derivative2=0 so numerator of derivative2=0;
    get_positive_real_roots_aberth(derivativenum,critpoints); //critical points are derivative2=0 so numerator of derivative2=0;

    if (critpoints.size()>0){
    	return 1;
    }else{
    	return 0;
    }
}