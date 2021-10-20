#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <polynomial.hpp>
#include <posstpfunc_cpp_boost_multiprT.hpp>

using namespace std;

namespace py=pybind11;

template <typename T, typename Tpoly, typename Tpolyc, typename polytype, typename thresholdtype>
class GRFCalculations{

	protected:
	    vector<T> num;
		vector<T> den;

    public:

		unsigned int precision;
		//unsigned int precision_poly;
		

		//GRFCalculations(unsigned int &precision, unsigned int &precision_poly) :  precision(precision), precision_poly(precision_poly) {};
	    //GRFCalculations() :  precision(precision) {};

        void set_precision(unsigned int prec){
	        precision=prec;
	    }

	    unsigned int get_precision(){
	        return precision;
	    }
	    
	    void fill_num_den_class(py::array_t<double> parsar, py::array_t<double>othervars);
	    
		void print_num_den(){
	    	cout.precision(precision);
	    	cout << num[0] <<"," << num[1] << "\n";
	    	cout << den[0] <<"," << den[1] << "\n";
	    	cout.flush();
	    }

	    double interfaceGRF(double varGRFval ){

	    	T result;
	    	T varGRFval_T=varGRFval;
            result=GRFatxonly<T>(num,den,varGRFval_T);

            return result.template convert_to<double>(); //ammend so it also works with long double T

        }

        py::array_t<double> interface_return_num_den(int ncoeffs){

        	py::array_t<double> resultpy = py::array_t<double>(2*ncoeffs);
		    py::buffer_info bufresultpy = resultpy.request();
		    double *ptrresultpy=(double *) bufresultpy.ptr;
		    

		    for (int i=0;i<num.size();i++){
		        ptrresultpy[i]=num[i].template convert_to<double>();
		    }
		    int j=num.size();
		    for (int i=0;i<den.size();i++){
		        ptrresultpy[i+j]=den[i].template convert_to<double>();
		    }
		    resultpy.resize({2,7});

            return resultpy;

        }

        py::array_t<double> interfaceps(bool verbose=false, double thresholdimag_=1e-15, bool writex05coeffs=false, bool absder=false, bool normalisefirst=true, string fnamecoeffs="filename.txt") {

        
            vector<double>result;
            vector<T> min_max(2);
            compute_min_maxGRF<T>(num,den,min_max);

            T Gmin=min_max[0];
            T Gmax=min_max[1];
            T midpoint=Gmin+0.5*(Gmax-Gmin);

            
            if  (Gmax<0){
                result={-1.0,-1.0,-1.0};
            }else{
            	thresholdtype thresholdimag=thresholdimag_;
                result=compute_pos_stp<T,Tpoly,Tpolyc,polytype,thresholdtype>(num,den,"aberth",verbose,midpoint,thresholdimag,writex05coeffs, absder, normalisefirst,fnamecoeffs);

            }

            py::array_t<double> resultpy = py::array_t<double>(3);
            py::buffer_info bufresultpy = resultpy.request();
            double *ptrresultpy=(double *) bufresultpy.ptr;
            ptrresultpy[0]=result[0];
            ptrresultpy[1]=result[1];
            ptrresultpy[2]=result[2];

            return  resultpy;
        }

        py::array_t<double> interfacemonotonic(double thresholdimag_=1e-15) {
    
            vector<double> result;
            thresholdtype thresholdimag=thresholdimag_;
            result=compute_monotonic<T,Tpoly,Tpolyc,polytype,thresholdtype>(num,den,thresholdimag); //return {-1} if derivative is 0, {-2} if no roots for the derivative of the GRF, -3 for each root out of the 10^-15,10^15 range, and the roots otherwise
            int n=result.size();
            py::array_t<double> resultpy = py::array_t<double>(n);
		    py::buffer_info bufresultpy = resultpy.request();
		    double *ptrresultpy=(double *) bufresultpy.ptr;
		    for (int i=0;i<n;i++){
		        //py::print("Result is",result[i]);
		        if (result[i]<-0.5){
		        ptrresultpy[i]=result[i];
		        }else{
		        if ((result[i]<pow(10.0,15))&&(result[i]>pow(10.0,-15))){
		            ptrresultpy[i]=result[i];
		        }else{
		            ptrresultpy[i]=-3;
		        }
		        }
		    }
            return resultpy;
        }

};
