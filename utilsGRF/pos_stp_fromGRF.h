#include <pybind11/pybind11.h>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>

using namespace std;

long double GRF(vector<long double>&num, vector<long double>&den, double x){
	long double num_sum=0;
	long double den_sum=0;
	int nnum=num.size();
	int nden=den.size();
	int i;
    if (nnum==nden){
        double xp;
        for (i=0;i<nnum;i++){
            xp=pow(x,i);
            num_sum+=num[i]*xp;
            den_sum+=den[i]*xp;
        }

    }else{
        for (i=0;i<nnum;i++){
            num_sum+=num[i]*pow(x,i);
        }
	 
	    for (i=0;i<nden;i++){
		    den_sum+=den[i]*pow(x,i);
	    }
    }

	return num_sum/den_sum;

}
vector<double> compute_pos_stp_fromGRF(vector<long double>&num, vector<long double>&den, bool verbose=false, int npoints=1000){
    int Nmax=70;
    vector<double> xvecmax_(Nmax);
    //vector<double> xvecmax(Nmax); //10^xvecmax_
    vector<long double> outmax(Nmax);
    double xval=-15;
    double xval10;
    int i;
    vector<double>result={-1.0,-1.0};
    long double maxGRF=GRF(num,den,pow(10,20));
    long double GRF01=maxGRF*0.02;
    long double GRF99=maxGRF*0.98;
    long double GRF95=maxGRF*0.95;
    long double GRF05=maxGRF*0.5;
    double x0=-15;
    double x1=-15;
    double x5val=-15;
    double GRFval;

    //First compute GRF at points from 10^-15 to 10^20 (70 points total) to get the range of the GRF
    for (i=0;i<Nmax;i++){
    	xval10=pow(10,xval);
    	xvecmax_[i]=xval;
    	//xvecmax[i]=xval10;
    	GRFval=GRF(num,den,xval10);
    	xval+=0.5;

        if ((GRFval>=GRF01)&(x0<-14)){
            x0=xvecmax_[i]-1;
        }
        if ((GRFval>=GRF99)&(x1<-14)){
            x1=xvecmax_[i]+1;
            break;
        }

    }

    if ((x0<-14) | (x1<-14)){
    	return result;
    }else{
    	double dx;
        int norders=(int)(x1-x0+0.5);
        int N=npoints*norders; //1000 points per order magnitude
        dx=(x1-x0)/N;

    	//vector<double>xvec_(N);
    	vector<double>xvec(N);
    	vector<double>out(N);
    	double x=x0;
    	for (i=0;i<N;i++){
    		//xvec_[i]=x;
    		xval10=pow(10,x);
    		xvec[i]=xval10;
    		out[i]=GRF(num,den,xval10);
    		x+=dx;
    	}

    	double halfX=-1;
    	long double der=0;
    	long double maxder=0;
    	double xmaxder=0;

    	for (i=0;i<N-1;i++){
    		der=(out[i+1]-out[i])/(xvec[i+1]-xvec[i]);
    		if (der>maxder){
    			maxder=der;
    		    xmaxder=xvec[i];
    		}else if(out[i]>GRF95){
                break;
            }
            if ((out[i]>=GRF05)&(halfX<0)){
                halfX=xvec[i];
            }
    	}
        result={xmaxder/halfX,(double)maxder*halfX};
        return result;

    }

}