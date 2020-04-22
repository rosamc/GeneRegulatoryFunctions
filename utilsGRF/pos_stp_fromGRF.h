#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include "solvenullspace.h"

using namespace std;

typedef double T;
typedef Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > MatrixXd;

T GRF(vector<T>&num, vector<T>&den, double x){
	T num_sum=0;
	T den_sum=0;
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


vector<double> compute_pos_stp_fromGRF(vector<T>&num, vector<T>&den, bool verbose=false, int npoints=1000){
    //this function uses the numerator and denominator of the ss, obtained from MTT
    int Nmax=70;
    vector<double> xvecmax_(Nmax);
    //vector<double> xvecmax(Nmax); //10^xvecmax_
    vector<T> outmax(Nmax);
    double xval=-15;
    double xval10;
    int i;
    vector<double>result={-1.0,-1.0};
    T maxGRF=GRF(num,den,pow(10,20));
    T GRF01=maxGRF*0.02;
    T GRF99=maxGRF*0.98;
    T GRF95=maxGRF*0.95;
    T GRF05=maxGRF*0.5;
    double x0=-15;
    double x1=-15;
    double x5val=-15;
    T GRFval;

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
    	T der=0;
    	T maxder=0;
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

T GRFfromnullspace(Ref<MatrixXd> L_, Ref<MatrixXd> Lx, double xval, vector<int> &indices, vector<double> &coefficients, bool doublecheck=false){
    //L is passed by copying it so I can modify its contents with the appropriate xval
    //Lx should contain a 1 at each position for which L has to be multiplied by xval
    if (L_.size()!=Lx.size()){
        cout << "L and Lx do not have matching sizes.";
        return -1.0;
    }
    int n=Lx.rows(); //Lx is square so this is number of rows/columns
    int i, j;
    T ss=0;
    T cs;
    MatrixXd L;
    L=L_.replicate(1,1);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            if (Lx(i,j)>0){
                L(i,j)=L_(i,j)*xval;
            }
        }
    }
    

    //diagonals
    for ( i=0;i<n;i++){ //for each col
        cs=0;
        for (j=0;j<n;j++){ //for each element
            cs+=L(j,i);
        }
        L(i,i)=-cs;
    }

    //cout << "After multiplying by x\n.";
    //cout << L;
    //cout << "\n";

    ss=ssfromnullspace(L,indices,coefficients,doublecheck);
    return ss;

}

vector<double> compute_pos_stp_fromGRF(Ref<MatrixXd> L, Ref<MatrixXd> Lx, vector<int> &indices, vector <double> &coefficients, bool verbose=false, int npoints=1000, bool doublecheck=false){
    //this function solves for the nullspace of the L to obtain the ss
    //L: Laplacian with only parameter values. This is fixed.
    //Lx: matrix with 1s at the positions that correspond to concentration dependent parameters
    //indices: indices of the nullspace vector to use for the ss computation
    //coefficients: coefficient that multiplies each of rho_i (with i in indices) to compute the ss. 
    int Nmax=70;
    vector<double> xvecmax_(Nmax);
    //vector<double> xvecmax(Nmax); //10^xvecmax_
    vector<T> outmax(Nmax);
    double xval=-15;
    double xval10;
    int i;
    vector<double>result={-1.0,-1.0};
    
    T maxGRF=GRFfromnullspace(L,Lx,pow(10,20),indices,coefficients);
    if (maxGRF<0){
        cout << "Failed when computing maxGRF\n";
        return result;
    }
    T GRF01=maxGRF*0.02;
    T GRF99=maxGRF*0.98;
    T GRF95=maxGRF*0.95;
    T GRF05=maxGRF*0.5;
    double x0=-15;
    double x1=-15;
    double x5val=-15;
    T GRFval;

    //First compute GRF at points from 10^-15 to 10^20 (70 points total) to get the range of the GRF
    for (i=0;i<Nmax;i++){
        xval10=pow(10,xval);
        xvecmax_[i]=xval;
        //xvecmax[i]=xval10;
        GRFval=GRFfromnullspace(L,Lx,xval10,indices,coefficients,doublecheck);
        xval+=0.5;
        if (GRFval<0){
        cout << "Failed when computing GRFval, " << xval10 << "\n";
        return result;
        }

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
        cout << "N is" << N;
        cout << "\n";

        //vector<double>xvec_(N);
        vector<double>xvec(N);
        vector<double>out(N);
        double x=x0;
        for (i=0;i<N;i++){
            //xvec_[i]=x;
            xval10=pow(10,x);
            xvec[i]=xval10;
            out[i]=GRFfromnullspace(L,Lx,xval10,indices,coefficients,doublecheck); // GRF(num,den,xval10);
            if (out[i]<0){
                cout << "Failed when computing out[i] at " << xval10 << "\n";
                return result;
            }
            x+=dx;
        }

        double halfX=-1;
        T der=0;
        T maxder=0;
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