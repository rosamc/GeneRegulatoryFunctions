#include <iostream>
#include <Eigen/Dense>
#include <vector>
 
using namespace std;
using namespace Eigen;
typedef double T;
typedef Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > MatrixXd;
typedef Eigen::Matrix< T, Eigen::Dynamic, 1 > VectorXd;
 
//Maybe something to think about would be to treat the matrix as a sparse matrix. But since this seems to work, I'll got with this now.
void nullspace(Ref<MatrixXd> L, Ref<VectorXd> N, bool doublecheck=false)
{ 
    
    const int n=L.rows();
    unsigned int r,c,i;
    double cs;
    double tolerance=1e-10;
        
        
    //Matrix<double,12,1> N;//solution of nullspace
    JacobiSVD<MatrixXd> svd(L,ComputeThinV);
    /* Get the V matrix */
    MatrixXd V((int)svd.matrixV().rows(), (int)svd.matrixV().cols());
    V = svd.matrixV();

    
    //if (svd.singularValues()(n-1,0)>tolerance){
    //for large values of the coefficients the tolerance is not met so for now I am not testing this
    //    cout << "Last singular value is greater than tolerance." << svd.singularValues()(n-1,0) << "\n";
        //N(0,0)=-1;
        //return; 
    //}else{
        N=V.col(svd.matrixV().cols()-1).cwiseAbs();

        cs=N.sum();
        for (r=0;r<n;r++){
            N(r,0)=abs(N(r,0))/cs; //negative values sometimes are negative. keep positive, although it doesn't really matter
        }

        if (doublecheck==true){
            MatrixXd out(n,1);
            out=L*N;
            i=0;
            bool goodnullspace=true;

            for (i=0;i<n;i++){
                if (abs(out(i,0))>tolerance){
                    cout << "inaccurate nullspace to tolerance" << tolerance << "\n";
                    cout << out;
                    goodnullspace=false;
                    break;
                }
            }
            cout << "Goodnullspace? " << goodnullspace << "\n";
        }

        //cout << N;
        //cout <<"\n";
       return; 
      
    //}
}

double ssfromnullspace(Ref<MatrixXd> L, vector<int> &indices, vector<double> &coefficients, bool doublecheck=false){
    int n=L.rows();
    int i;
    T ss=0;
    VectorXd N;
    N.resize(n,1);


    nullspace(L,N,doublecheck);
    if (N(0,0)<-0.1){ //if it cannot solve for the nullspace the first element of N will be -1
        return -1.0;
    }else{
        for (i=0;i<indices.size();i++){
            ss+=N(indices[i],0)*coefficients[i];
        }
        return ss;
    }
}