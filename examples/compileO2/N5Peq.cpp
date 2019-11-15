
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
namespace py=pybind11;

void GRF_N5P_x(py::array_t<double> parsar, vector<long double> &num, vector<long double> &den, py::array_t<double>othervars){
    typedef long double T;

    auto parsarbuf=parsar.request();
    double *pars=(double *) parsarbuf.ptr;
    T K1=pars[0];
    T K2=pars[1];
    T K3=pars[2];
    T K4=pars[3];
    T K5=pars[4];
    T w12=pars[5];
    T w13=pars[6];
    T w14=pars[7];
    T w15=pars[8];
    T w23=pars[9];
    T w24=pars[10];
    T w25=pars[11];
    T w34=pars[12];
    T w35=pars[13];
    T w45=pars[14];
    T w123=pars[15];
    T w124=pars[16];
    T w125=pars[17];
    T w134=pars[18];
    T w135=pars[19];
    T w145=pars[20];
    T w234=pars[21];
    T w235=pars[22];
    T w245=pars[23];
    T w345=pars[24];
    T w1234=pars[25];
    T w1235=pars[26];
    T w1245=pars[27];
    T w1345=pars[28];
    T w2345=pars[29];
    T w12345=pars[30];
    T wp1=pars[31];
    T wp2=pars[32];
    T wp3=pars[33];
    T wp4=pars[34];
    T wp5=pars[35];
    T wp12=pars[36];
    T wp13=pars[37];
    T wp14=pars[38];
    T wp15=pars[39];
    T wp23=pars[40];
    T wp24=pars[41];
    T wp25=pars[42];
    T wp34=pars[43];
    T wp35=pars[44];
    T wp45=pars[45];
    T wp123=pars[46];
    T wp124=pars[47];
    T wp125=pars[48];
    T wp134=pars[49];
    T wp135=pars[50];
    T wp145=pars[51];
    T wp234=pars[52];
    T wp235=pars[53];
    T wp245=pars[54];
    T wp345=pars[55];
    T wp1234=pars[56];
    T wp1235=pars[57];
    T wp1245=pars[58];
    T wp1345=pars[59];
    T wp2345=pars[60];
    T wp12345=pars[61];

    auto varsarbuf=othervars.request();
    double *varsar=(double *) varsarbuf.ptr;
    T P=varsar[0];
    vector<T> coeffs_1={1};
    vector<T> coeffs_2={0, K1};
    vector<T> coeffs_3={0, K2};
    vector<T> coeffs_4={0, K3};
    vector<T> coeffs_5={0, K4};
    vector<T> coeffs_6={0, K5};
    vector<T> coeffs_7={P};
    vector<T> coeffs_8={0, 0, K1*K2*w12};
    vector<T> coeffs_9={0, 0, K1*K3*w13};
    vector<T> coeffs_10={0, 0, K1*K4*w14};
    vector<T> coeffs_11={0, 0, K1*K5*w15};
    vector<T> coeffs_12={0, K1*P*wp1};
    vector<T> coeffs_13={0, 0, K2*K3*w23};
    vector<T> coeffs_14={0, 0, K2*K4*w24};
    vector<T> coeffs_15={0, 0, K2*K5*w25};
    vector<T> coeffs_16={0, K2*P*wp2};
    vector<T> coeffs_17={0, 0, K3*K4*w34};
    vector<T> coeffs_18={0, 0, K3*K5*w35};
    vector<T> coeffs_19={0, K3*P*wp3};
    vector<T> coeffs_20={0, 0, K4*K5*w45};
    vector<T> coeffs_21={0, K4*P*wp4};
    vector<T> coeffs_22={0, K5*P*wp5};
    vector<T> coeffs_23={0, 0, 0, K1*K2*K3*w123*w23};
    vector<T> coeffs_24={0, 0, 0, K1*K2*K4*w124*w24};
    vector<T> coeffs_25={0, 0, 0, K1*K2*K5*w125*w25};
    vector<T> coeffs_26={0, 0, K1*K2*P*w12*wp12};
    vector<T> coeffs_27={0, 0, 0, K1*K3*K4*w134*w34};
    vector<T> coeffs_28={0, 0, 0, K1*K3*K5*w135*w35};
    vector<T> coeffs_29={0, 0, K1*K3*P*w13*wp13};
    vector<T> coeffs_30={0, 0, 0, K1*K4*K5*w145*w45};
    vector<T> coeffs_31={0, 0, K1*K4*P*w14*wp14};
    vector<T> coeffs_32={0, 0, K1*K5*P*w15*wp15};
    vector<T> coeffs_33={0, 0, 0, K2*K3*K4*w234*w34};
    vector<T> coeffs_34={0, 0, 0, K2*K3*K5*w235*w35};
    vector<T> coeffs_35={0, 0, K2*K3*P*w23*wp23};
    vector<T> coeffs_36={0, 0, 0, K2*K4*K5*w245*w45};
    vector<T> coeffs_37={0, 0, K2*K4*P*w24*wp24};
    vector<T> coeffs_38={0, 0, K2*K5*P*w25*wp25};
    vector<T> coeffs_39={0, 0, 0, K3*K4*K5*w345*w45};
    vector<T> coeffs_40={0, 0, K3*K4*P*w34*wp34};
    vector<T> coeffs_41={0, 0, K3*K5*P*w35*wp35};
    vector<T> coeffs_42={0, 0, K4*K5*P*w45*wp45};
    vector<T> coeffs_43={0, 0, 0, 0, K1*K2*K3*K4*w1234*w234*w34};
    vector<T> coeffs_44={0, 0, 0, 0, K1*K2*K3*K5*w1235*w235*w35};
    vector<T> coeffs_45={0, 0, 0, K1*K2*K3*P*w123*w23*wp123};
    vector<T> coeffs_46={0, 0, 0, 0, K1*K2*K4*K5*w1245*w245*w45};
    vector<T> coeffs_47={0, 0, 0, K1*K2*K4*P*w124*w24*wp124};
    vector<T> coeffs_48={0, 0, 0, K1*K2*K5*P*w125*w25*wp125};
    vector<T> coeffs_49={0, 0, 0, 0, K1*K3*K4*K5*w1345*w345*w45};
    vector<T> coeffs_50={0, 0, 0, K1*K3*K4*P*w134*w34*wp134};
    vector<T> coeffs_51={0, 0, 0, K1*K3*K5*P*w135*w35*wp135};
    vector<T> coeffs_52={0, 0, 0, K1*K4*K5*P*w145*w45*wp145};
    vector<T> coeffs_53={0, 0, 0, 0, K2*K3*K4*K5*w2345*w345*w45};
    vector<T> coeffs_54={0, 0, 0, K2*K3*K4*P*w234*w34*wp234};
    vector<T> coeffs_55={0, 0, 0, K2*K3*K5*P*w235*w35*wp235};
    vector<T> coeffs_56={0, 0, 0, K2*K4*K5*P*w245*w45*wp245};
    vector<T> coeffs_57={0, 0, 0, K3*K4*K5*P*w345*w45*wp345};
    vector<T> coeffs_58={0, 0, 0, 0, 0, K1*K2*K3*K4*K5*w12345*w2345*w345*w45};
    vector<T> coeffs_59={0, 0, 0, 0, K1*K2*K3*K4*P*w1234*w234*w34*wp1234};
    vector<T> coeffs_60={0, 0, 0, 0, K1*K2*K3*K5*P*w1235*w235*w35*wp1235};
    vector<T> coeffs_61={0, 0, 0, 0, K1*K2*K4*K5*P*w1245*w245*w45*wp1245};
    vector<T> coeffs_62={0, 0, 0, 0, K1*K3*K4*K5*P*w1345*w345*w45*wp1345};
    vector<T> coeffs_63={0, 0, 0, 0, K2*K3*K4*K5*P*w2345*w345*w45*wp2345};
    vector<T> coeffs_64={0, 0, 0, 0, 0, K1*K2*K3*K4*K5*P*w12345*w2345*w345*w45*wp12345};
    T numdeg0=1*coeffs_12[0]+1*coeffs_16[0]+1*coeffs_19[0]+1*coeffs_21[0]+1*coeffs_22[0]+1*coeffs_26[0]+1*coeffs_29[0]+1*coeffs_31[0]+1*coeffs_32[0]+1*coeffs_35[0]+1*coeffs_37[0]+1*coeffs_38[0]+1*coeffs_40[0]+1*coeffs_41[0]+1*coeffs_42[0]+1*coeffs_45[0]+1*coeffs_47[0]+1*coeffs_48[0]+1*coeffs_50[0]+1*coeffs_51[0]+1*coeffs_52[0]+1*coeffs_54[0]+1*coeffs_55[0]+1*coeffs_56[0]+1*coeffs_57[0]+1*coeffs_59[0]+1*coeffs_60[0]+1*coeffs_61[0]+1*coeffs_62[0]+1*coeffs_63[0]+1*coeffs_64[0];
    T numdeg1=1*coeffs_12[1]+1*coeffs_16[1]+1*coeffs_19[1]+1*coeffs_21[1]+1*coeffs_22[1]+1*coeffs_26[1]+1*coeffs_29[1]+1*coeffs_31[1]+1*coeffs_32[1]+1*coeffs_35[1]+1*coeffs_37[1]+1*coeffs_38[1]+1*coeffs_40[1]+1*coeffs_41[1]+1*coeffs_42[1]+1*coeffs_45[1]+1*coeffs_47[1]+1*coeffs_48[1]+1*coeffs_50[1]+1*coeffs_51[1]+1*coeffs_52[1]+1*coeffs_54[1]+1*coeffs_55[1]+1*coeffs_56[1]+1*coeffs_57[1]+1*coeffs_59[1]+1*coeffs_60[1]+1*coeffs_61[1]+1*coeffs_62[1]+1*coeffs_63[1]+1*coeffs_64[1];
    T numdeg2=1*coeffs_26[2]+1*coeffs_29[2]+1*coeffs_31[2]+1*coeffs_32[2]+1*coeffs_35[2]+1*coeffs_37[2]+1*coeffs_38[2]+1*coeffs_40[2]+1*coeffs_41[2]+1*coeffs_42[2]+1*coeffs_45[2]+1*coeffs_47[2]+1*coeffs_48[2]+1*coeffs_50[2]+1*coeffs_51[2]+1*coeffs_52[2]+1*coeffs_54[2]+1*coeffs_55[2]+1*coeffs_56[2]+1*coeffs_57[2]+1*coeffs_59[2]+1*coeffs_60[2]+1*coeffs_61[2]+1*coeffs_62[2]+1*coeffs_63[2]+1*coeffs_64[2];
    T numdeg3=1*coeffs_45[3]+1*coeffs_47[3]+1*coeffs_48[3]+1*coeffs_50[3]+1*coeffs_51[3]+1*coeffs_52[3]+1*coeffs_54[3]+1*coeffs_55[3]+1*coeffs_56[3]+1*coeffs_57[3]+1*coeffs_59[3]+1*coeffs_60[3]+1*coeffs_61[3]+1*coeffs_62[3]+1*coeffs_63[3]+1*coeffs_64[3];
    T numdeg4=1*coeffs_59[4]+1*coeffs_60[4]+1*coeffs_61[4]+1*coeffs_62[4]+1*coeffs_63[4]+1*coeffs_64[4];
    T numdeg5=1*coeffs_64[5];
    T dendeg0=1*coeffs_1[0]+1*coeffs_2[0]+1*coeffs_3[0]+1*coeffs_4[0]+1*coeffs_5[0]+1*coeffs_6[0]+1*coeffs_7[0]+1*coeffs_8[0]+1*coeffs_9[0]+1*coeffs_10[0]+1*coeffs_11[0]+1*coeffs_12[0]+1*coeffs_13[0]+1*coeffs_14[0]+1*coeffs_15[0]+1*coeffs_16[0]+1*coeffs_17[0]+1*coeffs_18[0]+1*coeffs_19[0]+1*coeffs_20[0]+1*coeffs_21[0]+1*coeffs_22[0]+1*coeffs_23[0]+1*coeffs_24[0]+1*coeffs_25[0]+1*coeffs_26[0]+1*coeffs_27[0]+1*coeffs_28[0]+1*coeffs_29[0]+1*coeffs_30[0]+1*coeffs_31[0]+1*coeffs_32[0]+1*coeffs_33[0]+1*coeffs_34[0]+1*coeffs_35[0]+1*coeffs_36[0]+1*coeffs_37[0]+1*coeffs_38[0]+1*coeffs_39[0]+1*coeffs_40[0]+1*coeffs_41[0]+1*coeffs_42[0]+1*coeffs_43[0]+1*coeffs_44[0]+1*coeffs_45[0]+1*coeffs_46[0]+1*coeffs_47[0]+1*coeffs_48[0]+1*coeffs_49[0]+1*coeffs_50[0]+1*coeffs_51[0]+1*coeffs_52[0]+1*coeffs_53[0]+1*coeffs_54[0]+1*coeffs_55[0]+1*coeffs_56[0]+1*coeffs_57[0]+1*coeffs_58[0]+1*coeffs_59[0]+1*coeffs_60[0]+1*coeffs_61[0]+1*coeffs_62[0]+1*coeffs_63[0]+1*coeffs_64[0];
    T dendeg1=1*coeffs_2[1]+1*coeffs_3[1]+1*coeffs_4[1]+1*coeffs_5[1]+1*coeffs_6[1]+1*coeffs_8[1]+1*coeffs_9[1]+1*coeffs_10[1]+1*coeffs_11[1]+1*coeffs_12[1]+1*coeffs_13[1]+1*coeffs_14[1]+1*coeffs_15[1]+1*coeffs_16[1]+1*coeffs_17[1]+1*coeffs_18[1]+1*coeffs_19[1]+1*coeffs_20[1]+1*coeffs_21[1]+1*coeffs_22[1]+1*coeffs_23[1]+1*coeffs_24[1]+1*coeffs_25[1]+1*coeffs_26[1]+1*coeffs_27[1]+1*coeffs_28[1]+1*coeffs_29[1]+1*coeffs_30[1]+1*coeffs_31[1]+1*coeffs_32[1]+1*coeffs_33[1]+1*coeffs_34[1]+1*coeffs_35[1]+1*coeffs_36[1]+1*coeffs_37[1]+1*coeffs_38[1]+1*coeffs_39[1]+1*coeffs_40[1]+1*coeffs_41[1]+1*coeffs_42[1]+1*coeffs_43[1]+1*coeffs_44[1]+1*coeffs_45[1]+1*coeffs_46[1]+1*coeffs_47[1]+1*coeffs_48[1]+1*coeffs_49[1]+1*coeffs_50[1]+1*coeffs_51[1]+1*coeffs_52[1]+1*coeffs_53[1]+1*coeffs_54[1]+1*coeffs_55[1]+1*coeffs_56[1]+1*coeffs_57[1]+1*coeffs_58[1]+1*coeffs_59[1]+1*coeffs_60[1]+1*coeffs_61[1]+1*coeffs_62[1]+1*coeffs_63[1]+1*coeffs_64[1];
    T dendeg2=1*coeffs_8[2]+1*coeffs_9[2]+1*coeffs_10[2]+1*coeffs_11[2]+1*coeffs_13[2]+1*coeffs_14[2]+1*coeffs_15[2]+1*coeffs_17[2]+1*coeffs_18[2]+1*coeffs_20[2]+1*coeffs_23[2]+1*coeffs_24[2]+1*coeffs_25[2]+1*coeffs_26[2]+1*coeffs_27[2]+1*coeffs_28[2]+1*coeffs_29[2]+1*coeffs_30[2]+1*coeffs_31[2]+1*coeffs_32[2]+1*coeffs_33[2]+1*coeffs_34[2]+1*coeffs_35[2]+1*coeffs_36[2]+1*coeffs_37[2]+1*coeffs_38[2]+1*coeffs_39[2]+1*coeffs_40[2]+1*coeffs_41[2]+1*coeffs_42[2]+1*coeffs_43[2]+1*coeffs_44[2]+1*coeffs_45[2]+1*coeffs_46[2]+1*coeffs_47[2]+1*coeffs_48[2]+1*coeffs_49[2]+1*coeffs_50[2]+1*coeffs_51[2]+1*coeffs_52[2]+1*coeffs_53[2]+1*coeffs_54[2]+1*coeffs_55[2]+1*coeffs_56[2]+1*coeffs_57[2]+1*coeffs_58[2]+1*coeffs_59[2]+1*coeffs_60[2]+1*coeffs_61[2]+1*coeffs_62[2]+1*coeffs_63[2]+1*coeffs_64[2];
    T dendeg3=1*coeffs_23[3]+1*coeffs_24[3]+1*coeffs_25[3]+1*coeffs_27[3]+1*coeffs_28[3]+1*coeffs_30[3]+1*coeffs_33[3]+1*coeffs_34[3]+1*coeffs_36[3]+1*coeffs_39[3]+1*coeffs_43[3]+1*coeffs_44[3]+1*coeffs_45[3]+1*coeffs_46[3]+1*coeffs_47[3]+1*coeffs_48[3]+1*coeffs_49[3]+1*coeffs_50[3]+1*coeffs_51[3]+1*coeffs_52[3]+1*coeffs_53[3]+1*coeffs_54[3]+1*coeffs_55[3]+1*coeffs_56[3]+1*coeffs_57[3]+1*coeffs_58[3]+1*coeffs_59[3]+1*coeffs_60[3]+1*coeffs_61[3]+1*coeffs_62[3]+1*coeffs_63[3]+1*coeffs_64[3];
    T dendeg4=1*coeffs_43[4]+1*coeffs_44[4]+1*coeffs_46[4]+1*coeffs_49[4]+1*coeffs_53[4]+1*coeffs_58[4]+1*coeffs_59[4]+1*coeffs_60[4]+1*coeffs_61[4]+1*coeffs_62[4]+1*coeffs_63[4]+1*coeffs_64[4];
    T dendeg5=1*coeffs_58[5]+1*coeffs_64[5];
    num={numdeg0,numdeg1,numdeg2,numdeg3,numdeg4,numdeg5};
    den={dendeg0,dendeg1,dendeg2,dendeg3,dendeg4,dendeg5};
}

py::array_t<double> interfaceps_s_GRF_N5P_x(py::array_t<double> parsar, py::array_t<double> othervars ) {
    typedef long double T;

    vector<T> num;
    vector<T> den;
    vector<double>result;
    GRF_N5P_x(parsar,num,den,othervars);

    result=compute_pos_stp(num,den,"simple");
    py::array_t<double> resultpy = py::array_t<double>(2);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    ptrresultpy[0]=result[0];
    ptrresultpy[1]=result[1];

    return  resultpy;
    }

py::array_t<double> interfaceps_a_GRF_N5P_x(py::array_t<double> parsar, py::array_t<double> othervars ) {
    typedef long double T;

    vector<T> num;
    vector<T> den;
    vector<double>result;
    GRF_N5P_x(parsar,num,den,othervars);

    result=compute_pos_stp(num,den,"aberth");
    py::array_t<double> resultpy = py::array_t<double>(2);
    py::buffer_info bufresultpy = resultpy.request();
    double *ptrresultpy=(double *) bufresultpy.ptr;
    ptrresultpy[0]=result[0];
    ptrresultpy[1]=result[1];

    return  resultpy;
    }

int interfacemonotonic_GRF_N5P_x(py::array_t<double> parsar, py::array_t<double> othervars ) {
    typedef long double T;

    vector<T> num;
    vector<T> den;
    int result;
    GRF_N5P_x(parsar,num,den,othervars);

    result=compute_monotonic(num,den); //return either 1 or 0

    return  result;
    }

double interface_GRF_N5P_x(py::array_t<double> parsar, py::array_t<double> othervars, double varGRFval ) {

    typedef long double T;

    vector<T> num;
    vector<T> den;
    
    GRF_N5P_x(parsar,num,den,othervars);

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
    }

PYBIND11_MODULE(N5Peq, m) {
    m.def("interfaceps_s_GRF_N5P_x", &interfaceps_s_GRF_N5P_x, "A function which returns pos stp, roots with eigenvalues of companion matrix.");
    m.def("interfaceps_a_GRF_N5P_x", &interfaceps_a_GRF_N5P_x, "A function which returns pos stp, roots with aberth method.");
    m.def("interfacemonotonic_GRF_N5P_x", &interfacemonotonic_GRF_N5P_x, "A function which assessess whether GRF has a local maximum.");
    m.def("interface_GRF_N5P_x", &interface_GRF_N5P_x, " A function that returns GRF at a given input value.");
}
