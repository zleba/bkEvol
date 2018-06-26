#ifndef _KERNELS_
#define _KERNELS_

#include <iostream>
using namespace std;

inline int NotImpleneted()
{
    cout <<"The selected kernel with given regularisation not implemented" << endl;
    exit(1); 
    return 0;
}


struct Kernel {

#define printKernel(name) double name##__OffEps  (double l, double lp, double z); \
                          double name##__DiagEps (double l, double lp, double z); \
                          double name##__zDiagEps(double l, double lp, double z); \
                          double name##__OffSub  (double l, double lp, double z); \
                          double name##__DiagSub (double l, double lp, double z); \
                          double name##__zDiagSub(double l, double lp, double z);

    double mu2, eps;
    double rapMin, rapMax;
    int Nrap;
    bool putZero;
    double LnFreeze2;

    double alphaS(double l, double lp);
    double DGLAPterm(double l, double lp, double z);
    pair<double,double> GetKerPar(double l, double lp);

    printKernel(BFKLplain)
    printKernel(BFKL)
    printKernel(BFKL_res)
    printKernel(BFKL_res_kc_simp)
    printKernel(BFKL_res_kc_v_r_simp)
    printKernel(BFKL_res_kc_full)
    printKernel(BFKL_res_kc_v_r_full)
    printKernel(BFKL_res_DGLAP)
    printKernel(BFKL_res_kc_full_DGLAP)
};

#endif
