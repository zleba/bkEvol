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
                          double name##__OffSub  (double l, double lp, double z); \
                          double name##__DiagSub (double l, double lp, double z); \

    double mu2, eps;
    double rapMin, rapMax;
    int Nrap;
    bool putZero = false;
    double LnFreeze2;

    double alphaS(double l, double lp);
    pair<double,double> GetKerPar(double l, double lp);

    double DGLAPfactorOld(double l, double lp);
    double DGLAPfactorNew(double l, double lp);

    double DGLAPoffOld(double l, double lp, double z);
    double DGLAPoffNew(double l, double lp, double z);

    double DGLAPdiagOld(double l, double lp);
    double DGLAPdiagNew(double l, double lp);

    double getRapIndexReal(double z);
    int getRapIndex(double z);
    double PggModSing(double z);
    double Harmonic(double startRap, double stepSize, double a, int nStep);

    printKernel(BFKLplain)
    printKernel(BFKL)
    printKernel(BFKL_res)
    printKernel(BFKL_res_kc_simp)
    printKernel(BFKL_res_kc_v_r_simp)
    printKernel(BFKL_res_kc_full)
    printKernel(BFKL_res_kc_v_r_full)
    printKernel(BFKL_res_DGLAP)
    printKernel(BFKL_res_kc_full_DGLAP)
    printKernel(BFKL_res_kc_full_DGLAP_simp_kc)
};

#endif
