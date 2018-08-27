#ifndef _KERNELS_
#define _KERNELS_

#include <iostream>
#include <armadillo>
using namespace std;

inline int NotImpleneted()
{
    cout <<"The selected kernel with given regularisation not implemented" << endl;
    exit(1); 
    return 0;
}


struct Kernel {


#define printOne(name,tag) double name##__##tag##_Off  (double l, double lp, double z); \
                           double name##__##tag##_Diag (double l, double lp, double z);
#define printTwo(name)  printOne(name,Eps); printOne(name,Sub);



    arma::mat myCheck;

    double mu2, eps;
    double rapMin, rapMax;
    int Nrap;
    bool putZero = false;
    double LnFreeze2;
    vector<double> SqrtExpXi;

    double alphaS(double l, double lp);
    pair<double,double> GetKerPar(double l, double lp);

    double DGLAPfactorOld(double l, double lp);
    double DGLAPfactorNew(double l, double lp);
    double DGLAPfactorExtra(double l, double lp, double z, double a);

    double DGLAPoffOld(double l, double lp, double z);
    double DGLAPoffOldEps(double l, double lp, double z);
    double DGLAPoffNew(double l, double lp, double z);
    double DGLAPoffNewEps(double l, double lp, double z);
    double DGLAPoffExtra(double l, double lp, double z, double a);
    double DGLAPoffExtraEps(double l, double lp, double z, double a);

    double DGLAPdiagOldEps(double l, double lp);
    double DGLAPdiagOld(double l, double lp);
    double DGLAPdiagNew(double l, double lp);
    double DGLAPdiagNewEps(double l, double lp);
    double DGLAPdiagExtra(double l, double lp, double z, double a);
    double DGLAPdiagExtraEps(double l, double lp, double z, double a);

    double DGLAPextraCommon(double l, double lp, double z, double a);
    double DGLAPextraCommonEps(double l, double lp, double z, double a);


    double getRapIndexReal(double z) const;
    int getRapIndex(double z) const;
    double PggModSing(double z, double a = 1) const;
    double PggModSingEps(double z, double a = 1) const;
    double Harmonic(double startRap, double stepSize, double a, int nStep) const;
    double HarmonicEps(double startRap, double stepSize, double a, int nStep) const;

    printOne(BFKLplain,Sub)
    printTwo(BFKL)
    printTwo(BFKL_res)
    printTwo(BFKL_res_kc_simp)
    printTwo(BFKL_res_kc_v_r_simp)
    printTwo(BFKL_res_kc_full)
    printTwo(BFKL_res_kc_v_r_full)
    printOne(BFKL_res_DGLAP,Eps)
    printOne(BFKL_res_DGLAP,ZEps)
    printOne(BFKL_res_kc_full_DGLAP,Eps)
    printOne(BFKL_res_kc_full_DGLAP_simp_kc,Eps)
    printOne(BFKL_res_kc_full_DGLAP_simp_kc,ZEps)
    printOne(BFKL_res_kc_full_DGLAP_full_kc,Eps)

};

#endif
