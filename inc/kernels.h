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

class BFKLplain {
    public:
        static string name()    {return "BFKLplain";}
        static string formula() {return "";}
        static double OffEps(double l, double lp, double z) {return NotImpleneted();}
        static double DiagEps(double l, double lp, double z) {return NotImpleneted();}
        static double zDiagEps(double l, double lp, double z) {return NotImpleneted();}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


class BFKL {
    public:
        static string name()    {return "BFKL";}
        static string formula() {return "F(x,k^2)=F0(k^2)+ aS* int_x^1 dz/z * int d^2q/(pi*q^2) left{{cal F}left(frac{x}{z},|{bf k}+{bf q}|^2right)-thetaleft(k_perp^2-q_perp^2right){cal F}left(frac{x}{z},{k_perp^2}right)right}";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};

class BFKL_res {
    public:
        static string name()    {return "BFKL_res";}
        static string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


class BFKL_res_kc_simp {
    public:
        string name()    {return "BFKL_res_kc_simp";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


class BFKL_res_kc_v_r_simp {
    public:
        string name()    {return "BFKL_res_kc_v_r_simp";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


class BFKL_res_kc_full {
    public:
        string name()    {return "BFKL_res_kc_full";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


/*
class BFKL_res_kc_full_zcut: public Kernel {
    public:
        string name()    {return "BFKL_res_kc_full_zcut";}
        string formula() {return "";}
        double OffEps(double l, double lp, double z) {return 0}
        double DiagEps(double l, double lp, double z)  {return 0}

        double OffSub(double l, double lp, double z)  {return 0}
        double DiagSub(double l, double lp, double z)  {return 0}
};
*/

class BFKL_res_kc_v_r_full {
    public:
        string name()    {return "BFKL_res_kc_v_r_full";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z) {return 0;}

        static double OffSub(double l, double lp, double z);
        static double DiagSub(double l, double lp, double z);
        static double zDiagSub(double l, double lp, double z) {return 0;}
};


class BFKL_res_DGLAP {
    public:
        string name()    {return "BFKL_res_DGLAP";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z);

        static double OffSub(double l, double lp, double z) {return NotImpleneted();}
        static double DiagSub(double l, double lp, double z) {return NotImpleneted();}
        static double zDiagSub(double l, double lp, double z) {return NotImpleneted();}
};


class BFKL_res_kc_full_DGLAP {
    public:
        string name()    {return "BFKL_res_kc_full_DGLAP";}
        string formula() {return "";}
        static double OffEps(double l, double lp, double z);
        static double DiagEps(double l, double lp, double z);
        static double zDiagEps(double l, double lp, double z);

        static double OffSub(double l, double lp, double z) {return NotImpleneted();}
        static double DiagSub(double l, double lp, double z) {return NotImpleneted();}
        static double zDiagSub(double l, double lp, double z) {return NotImpleneted();}
};

#endif
