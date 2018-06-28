#include <cmath>
#include <cassert>
#include <utility>
#include <algorithm>
#include <algorithm>

#include "alphaSpline.h"
#include "kernels.h"

#include "Settings.h"

#include <iostream>

using namespace std;


//////////////////////////////////////
//Helpers
//////////////////////////////////////
static double IntegralCos(double a, double x)
{
    double p = sqrt(1 - a*a); 

    if(x == M_PI)
        return 2.*M_PI/2. / p; //In case of whole half cyrcle

    double res = 2*atan( (a+1)*tan(x/2) / p) / p;
    return res;
}

pair<double,double> Kernel::GetKerPar(double l, double lp)
{
    //cout << "Epsilon is " << eps << endl;

    double r = l / lp;
    double ker = 1./(2*M_PI) * 1./(1 + r*r);
    double par =  2*r / (1 + eps + r*r);
    return {ker, par};
}

//Only for l != lp
//Integral of  lp^2/2PI * 1/(l^2+lp^2 - 2*l*lp*cos(phi)) dphi
static double IntegralPhi(double l, double lp, double angle)
{
    assert(l != lp);

    double Int;
    if (angle < M_PI)
        Int = 2 * atan((l+lp)/(l-lp) * tan(angle/2)) / (l*l - lp*lp);
    else
        Int = M_PI / abs(l*l - lp*lp);

    return Int * lp*lp/ (2*M_PI);
}

//Only for l == lp
//Integral of  l^2/2PI * 1/(l^2+l^2 - 2*l*l*cos(phi)) dphi
//Measured from angle to M_PI
static double IntegralPhiDiag(double l, double lp, double angle)
{
    assert(l == lp);
    assert(angle <= M_PI);
    if(angle == M_PI) return 0;

    double Int = 1./tan(angle/2.);
    return Int / (4*M_PI);
}


static double GetAngle(double l, double lp, double cond)
{
    double Cos = (l*l + lp*lp - cond) / (2*l*lp);
    if(Cos <= -1) return M_PI;
    if(Cos >=  1) return 0;
    return acos(Cos);
}


double Kernel::alphaS(double l, double lp)
{
    //const double LnFreeze2 =  Settings::I().LnFreeze2;
    double LnQ2 = max(2*log(l), LnFreeze2); //Currently select l as scale

    //cout << "LnQ2 " << LnQ2 << endl;
    //cout << "LnFreeze2 " << LnFreeze2 <<" "<< alphaSpline::alphaS(LnQ2, 4)<<    endl;
    return alphaSpline::alphaS(LnQ2, 4) * 3./M_PI;
}



/*

///////////////
//BFKLplain
///////////////

double BFKLplain::OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(l == lp) return 0;
    double Exp = l*l / (lp*lp);
    double res = as / (abs(1 - Exp));
    return res;
}

double BFKLplain::DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double Exp = l*l / (lp*lp);
    double resSing = 0;
    if(l != lp) resSing = as * (- Exp/(abs(1. - Exp)) );
    double resReg  = as * (+ Exp/sqrt(4 + Exp*Exp));
    //resReg  = 0;
    return (resSing + resReg);
}


///////////////
//BFKL
///////////////


//Form of equation with phi
double BFKL::OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);
    return as * ker *  Int; 
}

//Form of bfkl with phi
double BFKL::DiagEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));
    
    double Int = 2 * IntegralCos(par, angleMax);
    return -as * ker * Int;

}


//Form of equation with phi
double BFKL::OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if (l == lp) return 0;
    double Int = 2 * IntegralPhi(l, lp, M_PI); //for uper and lower "hemisphere"
    return as * Int; 
}

//Form of bfkl with phi
double BFKL::DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;
    
    double Int;
    if (l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l,lp, angleMax);

    return -as * Int;
}


///////////////
//BFKL_res (with cutoff on mu)
///////////////

//Off-diagonal kernel with F(k+q)
double BFKL_res::OffEps(double l, double lp, double z)
{
    const double mu2 = Settings::I().mu2;

    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double BFKL_res::DiagEps(double l, double lp, double z)
{
    const double mu2 = Settings::I().mu2;
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));

    double Int = 2 * IntegralCos(par, angleMax);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return -as * ker * Int;

}


//Off-diagonal kernel with F(k+q)
double BFKL_res::OffSub(double l, double lp, double z)
{
    const double mu2 = Settings::I().mu2;
    double as = alphaS(l, lp);
    if (l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        return as * 2 * IntegralPhiDiag(l, lp, angleMin);
    }

    double Int = 2 * IntegralPhi(l, lp, M_PI);

    // theta(q2 - mu2) term
    if (pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return as * Int; 
}

//Diagonal kernel with F(k)
double BFKL_res::DiagSub(double l, double lp, double z)
{
    const double mu2 = Settings::I().mu2;
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;

    if(l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        if(angleMax <= angleMin) return 0;
        return -2*as * (IntegralPhiDiag(l, lp, angleMin) - IntegralPhiDiag(l, lp, angleMax));
    }

    double Int = 2 * IntegralPhi(l, lp, angleMax);
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return -as * Int;
}



////////////////////////////////////////////////
//BFKL_res_kc_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_simp::OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = M_PI;
    // theta(k2/q2 - z)
    double Cos = (lp*lp + l*l -l*l/z) / (2*l*lp);
    if(Cos >= 1) return 0;
    if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    //cout << ker << " "<< Int << " "<< angleMax <<" "<< Cos << " "<< l<<" "<<lp<<" "<< l*l/pow(l-lp,2) -z<<  endl;
    //assert(isfinite(as*ker*Int));

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double BFKL_res_kc_simp::DiagEps(double l, double lp, double z) {
    return BFKL::DiagEps(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_simp::OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l/z);
    if(angleMax == 0) return 0;

    double Int;
    if(l != lp)
        Int = 2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as * Int; 
}

//Diagonal kernel with F(k)
double BFKL_res_kc_simp::DiagSub(double l, double lp, double z) {
    return BFKL::DiagSub(l, lp, z);
}

////////////////////////////////////////////////
//BFKL_res_kc_v_r_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_v_r_simp::OffEps(double l, double lp, double z)
{
    return BFKL_res_kc_simp::OffEps(l, lp, z);
}

//Diagonal kernel with F(k)
double BFKL_res_kc_v_r_simp::DiagEps(double l, double lp, double z) {
    return BFKL_res_kc_simp::DiagEps(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_v_r_simp::OffSub(double l, double lp, double z)
{
    return BFKL_res_kc_simp::OffSub(l, lp, z);
}

//Diagonal kernel with F(k)
double BFKL_res_kc_v_r_simp::DiagSub(double l, double lp, double z) {
    return BFKL_res_kc_simp::DiagSub(l, lp, z);
}



////////////////////////////////////////////////
//BFKL_res_kc_full
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_full::OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    if(z == 1) return 0;
    double angleMax = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    return as * ker *  Int; 
}


//Diagonal kernel with F(k)
double BFKL_res_kc_full::DiagEps(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL::DiagEps(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_full::OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(z == 1) return 0;

    // theta(k2/q2 - z/(1-z))
    double angleMax = GetAngle(l,lp, l*l*(1-z)/z);

    double Int;
    if(l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as *  Int; 
}


//Diagonal kernel with F(k)
double BFKL_res_kc_full::DiagSub(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL::DiagSub(l, lp, z);
}


////////////////////////////////////////////////
//BFKL_res_kc_v_r_full
////////////////////////////////////////////////




//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_v_r_full::OffEps(double l, double lp, double z)
{
    return BFKL_res_kc_full::OffEps(l, lp, z);
}

//Diagonal kernel with F(k)
double BFKL_res_kc_v_r_full::DiagEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;
    if(z == 1) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax1 = acos(lp/(2.*l));
    
    double angleMax2 = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax2 = acos(Cos);

    double Int = 2 * IntegralCos(par, min(angleMax1, angleMax2));
    return -as * ker * Int;
}

//Off-diagonal kernel with F(k+q)
double BFKL_res_kc_v_r_full::OffSub(double l, double lp, double z)
{
    return BFKL_res_kc_full::OffSub(l, lp, z);
}

//Diagonal kernel with F(k)
double BFKL_res_kc_v_r_full::DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(z == 1) return 0;

    double angleMax1 = GetAngle(l,lp, l*l);
    if(angleMax1 == 0) return 0;

    double angleMax2 = GetAngle(l,lp, l*l*(1-z)/z);
    if(angleMax2 == 0) return 0;

    double Int;
    if (l != lp)
        Int = 2 * IntegralPhi(l, lp, min(angleMax1, angleMax2));
    else
        Int = -2 * IntegralPhiDiag(l, lp, min(angleMax1, angleMax2));

    return -as * Int;
}

*/


////////////////////////////////////////////////
//BFKL_res_DGLAP 
////////////////////////////////////////////////

//Return z/6*Pgg - 1
static double PggMod(double z)
{
    double reg = z*z*(1-z) - (2*z+1);
    //return reg;
    //reg = 0;
    if(z == 1) return reg;
    else       return reg + 1./(1-z);
}

double Kernel::DGLAPterm(double l, double lp, double z)
{
    const double mu2 = Settings::I().mu2;
    double as = alphaS(l, lp);
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return 0.;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }
    double res = PggMod(z) * as * ker * Int;

    return res;
}

static double Harmonic(double stepSize, int nStep)
{
    double rapNow = 0;
    double sum = 0;
    for(int i = 0; i < nStep; ++i, rapNow -= stepSize) {
        if(i == 0) continue;//skip z == 1
        double step = (i == nStep - 1) ? stepSize/2 : stepSize;
        double z = exp(rapNow);
        sum += -z/(1-z) * step;
        //sum += -1/(1-z) * step;
        //cout << "RADEK " << x <<" "<< z << endl;
    }
    return sum;
}

/*

//Off-diagonal kernel with F(k+q)
double BFKL_res_DGLAP::OffEps(double l, double lp, double z)
{
    double res = BFKL::OffEps(l, lp, z);
    res += DGLAPterm(l, lp, z);
    return res;
}

//Diagonal kernel with F(k)
double BFKL_res_DGLAP::DiagEps(double l, double lp, double z)
{
    return BFKL::DiagEps(l, lp, z);
}



//Off-diagonal z-diaginal kernel with F(k+q)
double BFKL_res_DGLAP::zDiagEps(double l, double lp, double x)
{
    const double mu2 = Settings::I().mu2;
    const double rapMin = Settings::I().rapMin;
    const double rapMax = Settings::I().rapMax;
    const double Nrap = Settings::I().Nrap;

    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);

    static double harVals[5000];
    static bool isInit = false;
    if(!isInit) {
        for(int nStep = 0; nStep < 5000; ++nStep)
            harVals[nStep] = Harmonic(stepSize, nStep);
        isInit = true;
    }


    double as = alphaS(l, lp);
    Settings::I().putZero = true;
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return 0;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }

    const int nf = 4;

    double rap = log(1/x);
    double nStepReal = (rap-rapMin)/(rapMax-rapMin) * (Nrap - 1) + 1;
    int nStep = round(nStepReal);
    assert(abs(nStepReal-nStep) < 1e-2);

    double sum = harVals[nStep];
    //cout << "" << sum / harVals[nStep] << endl;

    //assert(x != 1);
    //cout << "x val = " << x << endl;
    sum += (33 - 2*nf) / 36.0;

    if(nStep != 1)
        sum += log(1 - x/exp(-rapMin));
    else //For the begining?
        sum += log(1 - exp(-stepSize));

    double res = as * sum * ker * Int;

    return res;
}


////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP
////////////////////////////////////////////////


double BFKL_res_kc_full_DGLAP::OffEps(double l, double lp, double z)
{
    //Adding normal BFKL
    double res =  BFKL_res_kc_full::OffEps(l, lp, z);
    res += DGLAPterm(l, lp, z);

    return res;
}

double BFKL_res_kc_full_DGLAP::DiagEps(double l, double lp, double z)
{
    //return 0;
    return BFKL::DiagEps(l, lp, z);
}

//Off-diagonal z-diaginal kernel with F(k+q)
double BFKL_res_kc_full_DGLAP::zDiagEps(double l, double lp, double x)
{
    return BFKL_res_DGLAP::zDiagEps(l, lp, x);
}

*/


//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//New way
//////////////////////////////////////////////////////////////////////


//BFKLplain
///////////////

double Kernel::BFKLplain__OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    //assert(!isnan(as));
    if(l == lp) return 0;
    double Exp = l*l / (lp*lp);
    //assert(!isnan(Exp));
    double res = as / (abs(1 - Exp));
    //assert(!isnan(res));
    return res;
}

double Kernel::BFKLplain__DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double Exp = l*l / (lp*lp);
    double resSing = 0;
    if(l != lp) resSing = as * (- Exp/(abs(1. - Exp)) );
    double resReg  = as * (+ Exp/sqrt(4 + Exp*Exp));
    //resReg  = 0;
    return (resSing + resReg);
}

double Kernel::BFKLplain__zDiagSub(double, double, double) {return 0;}


double Kernel::BFKLplain__OffEps(double, double, double)
{return NotImpleneted(); }

double Kernel::BFKLplain__DiagEps(double, double, double)
{return NotImpleneted(); }

double Kernel::BFKLplain__zDiagEps(double, double, double)
{return NotImpleneted(); }




///////////////
//BFKL
///////////////


//Form of equation with phi
double Kernel::BFKL__OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    //assert(isfinite(as));
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);

    //assert(isfinite(Int));
    //assert(isfinite(ker));

    return as * ker *  Int; 
}

//Form of bfkl with phi
double Kernel::BFKL__DiagEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));
    
    double Int = 2 * IntegralCos(par, angleMax);


    return -as * ker * Int;

}

double Kernel::BFKL__zDiagEps(double l, double lp, double z) {return 0;}


//Form of equation with phi
double Kernel::BFKL__OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if (l == lp) return 0;
    double Int = 2 * IntegralPhi(l, lp, M_PI); //for uper and lower "hemisphere"
    return as * Int; 
}

//Form of bfkl with phi
double Kernel::BFKL__DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;
    
    double Int;
    if (l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l,lp, angleMax);

    return -as * Int;
}

double Kernel::BFKL__zDiagSub(double l, double lp, double z) {return 0;}

///////////////
//BFKL_res (with cutoff on mu)
///////////////

//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res__OffEps(double l, double lp, double z)
{

    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res__DiagEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));

    double Int = 2 * IntegralCos(par, angleMax);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return -as * ker * Int;

}
double Kernel::BFKL_res__zDiagEps(double l, double lp, double z) {return 0;}


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res__OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if (l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        return as * 2 * IntegralPhiDiag(l, lp, angleMin);
    }

    double Int = 2 * IntegralPhi(l, lp, M_PI);

    // theta(q2 - mu2) term
    if (pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return as * Int; 
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res__DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;

    if(l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        if(angleMax <= angleMin) return 0;
        return -2*as * (IntegralPhiDiag(l, lp, angleMin) - IntegralPhiDiag(l, lp, angleMax));
    }

    double Int = 2 * IntegralPhi(l, lp, angleMax);
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return -as * Int;
}

double Kernel::BFKL_res__zDiagSub(double l, double lp, double z) {return 0;}


////////////////////////////////////////////////
//BFKL_res_kc_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_simp__OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = M_PI;
    // theta(k2/q2 - z)
    double Cos = (lp*lp + l*l -l*l/z) / (2*l*lp);
    if(Cos >= 1) return 0;
    if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    //cout << ker << " "<< Int << " "<< angleMax <<" "<< Cos << " "<< l<<" "<<lp<<" "<< l*l/pow(l-lp,2) -z<<  endl;
    //assert(isfinite(as*ker*Int));

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_simp__DiagEps(double l, double lp, double z) {
    return BFKL__DiagEps(l, lp, z);
}

double Kernel::BFKL_res_kc_simp__zDiagEps(double l, double lp, double z) {return 0;}



//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_simp__OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double angleMax = GetAngle(l,lp, l*l/z);
    if(angleMax == 0) return 0;

    double Int;
    if(l != lp)
        Int = 2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as * Int; 
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_simp__DiagSub(double l, double lp, double z) {
    return BFKL__DiagSub(l, lp, z);
}

double Kernel::BFKL_res_kc_simp__zDiagSub(double l, double lp, double z) {return 0;}

////////////////////////////////////////////////
//BFKL_res_kc_v_r_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_simp__OffEps(double l, double lp, double z)
{
    return BFKL_res_kc_simp__OffEps(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_simp__DiagEps(double l, double lp, double z) {
    return BFKL_res_kc_simp__DiagEps(l, lp, z);
}

double Kernel::BFKL_res_kc_v_r_simp__zDiagEps(double l, double lp, double z) {return 0;}


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_simp__OffSub(double l, double lp, double z)
{
    return BFKL_res_kc_simp__OffSub(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_simp__DiagSub(double l, double lp, double z) {
    return BFKL_res_kc_simp__DiagSub(l, lp, z);
}

double Kernel::BFKL_res_kc_v_r_simp__zDiagSub(double l, double lp, double z) {return 0;}


////////////////////////////////////////////////
//BFKL_res_kc_full
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full__OffEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    if(z == 1) return 0;
    double angleMax = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    return as * ker *  Int; 
}


//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_full__DiagEps(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL__DiagEps(l, lp, z);
}

double Kernel::BFKL_res_kc_full__zDiagEps(double l, double lp, double z) {return 0;}

//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full__OffSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(z == 1) return 0;

    // theta(k2/q2 - z/(1-z))
    double angleMax = GetAngle(l,lp, l*l*(1-z)/z);

    double Int;
    if(l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as *  Int; 
}


//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_full__DiagSub(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL__DiagSub(l, lp, z);
}

double Kernel::BFKL_res_kc_full__zDiagSub(double l, double lp, double z) {return 0;}

////////////////////////////////////////////////
//BFKL_res_kc_v_r_full
////////////////////////////////////////////////




//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_full__OffEps(double l, double lp, double z)
{
    return Kernel::BFKL_res_kc_full__OffEps(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_full__DiagEps(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(lp > 2*l) return 0;
    if(z == 1) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax1 = acos(lp/(2.*l));
    
    double angleMax2 = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax2 = acos(Cos);

    double Int = 2 * IntegralCos(par, min(angleMax1, angleMax2));
    return -as * ker * Int;
}

double Kernel::BFKL_res_kc_v_r_full__zDiagEps(double l, double lp, double z) {return 0;}


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_full__OffSub(double l, double lp, double z)
{
    return Kernel::BFKL_res_kc_full__OffSub(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_full__DiagSub(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if(z == 1) return 0;

    double angleMax1 = GetAngle(l,lp, l*l);
    if(angleMax1 == 0) return 0;

    double angleMax2 = GetAngle(l,lp, l*l*(1-z)/z);
    if(angleMax2 == 0) return 0;

    double Int;
    if (l != lp)
        Int = 2 * IntegralPhi(l, lp, min(angleMax1, angleMax2));
    else
        Int = -2 * IntegralPhiDiag(l, lp, min(angleMax1, angleMax2));

    return -as * Int;
}

double Kernel::BFKL_res_kc_v_r_full__zDiagSub(double l, double lp, double z) {return 0;}


////////////////////////////////////////////////
//BFKL_res_DGLAP 
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_DGLAP__OffEps(double l, double lp, double z)
{
    double res = BFKL__OffEps(l, lp, z);
    res += DGLAPterm(l, lp, z);
    return res;
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_DGLAP__DiagEps(double l, double lp, double z)
{
    return BFKL__DiagEps(l, lp, z);
}



//Off-diagonal z-diaginal kernel with F(k+q)
double Kernel::BFKL_res_DGLAP__zDiagEps(double l, double lp, double x)
{

    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);

    static double harVals[5000];
    static bool isInit = false;
    if(!isInit) {
        for(int nStep = 0; nStep < 5000; ++nStep)
            harVals[nStep] = Harmonic(stepSize, nStep);
        isInit = true;
    }


    double as = alphaS(l, lp);
    putZero = true;
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return 0;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }

    const int nf = 4;

    double rap = log(1/x);
    double nStepReal = (rap-rapMin)/(rapMax-rapMin) * (Nrap - 1) + 1;
    int nStep = round(nStepReal);
    assert(abs(nStepReal-nStep) < 1e-2);

    double sum = harVals[nStep];
    //cout << "" << sum / harVals[nStep] << endl;

    //assert(x != 1);
    //cout << "x val = " << x << endl;
    sum += (33 - 2*nf) / 36.0;

    if(nStep != 1)
        sum += log(1 - x/exp(-rapMin));
    else //For the begining?
        sum += log(1 - exp(-stepSize));

    double res = as * sum * ker * Int;

    return res;
}

double Kernel::BFKL_res_DGLAP__OffSub(double l, double lp, double z)
{ return NotImpleneted(); }
double Kernel::BFKL_res_DGLAP__DiagSub(double l, double lp, double z)
{ return NotImpleneted(); }
double Kernel::BFKL_res_DGLAP__zDiagSub(double l, double lp, double x)
{ return NotImpleneted(); }




////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP
////////////////////////////////////////////////


double Kernel::BFKL_res_kc_full_DGLAP__OffEps(double l, double lp, double z)
{
    //Adding normal BFKL
    double res =  BFKL_res_kc_full__OffEps(l, lp, z);
    res += DGLAPterm(l, lp, z);

    return res;
}

double Kernel::BFKL_res_kc_full_DGLAP__DiagEps(double l, double lp, double z)
{
    //return 0;
    return BFKL__DiagEps(l, lp, z);
}

double Kernel::BFKL_res_kc_full_DGLAP__zDiagEps(double l, double lp, double x)
{ return BFKL_res_DGLAP__zDiagEps(l, lp, x); }

double Kernel::BFKL_res_kc_full_DGLAP__OffSub(double l, double lp, double z)
{ return NotImpleneted(); }
double Kernel::BFKL_res_kc_full_DGLAP__DiagSub(double l, double lp, double z)
{ return NotImpleneted(); }
double Kernel::BFKL_res_kc_full_DGLAP__zDiagSub(double l, double lp, double x)
{ return NotImpleneted(); }
