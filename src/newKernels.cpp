#include <cmath>
#include <cassert>
#include <utility>
#include <algorithm>
#include <tuple>

#include "alphaSpline.h"
#include "kernels.h"

#include <iostream>
#include <armadillo>

using namespace std;

const double regFactor = 0.05;


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
    double LnQ2 = max(2*log(l), LnFreeze2); //Currently select l as scale

    return 0.2;
    //cout << "LnQ2 " << LnQ2 << endl;
    //cout << "LnFreeze2 " << LnFreeze2 <<" "<< alphaSpline::alphaS(LnQ2, 4)<<    endl;
    return alphaSpline::alphaS(LnQ2, 4) * 3./M_PI;
}


double Kernel::getRapIndexReal(double z) const
{
    double rap = -log(z);
    double nStepReal = (rap-rapMin)/(rapMax-rapMin) * (Nrap - 1);
    return nStepReal;
}

int Kernel::getRapIndex(double z) const
{
    //double rap = -log(z);
    //double nStepReal = (rap-rapMin)/(rapMax-rapMin) * (Nrap - 1);

    double nStepReal = getRapIndexReal(z);
    int nStep = round(nStepReal);
    //cout << "Z is " << z <<" " << nStep  << " "<< nStepReal << endl;
    assert(abs(nStepReal-nStep) < 1e-2);
    
    return nStep;
}





////////////////////////////////////////////////
//BFKL_res_DGLAP 
////////////////////////////////////////////////

//Return z/6*Pgg - 1 (regular)
double PggModReg(double z)
{
    if(z <= 1)
        return z*z*(1-z) - (2*z+1);
    else
        return 0;
}


//Return z/6*Pgg - 1 (singular)
double Kernel::PggModSing(double z, double a) const
{
    //int zId = getRapIndex(z);
    double BreakID = getRapIndexReal(1./a);
    double dist = getRapIndex(z) - getRapIndexReal(1./a);
    if(dist < eps) return 0;

    double reg = 1./(1-a*z);
    if     (dist < 1 + eps) reg *= 2;
    else if(dist < 2 + eps) reg *= 0.5;

    return reg;
}

//Return z/6*Pgg - 1 (singular)
double Kernel::PggModSingEps(double z, double a) const
{
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    //int zId = getRapIndex(z);
    double BreakID = getRapIndexReal(1./a);
    double dist = getRapIndex(z) - getRapIndexReal(1./a);

    return 1./sqrt(pow(1-a*z,2)+eps) * pow(max(0.,1-a*z), regFactor);


    if(dist < -eps) { return 0; }
    if(dist < eps) {z*= exp(-stepSize); }

    double reg = 1./(1-a*z);
    if(dist < eps) reg *= 1.0;

    //cout << "z " << dist <<" "<<z <<" "<< reg <<  endl;

    return reg;
}







//Get the factor in front of the raw DGLAP term
double Kernel::DGLAPfactorOld(double l, double lp)
{
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
    return as * ker * Int;
}

double Kernel::DGLAPfactorNew(double l, double lp)
{
    double as = alphaS(l, lp);
    double C = lp*lp / (l*l);
    if(lp > l) return 0;
    return as * C;
}

double Kernel::DGLAPfactorExtra(double l, double lp, double z, double a)
{
    double as = alphaS(l, lp);
    if(lp <= l) return 0;
    if(z >= 1./a) return 0;
    
    return as;
}




double Kernel::DGLAPoffOld(double l, double lp, double z)
{
    double res = (PggModSing(z) + PggModReg(z)) * DGLAPfactorOld(l, lp);
    return res;
}

double Kernel::DGLAPoffOldEps(double l, double lp, double z)
{
    double res = (PggModSingEps(z) + PggModReg(z)) * DGLAPfactorOld(l, lp);
    return res;
}


double Kernel::DGLAPoffNew(double l, double lp, double z)
{
    double res = (PggModSing(z) + PggModReg(z)) * DGLAPfactorNew(l, lp);
    return res;
}

double Kernel::DGLAPoffNewEps(double l, double lp, double z)
{
    double res = (PggModSingEps(z) + PggModReg(z)) * DGLAPfactorNew(l, lp);
    return res;
}


double Kernel::DGLAPoffExtra(double l, double lp, double z, double a)
{
    double F = DGLAPfactorExtra(l, lp, z, a);
    if(F == 0) return 0;
    double res = (PggModSing(z,a) + PggModReg(a*z)) * F;
    return res;
}

double Kernel::DGLAPoffExtraEps(double l, double lp, double z, double a)
{
    double F = DGLAPfactorExtra(l, lp, z, a);
    if(F == 0) return 0;
    double res = (PggModSingEps(z,a) + PggModReg(a*z) *pow(max(0.,1-a*z), regFactor) ) * F;
    return res;
}





double Kernel::Harmonic(double startRap, double stepSize, double a, int nStep) const //a > 1 by definition
{
    //double rapNow = startRap;
    double sum = 0;
    for(int i = 0; i < nStep; ++i) {
        double rapNow = startRap + i*stepSize;
        double z = exp(-rapNow);
        sum += -a*z* PggModSing(z,a) * stepSize;
        //sum += -z/(1-z) * step; (old method)
        //cout << "RADEK " << x <<" "<< z << endl;
    }
    return sum;
}

double Kernel::HarmonicEps(double startRap, double stepSize, double a, int nStep) const //a > 1 by definition
{
    //double rapNow = startRap;
    double sum = 0;
    for(int i = 0; i < nStep; ++i) {
        double rapNow = startRap + i*stepSize;
        double z = exp(-rapNow);
        double factor = (i == 0) ? 0.5 : 1;
        sum += -a*z* PggModSingEps(z,a) * stepSize * factor;
        //sum += -z/(1-z) * step; (old method)
        //cout << "RADEK " << x <<" "<< z << endl;
    }
    cout << "HELENKA " << sum << endl;
    return sum;
}



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

double Kernel::BFKLplain__Sub_Off(double l, double lp, double z)
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

double Kernel::BFKLplain__Sub_Diag(double l, double lp, double z)
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
double Kernel::BFKL__Eps_Off(double l, double lp, double z)
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
double Kernel::BFKL__Eps_Diag(double l, double lp, double z)
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
double Kernel::BFKL__Sub_Off(double l, double lp, double z)
{
    double as = alphaS(l, lp);
    if (l == lp) return 0;
    double Int = 2 * IntegralPhi(l, lp, M_PI); //for uper and lower "hemisphere"
    return as * Int; 
}

//Form of bfkl with phi
double Kernel::BFKL__Sub_Diag(double l, double lp, double z)
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
double Kernel::BFKL_res__Eps_Off(double l, double lp, double z)
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
double Kernel::BFKL_res__Eps_Diag(double l, double lp, double z)
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


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res__Sub_Off(double l, double lp, double z)
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
double Kernel::BFKL_res__Sub_Diag(double l, double lp, double z)
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



////////////////////////////////////////////////
//BFKL_res_kc_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_simp__Eps_Off(double l, double lp, double z)
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
double Kernel::BFKL_res_kc_simp__Eps_Diag(double l, double lp, double z) {
    return BFKL__Eps_Diag(l, lp, z);
}




//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_simp__Sub_Off(double l, double lp, double z)
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
double Kernel::BFKL_res_kc_simp__Sub_Diag(double l, double lp, double z) {
    return BFKL__Sub_Diag(l, lp, z);
}


////////////////////////////////////////////////
//BFKL_res_kc_v_r_simp
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_simp__Eps_Off(double l, double lp, double z)
{
    return BFKL_res_kc_simp__Eps_Off(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_simp__Eps_Diag(double l, double lp, double z) {
    return BFKL_res_kc_simp__Eps_Diag(l, lp, z);
}



//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_simp__Sub_Off(double l, double lp, double z)
{
    return BFKL_res_kc_simp__Sub_Off(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_simp__Sub_Diag(double l, double lp, double z) {
    return BFKL_res_kc_simp__Sub_Diag(l, lp, z);
}



////////////////////////////////////////////////
//BFKL_res_kc_full
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full__Eps_Off(double l, double lp, double z)
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
double Kernel::BFKL_res_kc_full__Eps_Diag(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL__Eps_Diag(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full__Sub_Off(double l, double lp, double z)
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
double Kernel::BFKL_res_kc_full__Sub_Diag(double l, double lp, double z) {
    if(z == 1) return 0;
    return BFKL__Sub_Diag(l, lp, z);
}


////////////////////////////////////////////////
//BFKL_res_kc_v_r_full
////////////////////////////////////////////////




//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_v_r_full__Eps_Off(double l, double lp, double z)
{
    return Kernel::BFKL_res_kc_full__Eps_Off(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_full__Eps_Diag(double l, double lp, double z)
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
double Kernel::BFKL_res_kc_v_r_full__Sub_Off(double l, double lp, double z)
{
    return Kernel::BFKL_res_kc_full__Sub_Off(l, lp, z);
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_kc_v_r_full__Sub_Diag(double l, double lp, double z)
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



////////////////////////////////////////////////
//BFKL_res_DGLAP 
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_DGLAP__Eps_Off(double l, double lp, double z)
{
    double res = BFKL__Eps_Off(l, lp, z);
    res += DGLAPoffOld(l, lp, z);

    if(z == 1) {
        res +=  DGLAPdiagOld(l, lp);
    }


    return res;
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_DGLAP__Eps_Diag(double l, double lp, double z)
{
    return BFKL__Eps_Diag(l, lp, z);
}

double Kernel::DGLAPdiagOldEps(double l, double lp)
{
    //putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = HarmonicEps(0., stepSize, 1., 3000);

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorOld(l, lp);
    return 2*res / stepSize; //for correction for the rapIntegration
}


double Kernel::DGLAPdiagOld(double l, double lp)
{
    putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = Harmonic(0., stepSize, 1., 3000);

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorOld(l, lp);
    return 2*res / stepSize; //for correction for the rapIntegration
}

double Kernel::DGLAPdiagNew(double l, double lp)
{
    putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = Harmonic(0., stepSize, 1., 3000);


    //Cros-check - start
    int il = -1, ilp = -1;
    for(int i = 0; i < SqrtExpXi.size(); ++i) {
       if(SqrtExpXi[i] == l)   il = i;
       if(SqrtExpXi[i] == lp) ilp = i;
    }
    if(il != -1 && ilp != -1) {
        myCheck(il, ilp) += 2;
    }
    else {
        cout << "Weired" << endl;
    }
    //Cros-check - end

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorNew(l, lp);
    return 2*res / stepSize; //for correction for the rapIntegration
}

double Kernel::DGLAPdiagNewEps(double l, double lp)
{
    //putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = HarmonicEps(0., stepSize, 1., 3000);

    /*
    //Cros-check - start
    int il = -1, ilp = -1;
    for(int i = 0; i < SqrtExpXi.size(); ++i) {
       if(SqrtExpXi[i] == l)   il = i;
       if(SqrtExpXi[i] == lp) ilp = i;
    }
    if(il != -1 && ilp != -1) {
        myCheck(il, ilp) += 2;
    }
    else {
        cout << "Weired" << endl;
    }
    */
    //Cros-check - end

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorNew(l, lp);
    return 2*res / stepSize; //for correction for the rapIntegration
}







double Kernel::DGLAPdiagExtra(double l, double lp, double z, double a)
{
    putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = Harmonic(0., stepSize, a, 3000);

    //Cros-check
    int il = -1, ilp = -1;
    for(int i = 0; i < SqrtExpXi.size(); ++i) {
       if(SqrtExpXi[i] == l)   il = i;
       if(SqrtExpXi[i] == lp) ilp = i;
    }
    if(il != -1 && ilp != -1) {
        ++myCheck(il, ilp);
    }
    else {
        cout << "Weired" << endl;
    }

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorExtra(l, lp, z, a);
    return 2*res / stepSize; //for correction for the rapIntegration
}

double Kernel::DGLAPdiagExtraEps(double l, double lp, double z, double a)
{
    //putZero = true;
    const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
    static double harmSum = HarmonicEps(0., stepSize, a, 3000);

    /*
    //Cros-check
    int il = -1, ilp = -1;
    for(int i = 0; i < SqrtExpXi.size(); ++i) {
       if(SqrtExpXi[i] == l)   il = i;
       if(SqrtExpXi[i] == lp) ilp = i;
    }
    if(il != -1 && ilp != -1) {
        ++myCheck(il, ilp);
    }
    else {
        cout << "Weired" << endl;
    }
    */

    const int nf = 4;
    double sum = harmSum + (33 - 2*nf) / 36.0;
    double res = sum *  DGLAPfactorExtra(l, lp, z, a);
    return /*2**/res / stepSize; //for correction for the rapIntegration
}










//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_DGLAP__ZEps_Off(double l, double lp, double z)
{
    double res = BFKL__Eps_Off(l, lp, z);
    res += DGLAPoffOldEps(l, lp, z);

    if(z == 1) {
        res +=  DGLAPdiagOldEps(l, lp);
    }


    return res;
}

//Diagonal kernel with F(k)
double Kernel::BFKL_res_DGLAP__ZEps_Diag(double l, double lp, double z)
{
    return BFKL__Eps_Diag(l, lp, z);
}






////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP
////////////////////////////////////////////////


double Kernel::BFKL_res_kc_full_DGLAP__Eps_Off(double l, double lp, double z)
{
    //Adding normal BFKL
    double res =  BFKL_res_kc_full__Eps_Off(l, lp, z);
    res += DGLAPoffOld(l, lp, z);

    if(z == 1) {
        res += DGLAPdiagOld(l, lp);
    }


    return res;
}

double Kernel::BFKL_res_kc_full_DGLAP__Eps_Diag(double l, double lp, double z)
{
    //return 0;
    return BFKL__Eps_Diag(l, lp, z);
}




////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP_*_kc
////////////////////////////////////////////////


double Kernel::DGLAPextraCommon(double l, double lp, double z, double a)
{
    double res =  Kernel::BFKL_res_kc_v_r_full__Eps_Off(l, lp, z); //Original term

    if(lp <= l) {
        res += DGLAPoffNew(l, lp, z); //simple DGLAP

        if(z == 1) {
            res +=  DGLAPdiagNew(l, lp); 
        }
    }

    //Including of the last term
    if(lp > l && 0) {
        //double a = lp*lp / (l*l);
        double diff = getRapIndex(z) - getRapIndexReal(1./a);

        if(diff > eps) {
            res += DGLAPoffExtra(l, lp, z, a); //with real singularity

            if(diff < 1 + eps) {
                //cout << l << " "< lp << endl;
                res +=  DGLAPdiagExtra(l, lp, z, a); 
            }
        }
    }
    return res;
}


double Kernel::DGLAPextraCommonEps(double l, double lp, double z, double a)
{
    double res =  Kernel::BFKL_res_kc_v_r_full__Eps_Off(l, lp, z); //Original term

    if(lp <= l) {
        res += DGLAPoffNewEps(l, lp, z); //simple DGLAP

        if(z == 1) {
            res +=  DGLAPdiagNewEps(l, lp); 
        }
    }

    //Including of the last term
    if(lp > l && 1 &&  z < 0.1) {
        int idRap = getRapIndex(z);
        //double a = lp*lp / (l*l);
        double diff = getRapIndex(z) - getRapIndexReal(1./a);

        //if(diff > eps) {
        res += DGLAPoffExtraEps(l, lp, z, a); //with real singularity
        //}

        

        //Find the maximum the maximum
        double valMax = -1e20;
        int zMax = -1;
        const static double stepSize = (rapMax - rapMin) / (Nrap - 1);
        for(int i = 0; i <= Nrap; ++i) {
            double rapNow = rapMin + i*stepSize;
            double zNow = exp(-rapNow);
            double valNow = DGLAPoffExtraEps(l, lp, zNow, a);
            if(valNow > valMax) {
                zMax = zNow;
                valMax = valNow;
            }
        }
        if(z == zMax) {
            double factor = (idRap == 0 || idRap == Nrap) ? 2 : 1;
            res +=   factor * DGLAPdiagExtraEps(l, lp, z, a); 
        }


        /*
        if(diff > -1 + eps && diff < 1 - eps) {
            double f = 1 - abs(diff);
            double factor = (idRap == 0 || idRap == Nrap) ? 2 : 1;
            res +=  f * factor * DGLAPdiagExtraEps(l, lp, z, a); 
        }
        */

        /*
        if(diff > eps && diff < 1 + eps) {
            //cout << l << " "< lp << endl;
            double factor = (idRap == 0 || idRap == Nrap) ? 2 : 1;
            res +=  factor * DGLAPdiagExtraEps(l, lp, z, a); 
        }
        */
    }
    return res;
}



////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP_simp_kc
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full_DGLAP_simp_kc__Eps_Off(double l, double lp, double z)
{
    double a = lp*lp / (l*l);
    //cout << "Approximate constraint" << endl;
    return DGLAPextraCommon(l, lp, z, a);
}

double Kernel::BFKL_res_kc_full_DGLAP_simp_kc__Eps_Diag(double l, double lp, double z)
{
    //return 0;
    return Kernel::BFKL_res_kc_v_r_full__Eps_Diag(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full_DGLAP_simp_kc__ZEps_Off(double l, double lp, double z)
{
    double a = lp*lp / (l*l);
    //cout << "Approximate constraint" << endl;
    return DGLAPextraCommonEps(l, lp, z, a);
}

double Kernel::BFKL_res_kc_full_DGLAP_simp_kc__ZEps_Diag(double l, double lp, double z)
{
    //return 0;
    return Kernel::BFKL_res_kc_v_r_full__Eps_Diag(l, lp, z);
}




////////////////////////////////////////////////
//BFKL_res_kc_full_DGLAP_full_kc
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Kernel::BFKL_res_kc_full_DGLAP_full_kc__Eps_Off(double l, double lp, double z)
{
    double a = lp*lp / (l*l) + 1;
    //cout << "Full constraint" << endl;
    return DGLAPextraCommon(l, lp, z, a);
}

double Kernel::BFKL_res_kc_full_DGLAP_full_kc__Eps_Diag(double l, double lp, double z)
{
    //return 0;
    return Kernel::BFKL_res_kc_v_r_full__Eps_Diag(l, lp, z);
}

