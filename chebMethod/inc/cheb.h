#ifndef _CHEB_
#define _CHEB_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <cassert>
#include <functional>
#include <armadillo>
#include "nodes.h"

typedef double Double;

Double Integral61( std::function<Double(Double)> fun, Double xmin, Double xmax, Double &err );




//If defined between 0 and 1
inline double GetChebIntegral(int i, int j)
{
    auto Int = [](int k) {
        return (k%2==0) ? 1./(1.-k*k) : 0;
    };
    return (Int(i+j) + Int(i-j))/2;
}

inline arma::vec getChebCoef(int nCoef, std::function<double(double)> fun)
{
    arma::mat trMat = GetCoefs(nCoef); //zero coef must be divided by two
    arma::vec xi    = GetNodes(nCoef);
    arma::vec vals  = xi;
    for(int i = 0; i < nCoef; ++i)
        vals(i) = fun(xi(i));
    arma::vec res = trMat * vals;
    res(0)            /= 2;
    res(res.n_rows-1) /= 2;
    return res;
}


//Cheb Pol index for y and K
struct Indx {
    int iy, iK;
    Indx() {}
    Indx(int iy_, int yK_) : iy(iy_), iK(yK_) {}
};

inline Indx toLocal(int Ncheb, int I)
{
    return Indx(I/Ncheb, I%Ncheb);
}



double Cheb(int n, double x);
void IntegrateMC();
void IntegrateDummy(Indx idp, Indx idq);
double IntegrateClever(Indx idp, Indx idq);

const double rapMin = 0;
const double rapMax = -log(1e-6);

const double Kmin = log(1e-2);
const double Kmax = log(1e6);


inline double Rand(double a, double b) {return a + (b-a)*rand()/(RAND_MAX+0.);}


struct Point {
    double y;
    double K;
    Point(double y_, double K_) : y(y_), K(K_)
    {
        assert(rapMin <= y);
        assert(y <= rapMax);
        assert(Kmin <= K);
        if(K < Kmax + 1e-8)
            K = std::min(Kmax, K);
        assert(K <= Kmax);
    }
    Point() {}
    friend bool operator<(const Point& a, const Point &b);
    void rInit() {
        y  = Rand(rapMin,rapMax);
        K  = Rand(Kmin,Kmax);
    }
};

inline bool operator<(const Point& a, const Point &b)
{
    if(a.y < b.y) return true;
    if(a.y > b.y) return false;

    if(a.K < b.K) return true;
    else          return false;
}

double F(Indx indx, Point p);
double Fy(int n, double y);
double FK(int n, double K);




struct IntPoint {
    Point p,q;
    double ker;
    IntPoint(Point p_, Point q_, double ker_)
                : p(p_), q(q_), ker(ker_) {}
};

std::vector<IntPoint> KernelGeneral(Point p, Point q, double wR, double wS);

class Integrator {
    int N;
    std::vector<IntPoint> points;
    std::vector<int> breaks;
    std::map<double,double> pY, pK, qY, qK;

public:
    Integrator(int N_); 
    double Eval(Indx idp, Indx idq);
};



inline double abs(Point p) {return sqrt(p.y*p.y + p.K*p.K);}
inline Point operator+(Point a, Point b) { return Point(a.y+b.y, a.K+b.K);}
inline Point operator-(Point a, Point b) { return Point(a.y-b.y, a.K-b.K);}
inline Point operator*(double c, Point a) { return Point(c*a.y, c*a.K);}





#endif
