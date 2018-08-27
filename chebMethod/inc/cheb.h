#ifndef _CHEB_
#define _CHEB_

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

double Cheb(int n, double x);
void IntegrateMC();
void IntegrateDummy();
void IntegrateClever();

const double rapMin = 0;
const double rapMax = -log(1e-6);

const double Kmin = log(1e-2);
const double Kmax = log(1e6);


inline double Rand(double a, double b) {return a + (b-a)*rand()/(RAND_MAX+0.);}


struct Point {
    double y;
    double K;
    Point(double y_, double K_) : y(y_), K(K_) {}
    Point() {}
    void rInit() {
        y  = Rand(rapMin,rapMax);
        K  = Rand(Kmin,Kmax);
    }
};


struct IntPoint {
    Point p,q;
    double ker;
    IntPoint(Point p_, Point q_, double ker_)
                : p(p_), q(q_), ker(ker_) {}
};

//Cheb Pol index for y and K
struct Indx {
    int iy, iK;
    Indx() {}
    Indx(int iy_, int yK_) : iy(iy_), iK(yK_) {}
};

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
