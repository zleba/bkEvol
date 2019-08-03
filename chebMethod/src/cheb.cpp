#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
//#include <boost/math/special_functions/chebyshev.hpp>

#include "cheb.h"
#include "nodes.h"
#include <vector>
#include "kernels.h"

double eps = 1e-10;

Kernel *kerObj;


/*
double Cheb(int n, double x)
{
    if(x < -1 || x > 1) return 0;
    if(n == 0) return 1;
    else if(n == 1) return x;
    else return 2*x*Cheb(n-1,x) - Cheb(n-2,x);
}
*/

//Fast evaluation of Cheb (between 0 and 1)
double Cheb(int n, double x)
{
    if(x < 0 || x > 1)
        cout <<"Hela " <<  x << endl;
    assert(-eps <= x && x <= 1+eps);
    x = max( 0., x);
    x = min( 1., x);

    //x = 2*x-1;
    //assert(-1 <= x && x <= 1);
    //if(x < -1 || x > 1) return 0;
    if(n == 0) return 1;
    if(n == 1) return 2*x-1;

    int n2 = n/2;
    if(n % 2 == 0)  //even
        return 2*pow(Cheb(n2,x),2) - 1;
    else //odd
        return 2*Cheb(n2+1,x)*Cheb(n2,x) - (2*x-1);
}



using namespace std;

/*
double F(double x, double l2)
{
    double y = -log(x);
    return y;
}
*/

double FK(int n, double K)
{
    assert(Kmin <= K);
    assert(K <= Kmax);
    double KRen = (K-Kmin)/(Kmax-Kmin);
    assert(abs(KRen) <= 1 + eps);
    KRen = max( 0., KRen);
    KRen = min( 1., KRen);
    return Cheb(n, KRen);
}

double Fy(int n, double y)
{
    assert(rapMin <= y);
    assert(y <= rapMax);
    double yRen = (y - rapMin) / (rapMax - rapMin);
    assert(yRen <= 1);
    assert(yRen >= 0);
    return Cheb(n, yRen);
}


double F(Indx indx, Point p)
{
    return Fy(indx.iy, p.y) * FK(indx.iK, p.K);
    //double yRen = 2.*(p.y - (rapMin+rapMax)/2.) / (rapMax - rapMin);
    //double KRen = 2.*(p.K - (Kmin+Kmax)/2.) / (Kmax - Kmin);
    //double res = Cheb(indx.iy, yRen) * Cheb(indx.iK, KRen);
    //return res;
}

double KernelAdv(Indx idp, Indx idq, Point p, Point q) //Extended DGLAP term
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double k2 = exp(p.K);
    double l2 = exp(q.K);
    
    double zm = k2 / l2;
    double z = x/xp;

    if(z > 1) return 0;
    if(l2 < k2) return 0;
    if(z >= zm) return 0;
    
    q.y -= -log(zm);
    Point qq;
    qq.y = p.y + log(zm);
    qq.K = q.K;

    assert(z < zm);
    assert(x/zm < 1);
    double Fun = z/(zm - z) * (F(idq, q) - F(idq, qq));
    Fun -= -log(1-x/zm) * F(idq, qq);
    Fun *= F(idp, p);
    //Fun = 1;

    return Fun;
}

//Clasical DGLAP (regular + singular)
pair<double,double> KernelBasic(Indx idp, Indx idq, Point p, Point q)
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double k2 = exp(p.K);
    double l2 = exp(q.K);
    
    double z = x/xp;

    if(l2 > k2) return {0.,0.};
    if(z >= 1) return {0.,0.};
    
    Point qq;
    qq.y = p.y;
    qq.K = q.K;

    double Fl = F(idp, p);
    double Fr = F(idq, q);
    double Frsub = F(idq, qq);

    double FunReg  = Fl * z/(1 - z) * (Fr - Frsub);
    double FunSing = -log(1-x) * Fl*Frsub;
    assert(isfinite(FunSing));

    return {FunReg, FunSing};
}





vector<IntPoint> KernelGeneral(Point p, Point q, double wR, double wS)
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double l = exp(0.5*p.K);
    double lp = exp(0.5*q.K);
    
    double z = x/xp;

    //if(l > lp) return {};
    if(z > 1) return {};
    
    /*
    Point qq;
    qq.y = p.y;
    qq.K = q.K;
    */

    //double Fl = F(idp, p);
    //double Fr = F(idq, q);
    //double Frsub = F(idq, qq);

    //double FunReg = Fl * z/(1 - z) * (Fr - Frsub);
    //double FunSin = -log(1-x) * Fl*Frsub;

    //double kerReg = z/(1 - z);
    //double kerSin = -log(1-x);

    Point qq;
    qq.y = q.y;
    qq.K = p.K;

    double kerReg = kerObj->BFKL__Eps_Off(l, lp, z);
    double kerSin = kerObj->BFKL__Eps_Diag(l, lp, z);

    if(!isfinite(kerReg)) {
        cout <<"Helenka "  << l << " "<< lp <<" "<<  z << endl;
        cout << kerReg << endl;
        cout << kerSin << endl;
        assert(0);
    }
    assert(isfinite(kerReg));
    assert(isfinite(kerSin));

    return { IntPoint(p,q,  kerReg*wR),
             IntPoint(p,qq, kerSin*wR) };


    //return { IntPoint(p,q,  kerReg*wR),
             //IntPoint(p,qq,-kerReg*wR+kerSin*wS) };
}

pair<double,double> KernelSpecific(Indx idp, Indx idq, Point p, Point q)
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double l = exp(0.5*p.K);
    double lp = exp(0.5*q.K);
    
    double z = x/xp;

    assert(z <= 1+eps); 
    
    /*
    Point qq;
    qq.y = p.y;
    qq.K = q.K;
    */

    //double Fl = F(idp, p);
    //double Fr = F(idq, q);
    //double Frsub = F(idq, qq);

    //double FunReg = Fl * z/(1 - z) * (Fr - Frsub);
    //double FunSin = -log(1-x) * Fl*Frsub;

    //double kerReg = z/(1 - z);
    //double kerSin = -log(1-x);

    Point qq;
    qq.y = q.y;
    qq.K = p.K;


    double Fl = F(idp, p);
    double Fr = F(idq, q);
    double Frsub = F(idq, qq);

    double kerReg = kerObj->BFKL__Eps_Off(l, lp, z);
    double kerSin = kerObj->BFKL__Eps_Diag(l, lp, z);

    return {Fl*(kerReg*Fr + kerSin*Frsub), 0};
}








void IntegrateMC()
{
    Indx idp(2,1);
    Indx idq(2,1);

    const int N = 15000000;
    double sum  = 0;
    double sum2 = 0;
    for(int i = 0; i < N; ++i) {
        Point p, q;
        p.rInit();
        q.rInit();

        double res = KernelBasic(idp, idq, p, q).first;
        sum  += res;
        sum2 += res*res;
    }

    double val = sum / N;
    double err = sqrt(sum2/N - val*val)/sqrt(N);
    val *= pow((rapMax-rapMin)*(Kmax-Kmin),2);
    err *= pow((rapMax-rapMin)*(Kmax-Kmin),2);
    cout << val << " "<< err << endl;

}

#include "nodes.h"

//Simple method of integration (Trapezius rule)
void IntegrateDummy(Indx idp, Indx idq)
{

    double sum = 0;
    const int N = 80;
    for(int ip = 0; ip <= N; ++ip)
    for(int jp = 0; jp <= N; ++jp)
    for(int iq = 0; iq <= N; ++iq)
    for(int jq = 0; jq <= N; ++jq) {
        Point p(rapMin + (rapMax-rapMin)/N*ip, Kmin + (Kmax-Kmin)/N*jp);
        Point q(rapMin + (rapMax-rapMin)/N*iq, Kmin + (Kmax-Kmin)/N*jq);
        double res = KernelBasic(idp, idq, p, q).first;

        double fact = 1;
        fact *= (ip == 0 || ip == N) ? 0.5 : 1;
        fact *= (jp == 0 || jp == N) ? 0.5 : 1;
        fact *= (iq == 0 || iq == N) ? 0.5 : 1;
        fact *= (jq == 0 || jq == N) ? 0.5 : 1;
        sum += res * fact;
    }

    sum *= pow((rapMax-rapMin)/N*(Kmax-Kmin)/N, 2);
    cout << sum << endl;
}



double IntegrateClever(Indx idp, Indx idq)
{

    const int N = 60;
    //Nodes nodRap(N, rapMin,rapMax);
    //Nodes nodK(N, Kmin, Kmax);

    //for points between 0 and 1
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    //cout << sum(wgts) << endl;


    double sum = 0;
    for(int ip = 0; ip <= N; ++ip) 
    for(int jp = 0; jp <= N; ++jp) {
        Point p(rapMin + xi(ip)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(jp));
        double wBase = (rapMax-rapMin)*(Kmax-Kmin) * wgts(ip)*wgts(jp);

        //Skip integration over nothing
        if(ip == 0 || jp == 0) continue;

        for(int jq = 0; jq <= N; ++jq) {
            for(int iq = 0; iq <= N; ++iq) { //Integral over z
                Point q(rapMin + xi(iq)*(p.y-rapMin), Kmin + (p.K-Kmin)*xi(jq));
                double wAdd = (p.y-rapMin)*(p.K-Kmin) * wgts(iq)*wgts(jq);
                //double kerReg, kerSin;
                //tie(kerReg, kerSin) = KernelBasic(idp, idq, p, q);
                //double ker = kerReg + kerSin/(p.y-rapMin);
                
                sum += KernelSpecific(idp,idq,p,q).first * wBase * wAdd;
                //sum += ker* wBase * wAdd;
            }
        }
    }

    //uY*uK * ker * vY*vK

    //cout <<"Radek " <<  sum << endl;
    return sum;

}


Integrator::Integrator(int N_) : N(N_) 
{
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    kerObj = new Kernel();
    kerObj->mu2 = exp(Kmin);
    kerObj->eps = 1e-7;
    kerObj->LnFreeze2 = log(1);



    for(int ip = 0; ip <= N; ++ip) 
    for(int jp = 0; jp <= N; ++jp) {
        Point p(rapMin + xi(ip)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(jp));
        double wBase = (rapMax-rapMin)*(Kmax-Kmin) * wgts(ip)*wgts(jp);

        //Skip integration over nothing
        if(ip == 0 || jp == 0) continue;

        pY[p.y] = 0;
        pK[p.K] = 0;

        for(int jq = 0; jq <= N; ++jq) {
            for(int iq = 0; iq <= N; ++iq) { //Integral over z
                Point q(rapMin + xi(iq)*(p.y-rapMin), Kmin + (p.K-Kmin)*xi(jq));
                double wAdd = (p.y-rapMin)*(p.K-Kmin) * wgts(iq)*wgts(jq);
                double wgt = wBase * wAdd;

                auto res = KernelGeneral(p, q, wgt, wgt/(p.y-rapMin));
                for(auto &v : res)
                    points.push_back(v);

                qY[q.y] = 0;
                qK[q.K] = 0;
            }

        }
    }

    cout << "Sizes " << pK.size() <<" "<< pY.size() << " : "<< qK.size() <<" "<< qY.size() << endl;

    sort(points.begin(),points.end(), [](IntPoint &a, IntPoint &b) {
        if(a.p.y > b.p.y) return true;
        if(a.p.y < b.p.y) return false;
        if(a.p.K > b.p.K) return true;
        if(a.p.K < b.p.K) return false;

        if(a.q.y > b.q.y) return true;
        if(a.q.y < b.q.y) return false;
        if(a.q.K > b.q.K) return true;
        if(a.q.K < b.q.K) return false;
        return false;
    });

    double org = points[0].q.y;
    int c = 0;

    for(auto v : points) {
        if(v.q.y != org) {
            breaks.push_back(c);
            org = v.q.y;
            c = 0;
        }
        ++c;
    }
    breaks.push_back(c);
    cout << accumulate(breaks.begin(),breaks.end(),0)  << endl;
    cout << points.size()  << endl;
    assert(accumulate(breaks.begin(),breaks.end(),0) == points.size());


}

double Integrator::Eval(Indx idp, Indx idq)
{
    for(auto &m : pK) m.second = FK(idp.iK, m.first);
    for(auto &m : pY) m.second = Fy(idp.iy, m.first);

    for(auto &m : qK) m.second = FK(idq.iK, m.first);
    for(auto &m : qY) m.second = Fy(idq.iy, m.first);

    double sum = 0;


    for(auto &v : points) {
        //sum += v.ker * F(idp, v.p) * F(idq, v.q);
        sum += v.ker *  pY[v.p.y]*pK[v.p.K]* qY[v.q.y]*qK[v.q.K];
    }

    /*
    auto it = points.begin();
    for(auto b : breaks) {
        double fac = pY[it->p.y]*pK[it->p.K]* qY[it->q.y];
        double sumTemp = 0;
        for(int k = 0; k < b; ++k) {
            sum += it->ker * qK[it->q.K];
            ++it;
            if(it == points.end()) break;
        }
        sum += fac * sumTemp;
    }
    */

    return sum;

    //uY*uK * ker * vY*vK

}

