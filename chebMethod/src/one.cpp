#include <iostream>
#include "cheb.h"
#include "kernels.h"
#include "nodes.h"
#include "TH2D.h"

using namespace std;

//File to fill the kernel in the new way
//
//


struct IntP {
    Point q;
    double ker;
    IntP(Point q_, double ker_)
                : q(q_), ker(ker_) {}
};




struct Filler {
    Kernel *kerObj;
    void Init();
    arma::mat EvalAll(string kerName,  int Ncheb, long long nSplit, long long nNow);
    std::map<Point,double> GetIntegral(Point p);
    std::map<Point,double> GetIntegralFull(Point p);

    pair<double,double> KernelSpecific1D(Indx idq, Point p, Point q);
    map<Point,double> KernelGeneral1D(Point p, Point q, double wR, double wS);
    double IntegrateClever(int N, Point p, Indx idq);
    double IntegrateClever2(int N, Point p, Indx idq);

    double IntegrateCleverDGLAPadd(int N, Point p, Indx idq);
    std::map<Point,double>  GetIntegralDGLAPadd(Point p);

    double IntegrateCleverDGLAP(int N, Point p, Indx idq);
    std::map<Point,double>  GetIntegralDGLAP(Point p);
};





std::map<Point,double> Filler::GetIntegral(Point p)
{
    const int N = 256;
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    std::map<Point,double> points;
    for(int iq = 0; iq <= N; ++iq) { //Integral over z
        for(int jq = 0; jq <= N; ++jq) {
            Point q(rapMin + xi(iq)*(p.y-rapMin), Kmin + (Kmax-Kmin)*xi(jq));
            double wgt = (p.y-rapMin)*(Kmax-Kmin) * wgts(iq)*wgts(jq);

            auto res = KernelGeneral1D(p, q, wgt, wgt/(p.y-rapMin));
            for(auto &v : res)
                points[v.first] += v.second;
        }
    }
    return points;
}

//More advanced integration method
std::map<Point,double> Filler::GetIntegralFull(Point p) 
{
    const int N = 256;
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double k2 = exp(p.K);
    std::map<Point,double> points;
    for(int iq = 0; iq <= N; ++iq) { //Integral over z
        double yp = rapMin + xi(iq)*(p.y-rapMin);
        double z = exp(-yp);
        double q2 = (1/z - 1) * k2;


        //Vector with all breaking points
        vector<double> bp;
        bp.push_back(Kmin);
        bp.push_back(Kmax);
        double lMin = 2*log(sqrt(k2) - sqrt(q2));
        double lMax = 2*log(sqrt(k2) + sqrt(q2));
        if(Kmin < lMin && lMin < Kmax) {
            bp.push_back(lMin);
        }
        if(Kmin < lMax && lMax < Kmax) {
            bp.push_back(lMax);
        }
        if(log(4*k2) < Kmax) {
            bp.push_back(log(4*k2));
        }
        sort(bp.begin(), bp.end());

        for(int ib = 0; ib < bp.size()-1; ++ib) {
            double KminNow = bp[ib];
            double KmaxNow = bp[ib+1];
            for(int jq = 0; jq <= N; ++jq) {
                Point q(yp, KminNow + (KmaxNow-KminNow)*xi(jq));
                double wgt = (p.y-rapMin)*(KmaxNow-KminNow) * wgts(iq)*wgts(jq);

                auto res = KernelGeneral1D(p, q, wgt, wgt/(p.y-rapMin));
                for(auto &v : res)
                    points[v.first] += v.second;
            }
        }
    }
    return points;
}





double Filler::IntegrateClever(int N, Point p, Indx idq)
{
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double sum = 0;
    for(int iq = 0; iq <= N; ++iq) { //Integral over z
        for(int jq = 0; jq <= N; ++jq) {
            Point q(rapMin + xi(iq)*(p.y-rapMin), Kmin + (Kmax-Kmin)*xi(jq));
            double wgt = (p.y-rapMin)*(Kmax-Kmin) * wgts(iq)*wgts(jq);
            sum += KernelSpecific1D(idq,p,q).first * wgt;
        }
    }
    return sum;
}


double Preg(double z) { return   z*z*(1-z) - (2*z+1);}
double Psin(double z) { return   1/(1-z); }



double Filler::IntegrateCleverDGLAPadd(int N, Point p, Indx idq)
{
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double x  = exp(-p.y);
    double k2 = exp(p.K);
    if(x >= 0.5) return 0;

    const double as = 0.2;
    double sum = 0;
    for(int jq = 0; jq <= N; ++jq) { //Over kT
        double KmaxNow   = min(Kmax,log(k2*(1/x-1)));
        double KNow   = p.K + (KmaxNow-p.K)*xi(jq);
        double l2 = exp(KNow);

        double xp = x*(l2/k2+1);
        assert(xp < 1+kerObj->eps);
        xp = min(1.,xp);
        Point q2(-log(xp),log(l2));

        double wgtL = (KmaxNow-p.K)*wgts(jq);
        assert(KmaxNow-p.K>=0);
        if(wgtL == 0) continue;
        for(int iq = 0; iq <= N; ++iq) { //Integral over z
            double rapMaxNow = min(rapMax, max(rapMin,-log(xp)));
            double rapNow = rapMin + xi(iq)*(rapMaxNow-rapMin);
            Point q(rapNow, KNow);

            double wgtY = (rapMaxNow-rapMin)*wgts(iq);
            double z = exp(-q.y);
            if(xp/z > 1-kerObj->eps) continue; //TODO replace skip
            sum +=  wgtY*wgtL * as * ( Preg(xp/z)*F(idq, q)  + Psin(xp/z)*(F(idq,q) - F(idq,q2)) );
        }

        if(xp < 1-kerObj->eps) {
            const int nf = 4;
            double Fact =  (33 - 2*nf) / 36.0;
            sum += wgtL * as * F(idq,q2)* (-log(1/xp) + Fact) ;
        }

    }
    return sum;
}


std::map<Point,double>  Filler::GetIntegralDGLAPadd(Point p)
{
    std::map<Point,double> points;
    const int N = 256;
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double x  = exp(-p.y);
    double k2 = exp(p.K);
    if(x >= 0.5) return points;

    const double as = 0.2;
    for(int jq = 0; jq <= N; ++jq) { //Over kT
        double KmaxNow   = min(Kmax,log(k2*(1/x-1)));
        double KNow   = p.K + (KmaxNow-p.K)*xi(jq);
        double l2 = exp(KNow);

        double xp = x*(l2/k2+1);
        assert(xp < 1+kerObj->eps);
        xp = min(1.,xp);
        double yP = -log(xp);
        Point q2(yP,log(l2));

        double wgtL = (KmaxNow-p.K)*wgts(jq);
        assert(KmaxNow-p.K>=0);
        if(wgtL == 0) continue;
        for(int iq = 0; iq <= N; ++iq) { //Integral over z
            double rapMaxNow = min(rapMax, max(rapMin,yP));
            double rapNow = rapMin + xi(iq)*(rapMaxNow-rapMin);
            Point q(rapNow, KNow);

            double wgtY = (rapMaxNow-rapMin)*wgts(iq);
            double z = exp(-q.y);
            if(xp/z > 1-kerObj->eps) continue; //TODO replace skip
            //sum +=  wgtY*wgtL * as * ( Preg(xp/z)*F(idq, q)  + Psin(xp/z)*(F(idq,q) - F(idq,q2)) );
            points[q] +=  wgtY*wgtL * as * (Preg(xp/z) + Psin(xp/z)); //off diagonal
            points[q2]+= -wgtY*wgtL * as * Psin(xp/z); //off diagonal
        }

        if(xp < 1-kerObj->eps) {
            const int nf = 4;
            double Fact =  (33 - 2*nf) / 36.0;
            //sum += wgtL * as * F(idq,q2)* (-log(1/xp) + Fact);
            points[q2]+= wgtL * as * (-log(1/xp) + Fact);
        }
    }
    return points;
}

double Filler::IntegrateCleverDGLAP(int N, Point p, Indx idq)
{
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double x  = exp(-p.y);
    double k2 = exp(p.K);
    if(x >= 1-kerObj->eps) return 0;

    const double as = 0.2;
    double sum = 0;
    for(int jq = 0; jq <= N; ++jq) { //Over kT
        double KNow   = Kmin + (p.K - Kmin)*xi(jq);
        double l2 = exp(KNow);

        Point q2(-log(x),log(l2));

        double wgtL = (p.K - Kmin)*wgts(jq);
        assert(p.K - Kmin >= 0);
        if(wgtL == 0) continue;
        for(int iq = 0; iq <= N; ++iq) { //Integral over z
            double rapMaxNow = min(rapMax, max(rapMin,-log(x)));
            double rapNow = rapMin + xi(iq)*(rapMaxNow-rapMin);
            Point q(rapNow, KNow);

            double wgtY = (rapMaxNow-rapMin)*wgts(iq);
            double z = exp(-q.y);
            if(x/z > 1-kerObj->eps) continue; //TODO replace skip
            sum +=  wgtY*wgtL * as * l2/k2* ( Preg(x/z)*F(idq, q)  + Psin(x/z)*(F(idq,q) - F(idq,q2)) );
        }

        if(x < 1-kerObj->eps) {
            const int nf = 4;
            double Fact =  (33 - 2*nf) / 36.0;
            sum += wgtL * as * l2/k2* F(idq,q2)* (-log(1/x) + Fact);
        }
    }
    return sum;
}

std::map<Point,double>  Filler::GetIntegralDGLAP(Point p)
{
    std::map<Point,double> points;
    const int N = 256;
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double x  = exp(-p.y);
    double k2 = exp(p.K);
    if(x >= 1-kerObj->eps) return points;

    const double as = 0.2;
    double sum = 0;
    for(int jq = 0; jq <= N; ++jq) { //Over kT
        double KNow   = Kmin + (p.K - Kmin)*xi(jq);
        double l2 = exp(KNow);

        Point q2(-log(x),log(l2));

        double wgtL = (p.K - Kmin)*wgts(jq);
        assert(p.K - Kmin >= 0);
        if(wgtL == 0) continue;
        for(int iq = 0; iq <= N; ++iq) { //Integral over z
            double rapMaxNow = min(rapMax, max(rapMin,-log(x)));
            double rapNow = rapMin + xi(iq)*(rapMaxNow-rapMin);
            Point q(rapNow, KNow);

            double wgtY = (rapMaxNow-rapMin)*wgts(iq);
            double z = exp(-q.y);
            if(x/z > 1-kerObj->eps) continue; //TODO replace skip
            points[q]  += wgtY*wgtL * as * l2/k2* (Preg(x/z) + Psin(x/z));
            points[q2] += wgtY*wgtL * as * l2/k2* (- Psin(x/z));
        }

        if(x < 1-kerObj->eps) {
            const int nf = 4;
            double Fact =  (33 - 2*nf) / 36.0;
            points[q2] += wgtL * as *l2/k2* (-log(1/x) + Fact);
        }
    }
    return points;
}














double Filler::IntegrateClever2(int N, Point p, Indx idq)
{
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);  //basic nodes

    double sum = 0;
    for(int iq = 0; iq <= N; ++iq) { //Integral over z
        double wY = (p.y-rapMin)*wgts(iq);

        double Kb = p.K + log(4);
        Kb = min(Kmax, Kb); 

        for(int jq = 0; jq <= N; ++jq) {
            Point q(rapMin + xi(iq)*(p.y-rapMin), Kmin + (Kb-Kmin)*xi(jq));
            double wgt = wY*(Kb-Kmin) *wgts(jq);
            sum += KernelSpecific1D(idq,p,q).first * wgt;
        }
        if(Kb < Kmax) {
            for(int jq = 0; jq <= N; ++jq) {
                Point q(rapMin + xi(iq)*(p.y-rapMin), Kb + (Kmax-Kb)*xi(jq));
                double wgt = wY*(Kmax-Kb) *wgts(jq);
                sum += KernelSpecific1D(idq,p,q).first * wgt;
            }
        }

    }
    return sum;
}





arma::mat Filler::EvalAll(string kerName,  int Ncheb, long long nSplit, long long nNow)
{
    assert(Ncheb % 2 == 1);

    arma::vec wgts = GetWeights(Ncheb); //basic weights
    arma::vec xi   = GetNodes(Ncheb);  //basic nodes

     std::map<Point,double> (Filler::*fun) (Point);
     if(kerName == "simpleBFKL")
         fun = &Filler::GetIntegral;
     else if(kerName == "fullBFKL")
         fun = &Filler::GetIntegralFull;
     else if(kerName == "DGLAP")
         fun = &Filler::GetIntegralDGLAP;
     else if(kerName == "DGLAPadd")
         fun = &Filler::GetIntegralDGLAPadd;
     else
         assert(0);


    std::map<double,double> qY, qK;
    map<Point,double> points[100][100];
    long long nItem  = Ncheb * Ncheb;
    auto iStart = (nItem*nNow)/nSplit;
    auto iEnd   = (nItem*(nNow+1))/nSplit;

    #pragma omp parallel for
    for(int iy = 0; iy < Ncheb; ++iy) { //Loop over cheb nodes
        for(int ik = 0; ik < Ncheb; ++ik) {
            Point p(rapMin + xi(iy)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(ik));
            int iGlob = iy*Ncheb + ik;
            if(! (iStart <= iGlob && iGlob < iEnd)) continue;

            points[iy][ik] = (this->*fun)(p);
            //points[iy][ik] = GetIntegralDGLAP(p);
        }
        //cout << "Node " << iy <<"/" << Ncheb <<" evaluated" <<  endl;
    }

    for(int iy = 0; iy < Ncheb; ++iy)  //Loop over cheb nodes
    for(int ik = 0; ik < Ncheb; ++ik) {
        for(auto &point : points[iy][ik]) {
            qY[point.first.y] = 0;
            qK[point.first.K] = 0;
        }
    }


    cout << "Kernels evaluated" << endl;




    arma::mat convMat(Ncheb*Ncheb, Ncheb*Ncheb);



    for(int jy = 0; jy < Ncheb; ++jy) //Loop over cheb Polynomials
    for(int jk = 0; jk < Ncheb; ++jk) {
        int jGlob = jy*Ncheb + jk;

        //if(!(jy == 2 && jk == 2)) continue; //TODO cut

        for(auto &m : qK) m.second = FK(jk, m.first);
        for(auto &m : qY) m.second = Fy(jy, m.first);

        for(int iy = 0; iy < Ncheb; ++iy) //Loop over cheb nodes
        for(int ik = 0; ik < Ncheb; ++ik) {
            //if(!(iy == 1 && ik == 1)) continue; //TODO cut
            int iGlob = iy*Ncheb + ik;
            if(! (iStart <= iGlob && iGlob < iEnd)) continue;
            double sum = 0;
            for(auto &v : points[iy][ik])
                //sum += v.ker * F(idq, v.q);
                sum += v.second *  qY[v.first.y]*qK[v.first.K];
            convMat(iGlob,jGlob) = sum;

            //cout <<"Sum is " <<  sum << endl;
            /*
            Point p(rapMin + xi(iy)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(ik));
            double sumCl  = IntegrateClever(1024, p, Indx(jy,jk));
            double sumCl2 = IntegrateClever2(1024, p, Indx(jy,jk));
            double sumClDGadd = IntegrateCleverDGLAPadd(256, p, Indx(jy,jk));
            double sumClDG  = IntegrateCleverDGLAP(256, p, Indx(jy,jk));
            cout << "Radek "<<jy<<" "<<jk<<": " <<exp(-p.y)<< " " << exp(p.K)<<" : "<< sum << " "<< sumCl <<" "<<sumCl2<<" | "<< sumClDGadd <<" "<< sumClDG <<  endl;
            */

        }
        //cout << "Progress " << jGlob<<" "<< endl;
    }
    return convMat;
}





void Filler::Init()
{
    kerObj = new Kernel();
    kerObj->mu2 = exp(Kmin);
    kerObj->eps = 1e-7;
    kerObj->LnFreeze2 = log(1);
}



map<Point,double> Filler::KernelGeneral1D(Point p, Point q, double wR, double wS)
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double l = exp(0.5*p.K);
    double lp = exp(0.5*q.K);
    
    double z = x/xp;

    assert(z <= 1);

    Point qq;
    qq.y = q.y;
    qq.K = p.K;

    //double kerReg = kerObj->BFKL__Eps_Off(l, lp, z);
    //double kerSin = kerObj->BFKL__Eps_Diag(l, lp, z);
    double kerReg = kerObj->BFKL_res_kc_full__Eps_Off(l, lp, z);
    double kerSin = kerObj->BFKL_res_kc_full__Eps_Diag(l, lp, z);

    if(!isfinite(kerReg)) {
        cout <<"Helenka "  << l << " "<< lp <<" "<<  z << endl;
        cout << kerReg << endl;
        cout << kerSin << endl;
        assert(0);
    }
    assert(isfinite(kerReg));
    assert(isfinite(kerSin));

    map<Point,double> res;
    res[q]  += kerReg*wR;
    res[qq] += kerSin*wR;

    return res;
}

pair<double,double> Filler::KernelSpecific1D(Indx idq, Point p, Point q)
{
    double x  = exp(-p.y);
    double xp = exp(-q.y);

    double l = exp(0.5*p.K);
    double lp = exp(0.5*q.K);
    
    double z = x/xp;

    assert(z <= 1+kerObj->eps); 
    

    Point qq;
    qq.y = q.y;
    qq.K = p.K;

    double Fr = F(idq, q);
    double Frsub = F(idq, qq);

    //double kerReg = kerObj->BFKL__Eps_Off(l, lp, z);
    //double kerSin = kerObj->BFKL__Eps_Diag(l, lp, z);
    double kerReg = kerObj->BFKL_res_kc_full__Eps_Off(l, lp, z);
    double kerSin = kerObj->BFKL_res_kc_full__Eps_Diag(l, lp, z);

    return {kerReg*Fr + kerSin*Frsub, 0};
}

void saveAsHist(const arma::mat &m)
{
    
    TH2D *hMat = new TH2D("Mat", "Mat", m.n_rows, 0.5, 0.5+m.n_rows,    m.n_cols, 0.5, 0.5+m.n_cols);

    for(int i = 0; i < m.n_rows; ++i)
    for(int j = 0; j < m.n_cols; ++j)
        hMat->SetBinContent(i+1,j+1, m(i,j));
    hMat->SaveAs("mat.root");
}



int main(int argc, char **argv)
{
    int nCh = 33;

    long long nSplit = 1, nNow = 0;
    string kerName = argv[1];
    if(argc >= 3) nSplit = atoi(argv[2]);
    if(argc >= 4) nNow   = atoi(argv[3]);


    Filler filler;
    filler.Init();
    arma::mat evMat = filler.EvalAll(kerName, nCh,nSplit,nNow);

    //Transform to Nodes->Nodes
    arma::mat mCoef1d = GetCoefsCheb(nCh);
    arma::mat mCoef = arma::kron(mCoef1d, mCoef1d);
    arma::mat resMat = evMat * mCoef;


    saveAsHist(resMat);
    return 0;




    /*
    arma::vec val1(xi.n_rows), val2(xi.n_rows);
    for(int i = 0; i < xi.n_rows; ++i) {
        double x = xi(i);
        val1(i) = 5*Cheb(0,x) + Cheb(1,x) + 3*Cheb(2,x) + 8*Cheb(8,x);
        val2(i) = 2*Cheb(1,x) + 4*Cheb(2,x);
    }
    arma::vec val2d = arma::kron(val1, val2);
    */



    arma::vec xi   = GetNodes(nCh);
    //Check
    auto fun =  [](double x) {
        double K = Kmin + (Kmax-Kmin) * x;
        double kt2 = exp(K);
        return kt2*exp(-kt2);
    };

    arma::vec vals(xi.n_rows);
    for(int i = 0; i < xi.n_rows; ++i) {
        vals(i) = fun(xi(i));
        cout << "orgVal " << i << " "<< vals(i) << endl;
    }

    arma::vec resGl(nCh*nCh);

    for(int iGl = 0; iGl < nCh*nCh; ++iGl) {
        auto ind = toLocal(nCh, iGl);
        resGl(iGl) = vals(ind.iK);
    }

    /*
    cout << mCoef1d   << endl;
    cout << "papa" << endl;
    cout << mCoef1d -1  << endl;
    */

    cout <<  arma::solve(arma::eye<arma::mat>(resMat.n_rows,resMat.n_cols) - resMat, resGl) << endl;


}
