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
    arma::mat EvalAll(int Ncheb, long long nSplit, long long nNow);
    std::map<Point,double> GetIntegral(Point p);

    pair<double,double> KernelSpecific1D(Indx idq, Point p, Point q);
    map<Point,double> KernelGeneral1D(Point p, Point q, double wR, double wS);
    double IntegrateClever(int N, Point p, Indx idq);
    double IntegrateClever2(int N, Point p, Indx idq);
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





arma::mat Filler::EvalAll(int Ncheb, long long nSplit, long long nNow)
{
    assert(Ncheb % 2 == 1);

    arma::vec wgts = GetWeights(Ncheb); //basic weights
    arma::vec xi   = GetNodes(Ncheb);  //basic nodes


    std::map<double,double> qY, qK;
    map<Point,double> points[100][100];

    for(int iy = 0; iy < Ncheb; ++iy) //Loop over cheb nodes
    for(int ik = 0; ik < Ncheb; ++ik) {
        Point p(rapMin + xi(iy)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(ik));
        points[iy][ik] = GetIntegral(p);

        for(auto &point : points[iy][ik]) {
            qY[point.first.y] = 0;
            qK[point.first.K] = 0;
        }
    }
    cout << "Kernels evaluated" << endl;




    arma::mat convMat(Ncheb*Ncheb, Ncheb*Ncheb);

    long long nItem  = Ncheb * Ncheb;
    for(int jy = 0; jy < Ncheb; ++jy) //Loop over cheb Polynomials
    for(int jk = 0; jk < Ncheb; ++jk) {
        int jGlob = jy*Ncheb + jk;
        auto iStart = (nItem*nNow)/nSplit;
        auto iEnd   = (nItem*(nNow+1))/nSplit;
        if(! (iStart <= jGlob && jGlob < iEnd)) continue;

        for(auto &m : qK) m.second = FK(jk, m.first);
        for(auto &m : qY) m.second = Fy(jy, m.first);

        for(int iy = 0; iy < Ncheb; ++iy) //Loop over cheb nodes
        for(int ik = 0; ik < Ncheb; ++ik) {
            int iGlob = iy*Ncheb + ik;
            double sum = 0;
            for(auto &v : points[iy][ik])
                //sum += v.ker * F(idq, v.q);
                sum += v.second *  qY[v.first.y]*qK[v.first.K];
            convMat(iGlob,jGlob) = sum;

            Point p(rapMin + xi(iy)*(rapMax-rapMin), Kmin + (Kmax-Kmin)*xi(ik));
            //double sumCl  = IntegrateClever(2048*2, p, Indx(jy,jk));
            //double sumCl2 = IntegrateClever2(512, p, Indx(jy,jk));
            //cout << "Radek "<<jy<<" "<<jk<<": " <<exp(-p.y)<< " " << exp(p.K)<<" : "<< sum << " "<< sumCl <<" "<<sumCl2<< endl;
        }
        cout << "Progress " << jGlob<<" "<< endl;
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

    double kerReg = kerObj->BFKL__Eps_Off(l, lp, z);
    double kerSin = kerObj->BFKL__Eps_Diag(l, lp, z);

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
    if(argc >= 2) nSplit = atoi(argv[1]);
    if(argc >= 3) nNow   = atoi(argv[2]);


    Filler filler;
    filler.Init();
    arma::mat evMat = filler.EvalAll(nCh,nSplit,nNow);

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
