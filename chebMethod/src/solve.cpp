#include <iostream>
#include "TFile.h"
#include "TH2D.h"

#include <armadillo>
#include "nodes.h"
#include "cheb.h"
#include <iostream>

using namespace std;

arma::mat histo2mat(TH2D *h)
{
    arma::mat m(h->GetNbinsX(), h->GetNbinsY());

    for(int i = 0; i < m.n_rows; ++i)
    for(int j = 0; j < m.n_cols; ++j)
        m(i,j) = h->GetBinContent(i+1,j+1);

    return m;
}





arma::vec solve(const arma::mat &mEq, const arma::mat &vr)
{
    assert(mEq.n_rows == mEq.n_cols);
    assert(mEq.n_rows == vr.n_rows);
    //x = vr + mEq*x
    //(1-mEq)*x = vr

    arma::mat C(mEq.n_rows, mEq.n_cols);
    for(int i = 0; i < mEq.n_rows; ++i)
    for(int j = 0; j < mEq.n_cols; ++j) {
        auto indI = toLocal(mEq.n_rows, i);
        auto indJ = toLocal(mEq.n_rows, j);

        double resRap = GetChebIntegral(indI.iy, indJ.iy);
        double resK   = GetChebIntegral(indI.iK, indJ.iK);

        resRap *= rapMax-rapMin;
        resK   *= Kmax-Kmin;

        C(i,j) = resRap*resK;
    }

    return arma::solve(C - mEq, C*vr);
    //return arma::solve(1-mEq, vr);
}


//arma::mat GetCoefs(int oldSize, bool isInverse = false);


arma::vec test(int nCh2 = 21*21)
{
    int nCh = round(sqrt(nCh2));
    cout << "Start " << endl;
    arma::mat trMat = GetCoefs(nCh); //zero coef must be divided by two
    arma::mat trMatI= GetCoefs(nCh, true); //Inverse
    arma::vec xi    = GetNodes(nCh);
    arma::vec vals(nCh);
    //xi = 2*xi - 1;
    cout << "Coef calculated " << endl;

    auto fun =  [](double x) {
        double K = Kmin + (Kmax-Kmin) * x;
        double kt2 = exp(K);
        return kt2*exp(-kt2);
        //return pow(x,25);
        //return x*x*x*x -x*x + x + 2;
        //return rand()/(RAND_MAX+0.);
        //return x*exp(-120*x);
        x = 2*x-1;
        return Cheb(0,x)+Cheb(1,x)+Cheb(2,x)+
                         Cheb(3,x)+Cheb(4,x)
            ;
    };

    for(int i = 0; i < xi.n_rows; ++i) {
        vals(i) = fun(xi(i));
        cout << "orgVal " << i << " "<< vals(i) << endl;
    }

    arma::vec res  = trMat*vals;
    arma::vec resI = trMatI*res;

    //cout << "Helenka start" << endl;
    //cout << trMatI*trMat << endl;
    //cout << "Helenka end" << endl;

    res(0) /= 2;
    res(res.n_rows-1) /= 2;

    /*
    cout << "Result" << endl;
    cout << res << endl;
    cout << resI << endl;
    cout << vals << endl;
    for(int i = 0; i < xi.n_rows; ++i) {

        double sum = 0;
        for(int k = 0; k < nCh; ++k) {
           sum += res(k)*Cheb(k, (2*xi(i)-1)); 
        }


        cout << "Radek " << xi(i) <<" : "<< resI(i)<<" "<<vals(i) << "|"<< sum<< endl;
    }
    */

    //--nCh;
    arma::vec resGl(nCh*nCh);

    for(int iGl = 0; iGl < nCh*nCh; ++iGl) {
        auto ind = toLocal(nCh, iGl);
        if(ind.iy == 0)
            resGl(iGl) = res(ind.iK);
        else
            resGl(iGl) = 0;
    }

    return resGl;

    /*
    for(int i = 0; i < xi.n_rows; ++i) {
        double sum = 0;
        for(int k = 0; k < nCh; ++k) {
           sum += res(k)*Cheb(k, (2*xi(i)-1)); 
        }
        cout << xi(i) <<" "<< sum <<" "<<  fun(xi(i)) << endl;
    }
    */

    //Cheb(int n, double x)

}

double eval(arma::vec res, double x, double kt2)
{
    double y = log(1./x);
    double K = log(kt2);

    double yr = (y-rapMin)/(rapMax-rapMin);
    double kr =  (K-Kmin)/(Kmax - Kmin);

    int nCh = sqrt(res.n_rows);
    double sum = 0;
    for(int i = 0; i < res.n_rows; ++i) {
        auto ind = toLocal(nCh, i);

        sum += res(i) * Cheb(ind.iy, 2*yr-1) * Cheb(ind.iK, 2*kr-1);
    }
    return sum;
}



int main()
{
    //cout << v0 << endl;
    //return 0;

    TFile *file = TFile::Open("farm/results/bfkl21_1.root");

    TH2D *hMat = nullptr;
    file->GetObject("EqMat", hMat);

    arma::mat mEq = histo2mat(hMat);
    arma::vec v0 = test(mEq.n_rows);

    arma::vec res = solve(mEq, v0);

    cout << res << endl;

    arma::vec xi = GetNodes(round(sqrt(mEq.n_rows)));


    for(int iy = 0; iy < xi.n_rows; ++iy)
    for(int ikt = 0; ikt < xi.n_rows; ++ikt) {
        double y = rapMin + xi(iy)*(rapMax-rapMin);
        double K = Kmin + xi(ikt)*(Kmax-Kmin);

        double x   = exp(-y);
        double kt2 = exp(K);


        cout << x << " "<< kt2 <<" : "<< eval(res, x, kt2) <<" "<< eval(v0,x,kt2) << endl;
        //cout << x << " "<< kt2 <<" : "<< eval(res, x, kt2) <<" "<< v0(ikt)<< endl;

    }


    return 0;
}
