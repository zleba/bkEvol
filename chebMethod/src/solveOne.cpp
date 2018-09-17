#include <iostream>
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"

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


double eval(arma::vec res, double x, double kt2) //TODO let it work
{
    double y = log(1./x);
    double K = log(kt2);

    double yr = (y-rapMin)/(rapMax-rapMin);
    double kr =  (K-Kmin)/(Kmax - Kmin);
    int nCh = round(sqrt(res.n_rows));
    double sum = 0;
    return sum;
}


arma::vec getInput(int nCh, function<double(double)> fun)
{
    arma::vec xi   = GetNodes(nCh);

    arma::vec vals(xi.n_rows);
    for(int i = 0; i < xi.n_rows; ++i) {
        double K = Kmin + (Kmax-Kmin) * xi(i);
        double kt2 = exp(K);
        vals(i) = fun(kt2);
        //cout << "orgVal " << i << " "<< vals(i) << endl;
    }

    arma::vec resGl(nCh*nCh);

    for(int iGl = 0; iGl < nCh*nCh; ++iGl) {
        auto ind = toLocal(nCh, iGl);
        resGl(iGl) = vals(ind.iK);
    }
    return resGl;
}

arma::mat getMatrix(TString fName)
{
    TFile *file = TFile::Open(fName);

    TH2D *hMat = nullptr;
    file->GetObject("Mat", hMat);
    assert(hMat);

    arma::mat mEq = histo2mat(hMat);
    file->Close();
    return mEq;
}





int main()
{
    arma::mat mEqDP    = getMatrix("farm/results/_DGLAP.root");
    arma::mat mEqDPadd = getMatrix("farm/results/_DGLAPadd.root");
    arma::mat mEqBFKL = getMatrix("farm/results/_fullBFKL.root");
    arma::mat mEq =  mEqBFKL;
    //arma::mat mEq = mEqDP + mEqDPadd + mEqBFKL;

    int nCh = round(sqrt(mEq.n_rows));
    arma::vec v0 = getInput(nCh, [](double kt2){return kt2*exp(-kt2);});

    arma::vec res = arma::solve(arma::eye<arma::mat>(mEq.n_rows,mEq.n_cols) - mEq, v0);
    cout << res  << endl;

    arma::vec xi = GetNodes(nCh);

    for(int iy = 0; iy < xi.n_rows; ++iy)
    for(int ikt = 0; ikt < xi.n_rows; ++ikt) {
        double y = rapMin + xi(iy)*(rapMax-rapMin);
        double K = Kmin + xi(ikt)*(Kmax-Kmin);

        double x   = exp(-y);
        double kt2 = exp(K);

        int iGl = nCh*iy + ikt;
        //cout << x << " "<< kt2 <<" : "<< res(iGl) <<" "<< v0(ikt) << endl;
        cout << x << " "<< kt2 <<"  "<< res(iGl) << endl;
        //cout << x << " "<< kt2 <<" : "<< eval(res, x, kt2) <<" "<< v0(ikt)<< endl;

    }


    return 0;
}



