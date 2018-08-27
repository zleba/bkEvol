#include <iostream>
#include "TFile.h"
#include "TH2D.h"

#include <armadillo>
#include "nodes.h"

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
    //x = vr + mEq*x
    //(1-mEq)*x = vr

    return arma::solve(1-mEq, vr);
}


//arma::mat GetCoefs(int oldSize, bool isInverse = false);


void test()
{
    cout << "Start " << endl;
    arma::mat trMat = GetCoefs(11);
    arma::vec xi    = GetNodes(11);
    arma::vec vals(11);
    xi = 2*xi - 1;
    cout << "Coef calculated " << endl;

    auto fun =  [](double x) { return 2*x*x-1 + 2;};

    for(int i = 0; i < xi.n_rows; ++i)
        vals(i) = fun(xi(i));

    arma::vec res = trMat*vals;

    cout << "Result" << endl;
    cout << res << endl;

}


int main()
{
    test();
    return 0;

    TFile *file = TFile::Open("farm/histos/Sol98.root");

    TH2D *hMat = nullptr;
    file->GetObject("EqMat", hMat);

    arma::mat mEq = histo2mat(hMat);



    return 0;
}
