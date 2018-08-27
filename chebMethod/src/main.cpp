#include "cheb.h"
#include <iostream>

#include "TH2D.h"
#include "TFile.h"

using namespace std;

Indx toLocal(int Ncheb, int I)
{
    return Indx(I/Ncheb, I%Ncheb);
}

TH2D *FillEqMat(long long Ncheb, long long nSplit, long long nNow)
{
    Integrator integrator(50);

    long long pSize = Ncheb*Ncheb;
    long long nItem = pSize*pSize;

    TH2D *hEqMat = new TH2D("EqMat", "EqMat", pSize, 0.5, 0.5+pSize,    pSize, 0.5, 0.5+pSize);

    long long iGlob = -1;
    for(long long pI = 0; pI < pSize; ++pI) 
    for(long long qI = 0; qI < pSize; ++qI) {
        ++iGlob;
        auto iStart = (nItem*nNow)/nSplit;
        auto iEnd   = (nItem*(nNow+1))/nSplit;

        if(! (iStart <= iGlob && iGlob < iEnd)) continue;

        //cout << iGlob << endl;
        //continue;

        Indx idp(toLocal(Ncheb, pI));
        Indx idq(toLocal(Ncheb, qI));
        double val = integrator.Eval(idp,idq);
        hEqMat->SetBinContent(pI+1,qI+1, val);
        cout <<idp.iy <<" "<< idp.iK <<" : "<< idq.iy <<" "<< idq.iK << " <> "<< val << endl;

    }


    return hEqMat;
}

int main(int argc, char **argv)
{
    long long nSplit = 1, nNow = 0;
    
    if(argc >= 2) nSplit = atoi(argv[1]);
    if(argc >= 3) nNow   = atoi(argv[2]);

    TFile *file = TFile::Open("sol.root", "RECREATE");
    TH2D *hEqMat = FillEqMat(10, nSplit, nNow);
    file->Write();
    file->Close();
    

    return 0;

    //IntegrateMC();
    //IntegrateDummy();
    IntegrateClever();
    Integrator integrator(50);
    cout << integrator.Eval(Indx(2,1), Indx(2,1)) << endl;
    const int Ncheb = 4;
    for(int pI = 0; pI < Ncheb*Ncheb; ++pI) 
    for(int qI = 0; qI < Ncheb*Ncheb; ++qI) {
        Indx idp(toLocal(Ncheb, pI));
        Indx idq(toLocal(Ncheb, qI));
        cout <<idp.iy <<" "<< idp.iK <<" : "<< idq.iy <<" "<< idq.iK << " <> "<<   integrator.Eval(idp,idq) << endl;
    }
    /*
    double y = 0.9;
    double x = y;
    for(int i = 0; i < 100; ++i) {
        x = y / (1-log(x));
        cout <<i <<" "<<hexfloat<< x << endl;
    }
    */

    return 0;
}
