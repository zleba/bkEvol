#include <iostream>
#include <armadillo>
#include "cheb.h"
#include "nodes.h"

using namespace std;

//Final BK expression
arma::vec nonLinBKall(arma::mat nonKer, arma::vec fIn)
{
    fIn = square(fIn); //can be improved
    //Integrate over y
    arma::vec res = nonKer * fIn;
    return res;
}

//Get non-linear transformation
arma::vec nonLinBKTest(arma::vec fIn)
{
    //fIn.transform([](double x){return x*x;});

    int Ni = 21;
    int N = fIn.n_rows;
    arma::vec wgts = GetWeights(Ni); //basic weights
    arma::vec xi   = rapMin + (rapMax-rapMin)*GetNodes(N);

    arma::vec fOut(N);

    for(int i = 0; i < N; ++i) {
        double y = xi(i); //lover integration bound
        
        arma::vec xiNow = rapMin + (y-rapMin)*GetNodes(Ni);
        //cout << "Radek " << y / rapMax << endl;
        //arma::mat coef = GetCoefs(N, rapMin, y/rapMax);
        arma::mat coef = GetCoefsGeneric(N, Ni, rapMin, y/rapMax);

        arma::vec intVals =  coef *fIn;
        //Evaluate Kernel at these points
        fOut(i) = dot(wgts, intVals) * (y - rapMin);
        //cout << i <<" "<< rapMax << " "<< y << endl;
    }
    return fOut;
}

//get Integration matrix (kernel)
//It transform f -> int f
arma::mat nonLinBKmatTest(int N, int Nint)
{
    arma::vec wgts = GetWeights(Nint); //basic weights
    arma::vec xi   = rapMin + (rapMax-rapMin)*GetNodes(N);

    arma::mat mOut(N,N);

    for(int i = 0; i < N; ++i) {
        double y = xi(i); //lover integration bound
        
        arma::mat coef = GetCoefsGeneric(N, Nint, rapMin, y/rapMax);
        mOut.row(i) = wgts.t() * coef * (y - rapMin);
    }
    return mOut;
}




//get Integration matrix (kernel)
//It transform f -> int f
//Works just on the y-space
arma::mat nonLinBKmatBase(int N, int Nint)
{
    arma::vec wgts = GetWeights(Nint); //basic weights
    arma::vec xi   = rapMin + (rapMax-rapMin)*GetNodes(N);

    arma::mat mOut(N,N);

    for(int i = 0; i < N; ++i) {
        double y = xi(i); //lover integration bound
        
        arma::mat coef = GetCoefsGeneric(N, Nint, rapMin, y/rapMax);
        mOut.row(i) = wgts.t() * coef * (y - rapMin);
    }
    return mOut;
}

arma::mat nonLinBKmat(int N, int Nint)
{
    int nY = N;
    int nK = N;

    arma::mat yBase = nonLinBKmatBase(N, Nint);

    arma::mat mOut(nY*nK, nY*nK);

    for(int iy = 0; iy < nY; ++iy)
    for(int ik = 0; ik < nK; ++ik) {
        int iG = nK*iy + ik;
        for(int jy = 0; jy < nY; ++jy)
        for(int jk = 0; jk < nK; ++jk) {
            int jG = nK*jy + jk;
            mOut(iG, jG) = yBase(iy, jy);
        }
    }
    return mOut;
}








int main()
{
    arma::vec fIn(11);
    arma::vec xi   = rapMin + (rapMax-rapMin)*GetNodes(fIn.n_rows);
    for(int i = 0; i < fIn.n_rows; ++i)
        fIn(i) = 1;//xi(i)*xi(i);

    //fIn.transform([](double x){return 1;});

    arma::vec fOut = nonLinBKTest(fIn);
    cout << fOut << endl;

    cout << "Kernel based solution" << endl;
    cout << nonLinBKmatBase(11,21) * fIn << endl;


    xi.transform([](double x){return x;});
    cout << xi << endl;
    


    cout << "Hello" << endl;

    return 0;
}
