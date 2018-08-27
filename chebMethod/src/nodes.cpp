#include <armadillo>
#include <cassert>

//For points between 0 and 1
arma::vec GetWeights(int Size)
{
    //const int Size = 33;
    const int N = Size - 1;
    assert(N%2 == 0);

    //vector<vector<double>> coef(Size);
    //for(auto & c: coef) c.resize(Size);
    arma::mat coef(Size, Size);

    for(int k = 0; k <= N/2; ++k) {
        double s = 0;
        coef(2*k,N) = 1./N;
        coef(2*k,0) = 1./N ;

        coef(2*k,N/2) = 2./N * (2* ((k+1)%2) -1);

        for(int n = 1; n <= N/2-1; ++n)
            coef(2*k,n) = coef(2*k,N-n) = 2./N * cos(n*k*M_PI*2/N);
    }

    //vector<double> wgt(Size, 0.);
    arma::vec wgt(Size, arma::fill::zeros);

    for(int i = 0; i < Size; ++i) {
        wgt(i) += coef(0,i);
        wgt(i) += coef(N,i) / (1. - N*N);
        for(int k = 1; k <= N/2 - 1; ++k) {
           double w = 2./(1 - 4*k*k);
           wgt(i) += w*coef(2*k,i);
        }

        wgt(i) *= 0.5; //for interval (0,1)
    }
    return wgt;
    //for(int i = 0; i < Size; ++i)
        //cout << i <<" "<<setprecision(17)<< wgt[i] << endl;
}

arma::vec GetNodes(int Size)
{
    assert((Size - 1) % 2 == 0);
    arma::vec xi(Size, arma::fill::zeros);
    for(int i = 0; i < Size; ++i) {
      double Cos = cos(i /(Size-1.) * M_PI);
      xi(i) = (1-Cos)/2;
    }
    return xi;
}

