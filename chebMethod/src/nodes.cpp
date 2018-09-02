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

arma::mat GetCoefs(int oldSize, bool isInverse = false)
{
    const int N = oldSize - 1;
    assert(N%2 == 0);

    arma::mat coef(oldSize,oldSize);

    double mul = 1;
    double C = 1./N;
    if(isInverse == true) {C = 1./2; }

    //isInverse = false;
    for(int k = 0; k <= N; ++k) {
        double s = 0;
        if(!isInverse) {
            coef(k,N) = C;
            coef(k,0) = C * (k % 2 == 1 ? -1 : 1);
        }
        else {
            mul = k % 2 == 1 ? -1 : 1;
            coef(N-k, N) = C * mul;
            coef(N-k, 0) = C ;
        }

        for(int n = 1; n <= N-1; ++n) {
            double el = cos(n*k*M_PI / N) * 2.*C * mul;
            if(!isInverse) coef(k,N-n) = el;
            else           coef(N-k,N-n) = el;
        }
    }
    
    return coef;
}

arma::vec getPols(int Size, double x)
{
    arma::vec pol(Size);

    if(Size >= 1) pol(0) = 1;
    if(Size >= 2) pol(1) = 2*x-1;

    for(int i = 2; i < Size; ++i)
        pol(i) = 2*(2*x-1)*pol(i-1) - pol(i-2);
    return pol;
}



double evalPol(const arma::vec &polCoef, double x)
{
    arma::vec pols = getPols(polCoef.n_rows, x);
    return arma::dot(pols,polCoef);
}

arma::vec interpol(arma::vec xi, double x)
{
    arma::vec coefs(xi.n_rows);
    for(int i = 0; i < xi.n_rows; ++i) {
        double num = 1, den = 1;
        for(int j = 0; j < xi.n_rows; ++j)
            if(j != i) {
                num *= x     - xi(j);
                den *= xi(i) - xi(j);
            }
        coefs(i) = num/den;
    }
    return coefs;
}

double interpol(arma::vec xi, arma::vec vals, double x)
{
    arma::vec coefs = interpol(xi, x);
    return arma::dot(coefs,vals);
}



//Transformation from chebNodes between 0 and 1 to chebNodes between a and b
arma::mat GetCoefs(int Size, double a, double b)
{
    arma::vec xi   = a + (b-a)*GetNodes(Size);

    arma::mat polsAll(Size,Size);
    for(int i = 0; i < Size; ++i)
        polsAll.row(i) = getPols(Size,xi(i));

    arma::mat coef = GetCoefs(Size); //now we have cheb pol coef
    return polsAll*coef;
}
