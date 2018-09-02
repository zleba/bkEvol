#include "catch.hpp"
#include "cheb.h"
#include <iostream>
#include <armadillo>
using namespace std;

TEST_CASE( "Tchebishev pols are computed", "[Cheb]" ) {
    for(double x = 0; x <= 1; x+= 0.1) {
        CHECK( Cheb(0,0) == 1 );
        CHECK( Approx(Cheb(13,x)).margin(1e-10) == cos(13*acos(2*x-1)) );
        CHECK( Approx(Cheb(16,x)).margin(1e-10) == cos(16*acos(2*x-1)) );
   }
}

TEST_CASE( "Integration by Gauss-Conrod", "[Integral61]" ) {
    double err;
    CHECK( Approx(1./3) == Integral61([](double x){return x*x;},0, 1, err)   );
    CHECK( Approx(2)    == Integral61([](double x){return sin(x);},0, M_PI, err)   );
    CHECK( Approx(1-exp(-1))  == Integral61([](double x){return exp(-x);},0, 1, err)   );
}

TEST_CASE( "Integration by Chebishev-nodes", "[ChebIntegrationMethod]" ) {
    int N = 128;
    arma::vec wgts = GetWeights(N+1); //basic weights
    arma::vec xi   = GetNodes(N+1);
    REQUIRE (sum(wgts) == Approx(1));

    auto Int = [&](function<double(double)> fun) {
        double sum = 0;
        for(int i = 0; i < xi.n_rows; ++i)
            sum += fun(xi(i))*wgts(i);
        return sum;
    };

    double err;
    auto pol2   = [](double x){return x*x;};
    CHECK( Approx(Int(pol2)) == Integral61(pol2, 0, 1, err)   );
    auto polGen = [](double x){return 4*x*x*x- 2*x*x + 20*x -8;};
    CHECK( Approx(Int(polGen)) == Integral61(polGen, 0, 1, err)   );
    auto Sin = [](double x){return sin(x);};
    CHECK( Approx(Int(Sin)) == Integral61(Sin, 0, 1, err)   );

    auto Sin20 = [](double x){return sin(20*x);};
    CHECK( Approx(Int(Sin20)) == Integral61(Sin20, 0, 1, err)   );

}



TEST_CASE( "Get Chebishev integral basis", "[GetChebIntegral]" ) {
    double err;
    CHECK( Approx(GetChebIntegral(0,0)  ) == Integral61([](double x){return Cheb(0,x)*Cheb(0,x);},0.0, 1, err)   );
    CHECK( Approx(GetChebIntegral(1,3)  ) == Integral61([](double x){return Cheb(1,x)*Cheb(3,x);},0.0, 1, err)   );
    CHECK( Approx(GetChebIntegral(3,3)  ) == Integral61([](double x){return Cheb(3,x)*Cheb(3,x);},0, 1, err)   );
    CHECK( Approx(GetChebIntegral(15,5)  ) == Integral61([](double x){return Cheb(15,x)*Cheb(5,x);},0, 1, err)   );
}



TEST_CASE("Calculate Cheb Polynomials", "[ChebPolynomials]") {
    CHECK( getChebCoef(3, [](double x) {return Cheb(0,x);})(0) == Approx(1.).margin(1e-10) );
    CHECK( getChebCoef(3, [](double x) {return Cheb(0,x);})(1) == Approx(0.).margin(1e-10) );
    CHECK( getChebCoef(3, [](double x) {return Cheb(0,x);})(2) == Approx(0.).margin(1e-10) );


    CHECK( getChebCoef(5, [](double x) {return Cheb(0,x)+2*Cheb(1,x)+3*Cheb(2,x)+4*Cheb(3,x);})(0) == Approx(1.).margin(1e-10) );
    CHECK( getChebCoef(5, [](double x) {return Cheb(0,x)+2*Cheb(1,x)+3*Cheb(2,x)+4*Cheb(3,x);})(1) == Approx(2.).margin(1e-10) );
    CHECK( getChebCoef(5, [](double x) {return Cheb(0,x)+2*Cheb(1,x)+3*Cheb(2,x)+4*Cheb(3,x);})(2) == Approx(3.).margin(1e-10) );
    CHECK( getChebCoef(5, [](double x) {return Cheb(0,x)+2*Cheb(1,x)+3*Cheb(2,x)+4*Cheb(3,x);})(3) == Approx(4.).margin(1e-10) );

}

