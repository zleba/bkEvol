#include "catch.hpp"
#include "cheb.h"

TEST_CASE( "Tchebishev pols are computed", "[Cheb]" ) {
    for(double x = -1; x <= 1; x+= 0.1) {
        REQUIRE( Cheb(0,0) == 1 );
        REQUIRE( Approx(Cheb(1,x)) == x );
        REQUIRE( Approx(Cheb(2,x)) == 2*x*x -1 );
        REQUIRE( Approx(Cheb(3,x)) == 4*x*x*x - 3*x );
        REQUIRE( Approx(Cheb(10,x)) == 512*pow(x,10) - 1280*pow(x,8) +
                               1120*pow(x,6) - 400*pow(x,4) + 50*x*x - 1 );
   }
}

