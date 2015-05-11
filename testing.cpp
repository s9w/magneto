#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "helpers.h"

TEST_CASE( "linspace" ) {
    std::vector<double> result = linspace(0.0, 1.0, 1);
    REQUIRE( result.size() == 1 );

    result = linspace(0.0, 1.0, 3);
    REQUIRE( result.size() == 3 );
    REQUIRE( result[0] == Approx( 0.0 ) );
    REQUIRE( result[1] == Approx( 0.5 ) );
    REQUIRE( result[2] == Approx( 1.0 ) );
}

TEST_CASE( "normalspace" ) {
    std::vector<double> result = normalSpace(0.1, 4.6, 1, 2.2, 0.7);
    REQUIRE( result.size() == 1 );
    REQUIRE( result[0] == Approx( 0.1 ) );

    result = normalSpace(0.1, 4.6, 3, 2.2, 0.7);
    REQUIRE( result.size() == 3 );
    REQUIRE( result[0] == Approx( 0.1 ) );
    REQUIRE( result[1] == Approx( 2.20091813 ) );
    REQUIRE( result[2] == Approx( 4.6 ) );
}