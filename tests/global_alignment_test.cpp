#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "global_alignment.hpp"


TEST_CASE( "Factorials are computed", "[factorial]" ) {
        REQUIRE( Factorial(1) == 1 );
}
