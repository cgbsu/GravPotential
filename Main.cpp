#include "Testing.hpp"

using namespace GravitationalLensing;

int main( int argc, char** args )
{
    std::cout.precision( 16 );
    GalaxyGenerationTest();
    //Testing::DemoDerivatives();
    //TestRandomFunctionDerivative();
    //DumpGravitationalPotentialDerivativesForTest();
    return 0;
}
