#include "GalaxyDistribution.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <sstream>
#include <ctime>
namespace GravitationalLensing
{
	std::string TimePointToString( std::chrono::system_clock::time_point from );
	std::string FormatTimePointString( std::string timePoint );
	std::string TimeMarkedFileName( std::string appendName );
	void DumpGravitationalPotentialForTest();
	void DumpGravitationalPotentialDerivativesForTest();
	void TestRandomFunctionDerivative();
	void GeneralTest0();
}
