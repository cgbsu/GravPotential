#include "Testing.hpp"
#include <iomanip>
#include <ctime>
#include <chrono>

namespace GravitationalLensing
{
    std::string TimePointToString( std::chrono::system_clock::time_point from ) {
        auto now = std::chrono::system_clock::to_time_t( from );
        return std::string( std::ctime( &now ) );
    }

    std::string FormatTimePointString( std::string timePoint )
    {
        const size_t StringLengthConstant = timePoint.size();
        for( size_t i = 0; i < StringLengthConstant; ++i ) {
            if( std::isspace( timePoint[ i ] ) != 0 || timePoint[ i ] == ':' )
                timePoint.replace( i, 1, "_" );
        }
        return timePoint;
    }

    std::string TimeMarkedFileName( std::string appendName ) {
        return ( FormatTimePointString( 
                TimePointToString( std::chrono::system_clock::now() ) ) + appendName );
    }

    void DumpGravitationalPotentialForTest()
    {
        std::ofstream potentialData;
        potentialData.precision( 16 );
        potentialData.open( TimeMarkedFileName( "PotentialData.csv" ) );
        //Set X to value close to zero, so there is no divide by zero error.//
        Vector3 position{ 10e-30, 0.0, 0.0 };
        constexpr size_t AmountOfSteps = 1000;
        constexpr ScalerType StepConstant = Constants< ScalerType >::KiloparsecInMetersConstant;
        constexpr ScalerType ScaledRadiusConstant = 100.0 * Constants< ScalerType >::KiloparsecInMetersConstant;
        constexpr ScalerType TidalRadiusConstant = 3.0 * ScaledRadiusConstant;
        constexpr ScalerType MassConstant = 1e12;
        GravitationalPotential PotentialConstant( TidalRadiusConstant, ScaledRadiusConstant, MassConstant );
        for( size_t i = 0; i < AmountOfSteps; ++i )
        {
            ScalerType result = PotentialConstant( position );
            std::cout << position.x_ << ", " << result << "\n";
            potentialData << result << ",";
            position.x_ += StepConstant;
        }
        potentialData.close();
    }

    void DumpGravitationalPotentialDerivativesForTest()
    {
        std::ofstream potentialData;
        potentialData.precision( 16 );
        potentialData.open( TimeMarkedFileName( "2DDerivativePotentialData.csv" ) );
        const ScalerType StartConstant = 10e-30;
        //Set X to value close to zero, so there is no divide by zero error.//
        Vector3 position{ StartConstant, StartConstant, 0.0 };
        constexpr size_t AmountOfSteps = 1000;
        constexpr ScalerType StepConstant = Constants< ScalerType >::KiloparsecInMetersConstant;
        constexpr ScalerType ScaledRadiusConstant = 100.0 * Constants< ScalerType >::KiloparsecInMetersConstant;
        constexpr ScalerType TidalRadiusConstant = 3.0 * ScaledRadiusConstant;
        constexpr ScalerType MassConstant = 1e12 * 20e30;
        GravitationalPotential PotentialConstant( TidalRadiusConstant, ScaledRadiusConstant, MassConstant );
        std::function potential{ PotentialConstant };
        Deriver< decltype( potential ) > deriver( potential, .1, position );
        for( size_t i = 0; i < AmountOfSteps; ++i )
        {
            #ifdef _WIN32
                std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
            #endif
            position.y_ = StartConstant;
            for( size_t j = 0; j < AmountOfSteps; ++j )
            {
                ScalerType result = PotentialConstant( position );
                deriver.Derive( position );
                potentialData << position.x_ << ", " << position.y_ << ", " <<
                        result << ", " << deriver.x() << ", " << deriver.y() << ", " << deriver.x.y() << "\n";
                position.y_ += StepConstant;
            }
            position.x_ += StepConstant;
            std::cout << ( int ) ( ( ( double ) i / ( double ) AmountOfSteps ) * 100.0 ) << "%";
        }
        std::cout << "\nDONE\n";
        potentialData.close();
    }

    ScalerType WeirdFunction( Vector3 vector ) {
        return std::sin( vector.x_ ) + std::sin( vector.y_ ) + std::sin( vector.z_ );
    }

    void TestRandomFunctionDerivative()
    {
        const ScalerType StartConstant = 0.0; //10e-30;
        constexpr ScalerType StepConstant = .1;//Constants< ScalerType >::KiloparsecInMetersConstant;
        //Set X to value close to zero, so there is no divide by zero error.//
        Vector3 position{ StartConstant, StartConstant, 0.0 };
        constexpr size_t AmountOfSteps = 1000;
        Deriver deriver( &WeirdFunction, .1, position );
        for( size_t i = 0; i < AmountOfSteps; ++i )
        {
            std::cout << deriver.x() << "\n";
            deriver.Derive( position );
            position.x_ += StepConstant;
        }
    }

    void GeneralTest0()
    {
        std::cout << "Lens Radial Acceleration " << LensRadialAcceleration()( 1.0 / Constants< ScalerType, 0 >::HubbleNaughtConstant ) << "\n";
        std::cout << "Gravitational Potential " << GravitationalPotential( 1.0, 1.0, 1.0 )( Vector3{ 1.0, 1.0, 1.0 } ) << "\n";
        StatisticalPair DarkMatterCoeffiecentsConstant{ .average = 75.0, .standardDeviation = 5.0 };
        StatisticalPair GalaxyClusterDimensionsConstant{ .average = 10000.0, .standardDeviation = 10000.0 };
        StatisticalPair GalaxyClusterPositionConstant{ .average = 0.0, .standardDeviation = 100000.0 };
        std::random_device randomDevice;
        std::mt19937 generator{ randomDevice() };
        GalaxyCluster example = GenerateGalaxyCluster( DarkMatterCoeffiecentsConstant,
            GalaxyClusterDimensionsConstant, GalaxyClusterPositionConstant, generator );
        //    GalaxyCluster example{ .position_ = Vector3{.x_ = 0.0, .y_ = 0.0, .z_ = 0.0 },
        //            .dimensions_ = Vector3{.x_ = 10000.0, .y_ = 100000.0, .z_ = 5000.0 }, .darkMatterCoefficent_ = 70.0, .galaxies_ = GalaxyCluster::GalaxiesType() };
        StatisticalPair NumberOfGalaxiesConstant{ .average = 1000.0, .standardDeviation = 100.0 };
        StatisticalPair GalaxyDarkMatterCoefficentConstant{ .average = 75.0, .standardDeviation = 5.0 };
        ScalerType GalaxyPositionStandardDeviationCoefficentConstant = 1.0;
        PlaceGalaxies( NumberOfGalaxiesConstant, GalaxyPositionStandardDeviationCoefficentConstant,
            GalaxyDarkMatterCoefficentConstant, example, generator );
        std::stringstream fileName;
        fileName << "galaxyClusterData" << TimePointToString( std::chrono::system_clock::now() ) << ".csv";
        std::ofstream galaxyClusterData;
        galaxyClusterData.open( "Cluster.txt" );
        std::cout << "Galaxy Cluster With Dimensions x: " <<
            example.dimensions_.x_ << " y: " << example.dimensions_.y_ << " z: " << example.dimensions_.z_ <<
            " and dark matter coefficent " << example.darkMatterCoefficent_ << "\n";
        std::cout << "There are " << example.galaxies_.size() << " galaxies. Their positions are: \n";
        galaxyClusterData << example.galaxies_.size() << ", " << example.dimensions_.x_ << ", " << example.dimensions_.y_ << ", " << example.dimensions_.z_ << ", ";
        for( auto& galaxy : example.galaxies_ )
        {
            std::cout << "    x: " << galaxy.position_.x_ << " y: " <<
                galaxy.position_.y_ << " z: " << galaxy.position_.z_ <<
                " total dark matter percentage: " << galaxy.darkMatterCoefficent_ << "\n";
            galaxyClusterData << galaxy.position_.x_ << ", " << galaxy.position_.y_ << ", "
                << galaxy.position_.z_ << ", " << galaxy.darkMatterCoefficent_ << ", ";
        }
        galaxyClusterData.close();
    }
}
