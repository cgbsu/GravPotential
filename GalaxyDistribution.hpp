#include "GravitationalFunctions.hpp"
#include <random>

namespace GravitationalLensing
{
    //A percentage multiplied by this constant will be representated on a scale where 1.0 is 100%//
    const ScalarType PercentToDecimalConstant = .01;

    //Total mass assumed to be in Kg, TODO: Use physical units library?//
    struct Galaxy {
        Vector3 position_;
        ScalarType totalMass_;
    };

    template< typename StatisticalType >
    struct AbstractStatisticalPair {
        StatisticalType average;
        StatisticalType standardDeviation;
        std::normal_distribution< StatisticalType > ToDistributor() {
            return std::normal_distribution< StatisticalType >( average, standardDeviation );
        }
    };

    using StatisticalPair = const AbstractStatisticalPair< const ScalarType >;

    struct VectorDistributor {
        std::normal_distribution< ScalarType > x_, y_, z_;
    };

    VectorDistributor MakeDistributor( ScalarType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator );

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator );





    /*******************************************************************************************************************
    * TODO: Decouple the galaxy cluster data from the stuff that randomizes galaxies, moved randomization out of *******
    * functions into here because it was a lot of parameters and I wanted better dependancy injection for ease of use. *
    * This honestly is pretty bad though, there needs to be a good model where re-randomizing the clusters' parameters *
    * in particular dimensions, invalidates the galaxies
    *******************************************************************************************************************/
    template< template< typename GalaxyType > typename GalaxiesContainerType = std::vector >
    struct AbstractGalaxyCluster
    {

        constexpr static ScalarType GalaxyMassInSolarMassesLowerBoundConstant = 10e11;
        constexpr static ScalarType GalaxyMassInSolarMassesUpperBoundConstant = 10e13;
        constexpr static ScalarType GalaxyMassInSolarMassesLowerBoundExponentConstant = 12;
        constexpr static ScalarType GalaxyMassInSolarMassesUpperBoundExponentConstant = 14;
        constexpr static ScalarType DefaultGalaxyClusterDimensionConstant = 1500.0 * Constants::KiloparsecInMetersConstant;
        constexpr static StatisticalPair DefaultGalaxyClusterDimensionDistribution{ 
                .average = DefaultGalaxyClusterDimensionConstant, 
                .standardDeviation = DefaultGalaxyClusterDimensionConstant };
        constexpr static ScalarType DefaultGalaxyPositionStandardDeviationConstant = DefaultGalaxyDimensionConstant / 2.0;
        using GalaxiesType = GalaxiesContainerType< Galaxy >;
        constexpr static ScalarType DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant = 10e15;
        constexpr static ScalarType DefaultGalaxyPositionStandardDeviationConstant = 1.0;
        constexpr static Vector3 DefaultPositionConstant = Vector3{ 0.0, 0.0, 0.0 };
        constexpr static StatisticalPair DefaultGalaxyTotalMassDistributionConstant{
                .average = GalaxyMassInSolarMassesLowerBoundConstant,
                .standardDeviation = GalaxyMassInSolarMassesUpperBoundConstant };
        /* Its generally just a good idea when dealing with sizes of memory, **********
        * especially memory that will be moved around a lot, to uses multiples of 8. */
        constexpr static size_t DefaultNumberOfGalaxiesConstant = 1024;

        Vector3 position_, dimensions_;
        ScalarType clusterTotalAmbientMass_;
        GalaxiesType galaxies_;

        explicit AbstractGalaxyCluster( GalaxiesType galaxies, Vector3 dimensions, 
                ScalarType clusterTotalAmbientMass = DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant,
                Vector3 position = DefaultPositionConstant ) : position_( position ), dimensions_( dimensions ),
                        clusterTotalAmbientMass_( clusterTotalAmbientMass_ ), 
                        galaxies_( galaxies ), generator( std::random_device{}() ) {
            galaxyTotalMassDistributor_ = galaxyTotalMassDistribution_.ToDistributor();
        }
        explicit AbstractGalaxyCluster(
                bool singleGalaxyCluster = true,
                StatisticalPair galaxyTotalMassDistribution = DefaultGalaxyTotalMassDistributionConstant,
                size_t numberOfGalaxies = DefaultNumberOfGalaxiesConstant,
                StatisticalPair clusterDimensionsDistribution = DefaultGalaxyClusterDimensionDistributionConstant,
                ScalarType clusterTotalAmbientMass = DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant,
                StatisticalPair clusterPositionDistribution = DefaultPositionConstant ) :
                        galaxyTotalMassDistribution_( galaxyTotalMassDistribution ),
                        clusterDimensionsDistribution_( clusterDimensionsDistribution ),
                        clusterTotalAmbientMass_( clusterTotalAmbientMass ),
                        clusterPositionDistribution_( clusterPositionDistribution ),
                        generator( std::random_device{}() ),
        {
            galaxyTotalMassDistributor_ = galaxyTotalMassDistribution_.ToDistributor();
            RandomizeGalaxyClusterParameters();
            PlaceGalaxies( numberOfGalaxies, false );
        }

        //( StatisticalPair NumberOfGalaxiesConstant, ScalarType GalaxyPositionStandardDeviationCoefficentConstant,
        //    StatisticalPair GalaxyDarkMatterCoefficentConstant, GalaxyCluster& galaxyCluster, std::mt19937& generator )
        void PlaceGalaxies( const size_t NumberOfGalaxiesToGenerateConstant, const bool ResetGalaxiesConstant )
        {
            auto galaxyPositionDistributor = MakeDistributor( galaxyPositionStandardDeviation_, dimensions_, generator );
            if( ResetGalaxiesConstant == true ) {
                galaxies_.clear();
            }
            for( size_t i = 0; i < NumberOfGalaxiesToGenerateConstant; ++i )
            {
                galaxyCluster.galaxies_.push_back( Galaxy{
                        .position_ = RandomPointWithinBounds( galaxyPositionDistributor, generator )
                        .totalMass_ = ProduceRandomGalaxyTotalMass() * Constants::SolarMassInKilogramsConstant } );
            }
        }

        StatisticalPair GetGalaxyTotalMassDistribution() {
            return galaxyTotalMassDistribution_;
        }
        StatisticalPair GetClusterDimensionsDistribution() {
            return clusterDimensionsDistribution_;
        }
        StatisticalPair GetClusterPositionDistribution() {
            return clusterPositionDistribution_;
        }
        ScalarType GetGalaxyPositionStandardDeviation() {
            return galaxyPositionStandardDeviation_;
        }
        ScalarType SetGalaxyPositionStandardDeviation( ScalarType galaxyPositionStandardDeviation ) {
            galaxyPositionStandardDeviation_ = galaxyPositionStandardDeviation;
        }
        protected:
            std::mt19937 generator;
            StatisticalPair galaxyTotalMassDistribution_, clusterDimensionsDistribution_, 
                    clusterPositionDistribution_;
            ScalarType galaxyPositionStandardDeviation_;
            std::normal_distribution< ScalarType > galaxyTotalMassDistributor_;

            ScalarType ProduceRandomGalaxyTotalMassInSolarMasses() {
                //TODO: After MathFunctions is fixed, replace this with MathFucntions abs.//
                return std::abs( galaxyTotalMassDistributor_( generator ) );
            }

            //GalaxyCluster GenerateGalaxyCluster( StatisticalPair DarkMatterCoeffiecentsConstant,
            //    StatisticalPair GalaxyClusterDimensionsConstant, StatisticalPair GalaxyClusterPositionConstant, std::mt19937& generator )

            //According to the distributions//
            void RandomizeGalaxyClusterParameters()
            {
                //Selected to sound somewhat cool.//
                auto dimensionalGenerator = clusterDimensionsDistribution_.ToDistributor();
                auto positionGenrator = clusterPositionDistribution_.ToDistributor();
                position_ = Vector3{ positionGenrator( generator ), positionGenrator( generator ),
                        positionGenrator( generator ) };
                dimensions_ = Vector3{ std::abs( dimensionalGenerator( generator ) ),
                        std::abs( dimensionalGenerator( generator ) ), std::abs( dimensionalGenerator( generator ) ) };
            }
    };
}
