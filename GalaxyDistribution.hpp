#include "GravitationalFunctions.hpp"
#include <random>

namespace GravitationalLensing
{
    //A percentage multiplied by this constant will be representated on a scale where 1.0 is 100%//
    const ScalarType PERCENT_TO_DECIMAL_CONSTANT = .01;

    struct Galaxy {
        Vector3 position_;
        ScalarType totalMass_, darkMatterCoefficent_;
    };



    template< typename StatisticalType >
    struct AbstractStatisticalPair {
        StatisticalType average;
        StatisticalType standardDeviation;
        std::normal_distribution< StatisticalType > ToDistribution() {
            return std::normal_distribution< StatisticalType >( average, standardDeviation );
        }
    };

    using StatisticalPair = const AbstractStatisticalPair< const ScalarType >;

    struct VectorDistributor {
        std::normal_distribution< ScalarType > x_, y_, z_;
    };

    VectorDistributor MakeDistributor( ScalarType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator );

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator );

    template< template< typename GalaxyType > typename GalaxiesContainerType = std::vector >
    struct GalaxyCluster
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
        constexpr static StatisticalPair DefaultGalaxyDarkMatterCoefficentDistribution{ 75.0, 5.0 };
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
        ScalarType darkMatterCoefficent_, clusterTotalAmbientMass_;
        GalaxiesType galaxies_;
        std::mt19937 generator;

        explicit GalaxyCluster( GalaxiesType galaxies, 
                ScalarType darkMatterCoefficent, Vector3 dimensions, 
                ScalarType clusterTotalAmbientMass = DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant,
                Vector3 position = DefaultPositionConstant ) : position_( position ), dimensions_( dimensions ),
                        darkMatterCoefficent_( darkMatterCoefficent ), clusterTotalAmbientMass_( clusterTotalAmbientMass_ ), 
                        galaxies_( galaxies ), generator( std::random_device{}() ) {
        }
        explicit GalaxyCluster(
                bool singleGalaxyCluster = true,
                StatisticalPair galaxyDarkMatterCoefficientDistribution = DefaultGalaxyDarkMatterCoefficentDistribution,
                StatisticalPair galaxyTotalMassDistribution = DefaultGalaxyTotalMassDistributionConstant,
                size_t numberOfGalaxies = DefaultNumberOfGalaxiesConstant, 
                StatisticalPair clusterDimensionsDistribution = DefaultGalaxyClusterDimensionDistributionConstant,
                ScalarType clusterTotalAmbientMass = DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant,
                StatisticalPair clusterPositionDistribution = DefaultPositionConstant ) : 
                        galaxyDarkMatterCoefficentDistribution_( galaxyDarkMatterCoefficientDistribution ), 
                        galaxyTotalMassDistribution_( galaxyTotalMassDistribution ), 
                        clusterDimensionsDistribution_( clusterDimensionsDistribution ), 
                        clusterTotalAmbientMass_( clusterTotalAmbientMass ), 
                        clusterPositionDistribution_( clusterPositionDistribution ), 
                        generator( std::random_device{}( ) )
        {

        }

        //( StatisticalPair NumberOfGalaxiesConstant, ScalarType GalaxyPositionStandardDeviationCoefficentConstant,
        //    StatisticalPair GalaxyDarkMatterCoefficentConstant, GalaxyCluster& galaxyCluster, std::mt19937& generator )
        void PlaceGalaxies( const size_t NumberOfGalaxiesToGenerateConstant, const bool ResetGalaxiesConstant = true )
        {
            auto distributor = MakeDistributor( galaxyPositionStandardDeviation_, dimensions_, generator );
            std::normal_distribution darkMatterDistribution = darkMatterCoefficent_.ToDistribution();
            if( ResetGalaxiesConstant == true ) {
                galaxies_.clear();
            }
            for( size_t i = 0; i < NumberOfGalaxiesToGenerateConstant; ++i )
            {
                galaxyCluster.galaxies_.push_back( Galaxy newGalaxy{
                        .darkMatterCoefficent_ = std::abs( PERCENT_TO_DECIMAL_CONSTANT * darkMatterDistribution( generator ) ),
                        .position_ = RandomPointWithinBounds( distributor, generator ) } );
            }
        }
        //Encapsulation.//
        StatisticalPair GetGalaxyDarkMatterCoefficentDistribution() {
            return galaxyDarkMatterCoefficentDistribution_;
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
        protected:
            StatisticalPair galaxyDarkMatterCoefficentDistribution_, 
                    galaxyTotalMassDistribution_, clusterDimensionsDistribution_, 
                    clusterPositionDistribution_;
            ScalarType galaxyPositionStandardDeviation_;


            //GalaxyCluster GenerateGalaxyCluster( StatisticalPair DarkMatterCoeffiecentsConstant,
            //    StatisticalPair GalaxyClusterDimensionsConstant, StatisticalPair GalaxyClusterPositionConstant, std::mt19937& generator )

            //According to the distributions//
            void RandomizeGalaxyClusterParameters()
            {
                //Selected to sound somewhat cool.//
                auto dimensionalGenerator = clusterDimensionsDistribution_.ToDistribution();
                auto positionGenrator = clusterPositionDistribution_.ToDistribution();
                Vector3 position{ positionGenrator( generator ), positionGenrator( generator ),
                        positionGenrator( generator ) };
                Vector3 dimensions{ std::abs( dimensionalGenerator( generator ) ),
                        std::abs( dimensionalGenerator( generator ) ), std::abs( dimensionalGenerator( generator ) ) };
                return GalaxyCluster{ .position_ = position, .dimensions_ = dimensions, .darkMatterCoefficent_ = std::abs(
                        PERCENT_TO_DECIMAL_CONSTANT * galaxyDarkMatterCoefficentDistribution_.ToDistribution()( generator ) ), .galaxies_ = GalaxyCluster::GalaxiesType() };
            }
    };
}
