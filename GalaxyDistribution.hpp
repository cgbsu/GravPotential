#include "GravitationalFunctions.hpp"

#ifndef GRAVITATIONAL_LENSING_SIMULATION_GALAXY_DISTRIBUTION_HEARDER_H
#define GRAVITATIONAL_LENSING_SIMULATION_GALAXY_DISTRIBUTION_HEARDER_H
#include <random>

namespace GravitationalLensing
{
    //A percentage multiplied by this constant will be representated on a scale where 1.0 is 100%//
    const ScalarType PercentToDecimalConstant = .01;

    //Total mass assumed to be in Kg, TODO: Use physical units library?//
    struct Galaxy {
        Vector3 position_;
        ScalarType totalMass_, radius_;
    };

    template< typename StatisticalType >
    struct AbstractStatisticalPair
    {
        explicit AbstractStatisticalPair( StatisticalType average, StatisticalType standardDeviation ) : 
                average_( average ), standardDeviation_( standardDeviation ) {
        }
        AbstractStatisticalPair( const AbstractStatisticalPair< StatisticalType >& other ) = default;

        AbstractStatisticalPair< StatisticalType >& operator=( const AbstractStatisticalPair< StatisticalType >& other ) {
            average_ = other.GetAverage();
            standardDeviation_ = other.GetStandardDeviation();
            return *this;
        }

        std::normal_distribution< StatisticalType > ToDistributor() {
            return std::normal_distribution< StatisticalType >( average_, standardDeviation_ );
        }
        StatisticalType GetAverage() const {
            return average_;
        }
        StatisticalType GetStandardDeviation() const {
            return standardDeviation_;
        }
        void SetAverage( StatisticalType average ) {
            average_ = average;
        }
        void SetStandardDeviation( StatisticalType standardDeviation ) {
            standardDeviation_ = standardDeviation;
        }
        protected: 
            StatisticalType average_, standardDeviation_;
    };

    using StatisticalPair = AbstractStatisticalPair< ScalarType >;

    struct VectorDistributor {
        std::normal_distribution< ScalarType > x_, y_, z_;
    };

    VectorDistributor MakeDistributor( ScalarType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator );

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator );

    struct GalaxyClusterGenerationParameters
    {
        constexpr static ScalarType DefaultGalaxyClusterRadiusInMetersAverageConstant =
            1500.0 * Constants< ScalarType >::KiloparsecInMetersConstant;
        constexpr static ScalarType DefaultGalaxyClusterRaduisInMetersStandardDeviationConstant =
            100.0 * Constants< ScalarType >::KiloparsecInMetersConstant;

        constexpr static ScalarType DefaultGalaxyClusterMassInSolarMassesAverageConstant = 10e15;
        constexpr static ScalarType DefaultGalaxyClusterMassInSolarMassesStandardDeviationConstant = 10e14;

        constexpr static ScalarType DefaultGalaxyClusterUniformPositionConstant = 0.0;
        constexpr static ScalarType DefaultClusterPositionDistributionStandardDeviation =
            1e6 * Constants< ScalarType >::KiloparsecInMetersConstant;

        GalaxyClusterGenerationParameters(
                StatisticalPair clusterPositionDistribution = StatisticalPair{ 
                        DefaultGalaxyClusterUniformPositionConstant, 
                        DefaultClusterPositionDistributionStandardDeviation }, 
                StatisticalPair clusterRadiusDistribution = StatisticalPair{
                        DefaultGalaxyClusterRadiusInMetersAverageConstant,
                        DefaultGalaxyClusterRaduisInMetersStandardDeviationConstant },
                StatisticalPair clusterAmbientMassDistribution = StatisticalPair{
                        DefaultGalaxyClusterMassInSolarMassesAverageConstant,
                        DefaultGalaxyClusterMassInSolarMassesStandardDeviationConstant } ) :
                clusterRadiusDistribution_( clusterRadiusDistribution ), 
                clusterAmbientMassDistribution_( clusterAmbientMassDistribution ),
                clusterPositionDistribution_( clusterPositionDistribution ) {
        }

        StatisticalPair GetClusterRadiusDistribution() {
            return clusterRadiusDistribution_;
        }
        StatisticalPair GetClusterAmbientMassDistribution() {
            return clusterAmbientMassDistribution_;
        }
        StatisticalPair GetClusterPositionDistribution() {
            return clusterPositionDistribution_;
        }
        void SetClusterRadiusDistribution( StatisticalPair clusterRadiusDistribution ) {
            clusterRadiusDistribution_ = clusterRadiusDistribution_;
        }
        void SetClusterAmbientMassDistribution( StatisticalPair clusterAmbientMassDistribution ) {
            clusterAmbientMassDistribution_ = clusterAmbientMassDistribution;
        }
        void SetClusterPositionDistribution( StatisticalPair clusterPositionDistribution ) {
            clusterPositionDistribution_ = clusterPositionDistribution;
        }
    protected:
            StatisticalPair clusterRadiusDistribution_, clusterAmbientMassDistribution_, 
                    clusterPositionDistribution_;
    };

    struct GalaxyGenerationParameters
    {
        constexpr static ScalarType DefaultGalaxyMassInSolarMassesLowerBoundConstant = 10e11;
        constexpr static ScalarType DefaultGalaxyMassInSolarMassesUpperBoundConstant = 10e13;
        constexpr static ScalarType DefaultGalaxyRadiusInMetersAverageConstant = 
                150.0 * Constants< ScalarType >::KiloparsecInMetersConstant;
        constexpr static ScalarType DefaultGalaxyRadiusInMetersStandardDeviationConstant = 
                10 * Constants< ScalarType >::KiloparsecInMetersConstant;

        GalaxyGenerationParameters( 
                ScalarType clusterRadius, 
                StatisticalPair galaxyRadiusDistribution = StatisticalPair{ 
                        DefaultGalaxyRadiusInMetersAverageConstant, 
                        DefaultGalaxyRadiusInMetersStandardDeviationConstant }, 
                ScalarType galaxyMassInSolarMassesLowerBound = DefaultGalaxyMassInSolarMassesLowerBoundConstant,
                ScalarType galaxyMassInSolarMassesUpperBound = DefaultGalaxyMassInSolarMassesUpperBoundConstant ) : 
                        galaxyRadiusDistribution_( galaxyRadiusDistribution ), 
                        galaxyMassInSolarMassesLowerBound_( galaxyMassInSolarMassesLowerBound ), 
                        galaxyMassInSolarMassesUpperBound_( galaxyMassInSolarMassesUpperBound ) {
            SetGalaxyPositionStandardDeviationBasedOnClusterRadius( clusterRadius );
        }

        void SetGalaxyPositionStandardDeviationBasedOnClusterRadius( ScalarType clusterRadius ) {
            galaxyPositionStandardDeviation_ = clusterRadius / 2.0;
        }

        //------Encapsulation-----//
        StatisticalPair GetGalaxyRadiusDistribution() {
            return galaxyRadiusDistribution_;
        }
        ScalarType GetGalaxyMassInSolarMassesLowerBound() {
            return galaxyMassInSolarMassesLowerBound_;
        }
        ScalarType GetGalaxyMassInSolarMassesUpperBound() {
            return galaxyMassInSolarMassesUpperBound_;
        }
        ScalarType GetGalaxyPositionStandardDeviation() {
            return galaxyPositionStandardDeviation_;
        }
        void SetGalaxyRadiusDistribution( StatisticalPair galaxyRadiusDistribution ) {
            galaxyRadiusDistribution_ = galaxyRadiusDistribution;
        }
        void SetGalaxyPositionStandardDeviation( ScalarType galaxyPositionStandardDeviation ) {
            galaxyPositionStandardDeviation_ = galaxyPositionStandardDeviation;
        }
        void SetGalaxyMassInSolarMassesLowerBound( ScalarType galaxyMassInSolarMassesLowerBound ) {
            galaxyMassInSolarMassesLowerBound_ = galaxyMassInSolarMassesLowerBound;
        }
        void SetGalaxyMassInSolarMassesUpperBound( ScalarType galaxyMassInSolarMassesUpperBound ) {
            galaxyMassInSolarMassesUpperBound_ = galaxyMassInSolarMassesUpperBound;
        }
        protected:
            StatisticalPair galaxyRadiusDistribution_;
            ScalarType galaxyPositionStandardDeviation_, galaxyMassInSolarMassesLowerBound_, 
                    galaxyMassInSolarMassesUpperBound_;
    };

    template< template< typename GalaxyType > typename GalaxiesContainerType >
    struct AbstractGalaxyCluster
    {
        using GalaxiesType = GalaxiesContainerType< Galaxy >;
        AbstractGalaxyCluster( 
                ScalarType radius = GalaxyClusterGenerationParameters::DefaultGalaxyClusterRadiusInMetersAverageConstant,
                ScalarType totalAmbientMass = GalaxyClusterGenerationParameters::DefaultGalaxyClusterMassInSolarMassesAverageConstant,
                Vector3 position = Vector3{
                        GalaxyClusterGenerationParameters::DefaultGalaxyClusterUniformPositionConstant,
                        GalaxyClusterGenerationParameters::DefaultGalaxyClusterUniformPositionConstant,
                        GalaxyClusterGenerationParameters::DefaultGalaxyClusterUniformPositionConstant } ) :
                position_( position ), radius_( radius ), totalAmbientMass_( totalAmbientMass ) {
        }
        /*AbstractGalaxyCluster( AbstractGalaxyCluster&& other ) : 
                postion_( other.position_ ), radius_( other.radius_ ), 
                totalAmbientMass_( totalAmbientMass_ ), 
                galaxies_( std::move( other.galalxies ) ) {
        }*/
        Vector3 position_;
        ScalarType radius_, totalAmbientMass_;
        GalaxiesType galaxies_;
    };

    using TestingGalaxyCluster = AbstractGalaxyCluster< std::vector >;

    template< template< typename GalaxyType > typename GalaxiesContainerType >
    AbstractGalaxyCluster< GalaxiesContainerType > GenerateCluster(
            std::mt19937& generator, GalaxyClusterGenerationParameters parameters, bool singleCluster = true ) {
        ScalarType radius = parameters.GetClusterRadiusDistribution().ToDistributor()( generator );
        //Convert to kilograms.//
        ScalarType totalAmbientMass = parameters.GetClusterAmbientMassDistribution().ToDistributor()( generator ) * Constants< ScalarType >::SolarMassInKilogramsConstant;
        if( singleCluster == true ) {
            return AbstractGalaxyCluster< GalaxiesContainerType >{ radius, totalAmbientMass };
        }
        auto positionDistribution = parameters.GetClusterPositionDistribution().ToDistributor();
        return AbstractGalaxyCluster< GalaxiesContainerType >{ radius, totalAmbientMass,
                Vector3{ positionDistribution( generator ), positionDistribution( generator ), positionDistribution( generator ) } };
    }

    void GenerateGalaxies( const size_t NumberOfGalaxiesConstant, std::mt19937& generator, TestingGalaxyCluster& cluster, GalaxyGenerationParameters parameters );

    /*template< template< typename GalaxyType > typename GalaxiesContainerType = std::vector >
    struct AbstractGalaxyCluster
    {

        constexpr static ScalarType GalaxyMassInSolarMassesLowerBoundConstant = 10e11;
        constexpr static ScalarType GalaxyMassInSolarMassesUpperBoundConstant = 10e13;
        constexpr static ScalarType GalaxyMassInSolarMassesLowerBoundExponentConstant = 12;
        constexpr static ScalarType GalaxyMassInSolarMassesUpperBoundExponentConstant = 14;


        using GalaxiesType = GalaxiesContainerType< Galaxy >;
        constexpr static ScalarType DefaultGalaxyClusterTotalAmbientMassInSolarMassesConstant = 10e15;
        constexpr static ScalarType DefaultGalaxyPositionStandardDeviationConstant = 1.0;
        constexpr static Vector3 DefaultPositionConstant = Vector3{ 0.0, 0.0, 0.0 };
        constexpr static StatisticalPair DefaultGalaxyTotalMassDistributionConstant{
                .average = GalaxyMassInSolarMassesLowerBoundConstant,
                .standardDeviation = GalaxyMassInSolarMassesUpperBoundConstant };
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
    };*/
}
#endif 
