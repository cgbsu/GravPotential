#include "GalaxyDistribution.hpp"
namespace GravitationalLensing
{
    void GenerateGalaxies( const size_t NumberOfGalaxiesConstant, std::mt19937& generator, TestingGalaxyCluster& cluster, GalaxyGenerationParameters parameters )
    {
        std::normal_distribution< double > galaxyPositionDistributor( 0.0, parameters.GetGalaxyPositionStandardDeviation() );
        auto radiusDistributor = parameters.GetGalaxyRadiusDistribution().ToDistributor();
        const ScalarType MaximumGalaxyMassExponentConstant = MathFunctions< ScalarType >::Log10Constant( parameters.GetGalaxyMassInSolarMassesUpperBound() );
        const ScalarType MinimumGalaxyMassExponentConstant = MathFunctions< ScalarType >::Log10Constant( parameters.GetGalaxyMassInSolarMassesLowerBound() );
        const ScalarType MassExponentStandardDeviationConstant = MaximumGalaxyMassExponentConstant - MinimumGalaxyMassExponentConstant;
        std::normal_distribution< double > galaxyMassExponentDistribution(
            MaximumGalaxyMassExponentConstant - MassExponentStandardDeviationConstant,
            MassExponentStandardDeviationConstant );
        std::normal_distribution< double > significand( 10.0 );
        for( size_t i = 0; i < NumberOfGalaxiesConstant; ++i )
        {
            ScalarType massExponent = MinimumGalaxyMassExponentConstant - 1.0;
            while( massExponent < MinimumGalaxyMassExponentConstant ) {
                massExponent = MathFunctions< ScalarType >::ModuloConstant(
                    MathFunctions< ScalarType >::RoundUpConstant( galaxyMassExponentDistribution( generator ) ),
                    MaximumGalaxyMassExponentConstant + 1 );
            }
            cluster.galaxies_.push_back(
                    Galaxy{
                            .position_ = Vector3{
                                    cluster.position_.x_ + galaxyPositionDistributor( generator ),
                                    cluster.position_.y_ + galaxyPositionDistributor( generator ),
                                    cluster.position_.z_ + galaxyPositionDistributor( generator ) },
                            .totalMass_ = MathFunctions< ScalarType >::ModuloConstant( significand( generator ), 10.0 ) * 
                                    MathFunctions< ScalarType >::RaiseConstant( 10.0, massExponent ),
                            .radius_ = radiusDistributor( generator )
                    } );
        }
    }


    VectorDistributor MakeDistributor( ScalarType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator )
    {
        return VectorDistributor{ .x_ = std::normal_distribution< ScalarType >{ bounds.x_ / 2.0, ( bounds.x_ / 2.0 ) * PlacementStandardDeviationConstant },
                .y_ = std::normal_distribution< ScalarType >{ bounds.y_ / 2.0, ( bounds.y_ / 2.0 ) * PlacementStandardDeviationConstant },
                .z_ = std::normal_distribution< ScalarType >{ bounds.z_ / 2.0, ( bounds.z_ / 2.0 ) * PlacementStandardDeviationConstant } };
    }

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator )
    {
        return Vector3{ static_cast< std::normal_distribution< ScalarType > >( distributor.x_ )( generator ),
                static_cast< std::normal_distribution< ScalarType > >( distributor.y_ )( generator ),
                static_cast< std::normal_distribution< ScalarType > >( distributor.z_ )( generator ) };
    }

/*  For refrence  
    void PlaceGalaxies( StatisticalPair NumberOfGalaxiesConstant, ScalarType GalaxyPositionStandardDeviationCoefficentConstant,
        StatisticalPair GalaxyDarkMatterCoefficentConstant, GalaxyCluster& galaxyCluster, std::mt19937& generator )
    {
        auto distributor = MakeDistributor( GalaxyPositionStandardDeviationCoefficentConstant, galaxyCluster.dimensions_, generator );
        std::normal_distribution darkMatterDistribution(
            GalaxyDarkMatterCoefficentConstant.average, GalaxyDarkMatterCoefficentConstant.standardDeviation );
        const size_t AmountOfGalaxiesConstant = ( size_t ) std::normal_distribution(
            NumberOfGalaxiesConstant.average, NumberOfGalaxiesConstant.standardDeviation )( generator );
        for( size_t i = 0; i < AmountOfGalaxiesConstant; ++i )
        {
            Galaxy newGalaxy;
            newGalaxy.darkMatterCoefficent_ = std::abs( PercentToDecimalConstant * darkMatterDistribution( generator ) );
            newGalaxy.position_ = RandomPointWithinBounds( distributor, generator );
            galaxyCluster.galaxies_.push_back( newGalaxy );
        }
    }

    GalaxyCluster GenerateGalaxyCluster( StatisticalPair DarkMatterCoeffiecentsConstant,
        StatisticalPair GalaxyClusterDimensionsConstant, StatisticalPair GalaxyClusterPositionConstant, std::mt19937& generator )
    {
        //Selected to sound somewhat cool.//
        std::normal_distribution< ScalarType > dimensionalGenerator( GalaxyClusterDimensionsConstant.average, GalaxyClusterDimensionsConstant.standardDeviation );
        std::normal_distribution< ScalarType > darkMatterCoefficentGenerator( DarkMatterCoeffiecentsConstant.average, DarkMatterCoeffiecentsConstant.standardDeviation );
        std::normal_distribution< ScalarType > positionGenrator( GalaxyClusterPositionConstant.average, GalaxyClusterPositionConstant.standardDeviation );
        Vector3 position{ positionGenrator( generator ), positionGenrator( generator ),
                positionGenrator( generator ) };
        Vector3 dimensions{ std::abs( dimensionalGenerator( generator ) ),
                std::abs( dimensionalGenerator( generator ) ), std::abs( dimensionalGenerator( generator ) ) };
        return GalaxyCluster{ .position_ = position, .dimensions_ = dimensions, .darkMatterCoefficent_ = std::abs( 
                PercentToDecimalConstant * darkMatterCoefficentGenerator( generator ) ), .galaxies_ = GalaxyCluster::GalaxiesType() };
    }
    */
}
