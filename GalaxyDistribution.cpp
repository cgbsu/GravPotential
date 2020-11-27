#include "GalaxyDistribution.hpp"

namespace GravitationalLensing
{
    VectorDistributor MakeDistributor( ScalerType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator )
    {
        return VectorDistributor{ .x_ = std::normal_distribution< ScalerType >{ bounds.x_ / 2.0, ( bounds.x_ / 2.0 ) * PlacementStandardDeviationConstant },
                .y_ = std::normal_distribution< ScalerType >{ bounds.y_ / 2.0, ( bounds.y_ / 2.0 ) * PlacementStandardDeviationConstant },
                .z_ = std::normal_distribution< ScalerType >{ bounds.z_ / 2.0, ( bounds.z_ / 2.0 ) * PlacementStandardDeviationConstant } };
    }

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator )
    {
        return Vector3{ static_cast< std::normal_distribution< ScalerType > >( distributor.x_ )( generator ),
                static_cast< std::normal_distribution< ScalerType > >( distributor.y_ )( generator ),
                static_cast< std::normal_distribution< ScalerType > >( distributor.z_ )( generator ) };
    }

    void PlaceGalaxies( StatisticalPair NumberOfGalaxiesConstant, ScalerType GalaxyPositionStandardDeviationCoefficentConstant,
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
            newGalaxy.darkMatterCoefficent_ = std::abs( .01 * darkMatterDistribution( generator ) );
            newGalaxy.position_ = RandomPointWithinBounds( distributor, generator );
            galaxyCluster.galaxies_.push_back( newGalaxy );
        }
    }

    GalaxyCluster GenerateGalaxyCluster( StatisticalPair DarkMatterCoeffiecentsConstant,
        StatisticalPair GalaxyClusterDimensionsConstant, StatisticalPair GalaxyClusterPositionConstant, std::mt19937& generator )
    {
        //Selected to sound somewhat cool.//
        std::normal_distribution< ScalerType > dimensionalGenerator( GalaxyClusterDimensionsConstant.average, GalaxyClusterDimensionsConstant.standardDeviation );
        std::normal_distribution< ScalerType > darkMatterCoefficentGenerator( DarkMatterCoeffiecentsConstant.average, DarkMatterCoeffiecentsConstant.standardDeviation );
        std::normal_distribution< ScalerType > positionGenrator( GalaxyClusterPositionConstant.average, GalaxyClusterPositionConstant.standardDeviation );
        Vector3 position{ positionGenrator( generator ), positionGenrator( generator ),
                positionGenrator( generator ) };
        Vector3 dimensions{ std::abs( dimensionalGenerator( generator ) ),
                std::abs( dimensionalGenerator( generator ) ), std::abs( dimensionalGenerator( generator ) ) };
        return GalaxyCluster{ .position_ = position, .dimensions_ = dimensions, .darkMatterCoefficent_ = std::abs( .01 * darkMatterCoefficentGenerator( generator ) ),
                .galaxies_ = GalaxyCluster::GalaxiesType() };
    }
}
