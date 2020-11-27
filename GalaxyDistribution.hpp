#include "GravitationalFunctions.hpp"

namespace GravitationalLensing
{
    struct Galaxy {
        Vector3 position_;
        ScalerType darkMatterCoefficent_;
    };

    struct GalaxyCluster
    {
        using GalaxiesType = std::vector< Galaxy >;
        Vector3 position_, dimensions_;
        ScalerType darkMatterCoefficent_;
        GalaxiesType galaxies_;
    };

    struct VectorDistributor {
        std::normal_distribution< ScalerType > x_, y_, z_;
    };

    template< typename StatisticalType >
    struct AbstractStatisticalPair {
        StatisticalType average;
        StatisticalType standardDeviation;
    };

    using StatisticalPair = const AbstractStatisticalPair< const ScalerType >;

    VectorDistributor MakeDistributor( ScalerType PlacementStandardDeviationConstant, const Vector3& bounds, std::mt19937& generator );

    Vector3 RandomPointWithinBounds( const VectorDistributor& distributor, std::mt19937& generator );

    void PlaceGalaxies( StatisticalPair NumberOfGalaxiesConstant, ScalerType GalaxyPositionStandardDeviationCoefficentConstant,
        StatisticalPair GalaxyDarkMatterCoefficentConstant, GalaxyCluster& galaxyCluster, std::mt19937& generator );

    using GalaxyClusterArrayType = std::vector< GalaxyCluster >;

    GalaxyCluster GenerateGalaxyCluster( StatisticalPair DarkMatterCoeffiecentsConstant,
        StatisticalPair GalaxyClusterDimensionsConstant, StatisticalPair GalaxyClusterPositionConstant, std::mt19937& generator );
}
