#include "GravitationalFunctions.hpp"

namespace GravitationalLensing
{
    ScalerType Magnitude( Vector3 vector ) {
        return MathFunctions< ScalerType >::SquareRootConstant( ( vector.x_ * vector.x_ ) +
                ( vector.y_ * vector.y_ ) + ( vector.z_ * vector.z_ ) );
    }
}
