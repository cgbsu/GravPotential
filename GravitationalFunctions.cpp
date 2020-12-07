#include "GravitationalFunctions.hpp"

namespace GravitationalLensing
{
    ScalarType Magnitude( Vector3 vector ) {
        return MathFunctions< ScalarType >::SquareRootConstant( ( vector.x_ * vector.x_ ) +
                ( vector.y_ * vector.y_ ) + ( vector.z_ * vector.z_ ) );
    }
}
