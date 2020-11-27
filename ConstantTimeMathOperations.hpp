#include "Derivatives.hpp"

#include <functional>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>


namespace GravitationalLensing
{
    template< ScalerType BaseConstant, unsigned int CompileTimeRaiseerConstant >
    struct CompileTimeRaise
    {
        ScalerType constexpr operator()() {
            if constexpr( CompileTimeRaiseerConstant > 0 )
                return BaseConstant * CompileTimeRaise< BaseConstant, CompileTimeRaiseerConstant - 1 >()( );
            return 1;
        }
    };


    struct ConstantRaise
    {
        const ScalerType BaseConstant;
        const unsigned int PowerConstant;
        const ScalerType CurrentAccumulationConstant;
        ConstantRaise( const ScalerType BaseConstant_, const unsigned int PowerConstant_ ) :
            BaseConstant( BaseConstant_ ), PowerConstant( PowerConstant_ ), CurrentAccumulationConstant( 1.0 ) {}
        ConstantRaise( const ScalerType BaseConstant_, const unsigned int PowerConstant_, const ScalerType CurrentAccumulationConstant_ ) :
            BaseConstant( BaseConstant_ ), PowerConstant( PowerConstant_ ), CurrentAccumulationConstant( CurrentAccumulationConstant_ ) {}
        ScalerType constexpr operator()()
        {
            if( PowerConstant > 0 )
                return ConstantRaise( BaseConstant, PowerConstant - 1, CurrentAccumulationConstant * BaseConstant )( );
            return CurrentAccumulationConstant;
        }
    };

    constexpr unsigned int DefaultConstantSetConstant = 0;

    template< typename ScalerType = ScalerType, unsigned int ConstantSet = DefaultConstantSetConstant >
    struct Constants
    {
        /******************************************
        * Want to use a units library, there were *
        * complications when use one. *************
        ******************************************/
        constexpr static ScalerType HubbleNaughtConstant = 70.0 / 3.0856776e19; // //70.f; // kM / s Mpc
        constexpr static ScalerType OmegaSubMConstant = 0.3; // Unitless
        constexpr static ScalerType OmegaSubLambdaConstant = 1.0 - OmegaSubMConstant; // Unitless
        constexpr static ScalerType GravitationalConstant = 6.67408e-11; // m^3 / kg s^2 
        constexpr static ScalerType GravitationalScintificNotationDigitsConstant = 6.67408; // m^3 / kg s^2 
        constexpr static ScalerType GravitationalScintificNotationCompileTimeRaiserOf10Constant = 1e-11; // m^3 / kg s^2 
        constexpr static ScalerType PiConstant = ( ScalerType ) M_PI;
        constexpr static ScalerType KiloparsecInMetersConstant = 30856775814913700000.0;
    };

    template< typename ScalerType = ScalerType >
    struct MathFunctions
    {
        //Normally a cast would do the trick, that is the way that this is suppose to be done, but it wont compile with MSVC for some reason.//
        static const std::function< ScalerType( ScalerType ) > SquareRootConstant;
        static const std::function< ScalerType( ScalerType ) > HyperbolicSineConstant;
        static const std::function< ScalerType( ScalerType ) > ArcTanConstant;
        static const std::function< ScalerType( ScalerType ) > NaturalLogConstant;
        static const std::function< ScalerType( ScalerType, ScalerType ) > RaiseConstant;
    };
    template< typename ScalerType = ScalerType >
    const std::function< ScalerType( ScalerType ) > MathFunctions< ScalerType >::SquareRootConstant =
        static_cast< ScalerType( * )( ScalerType ) >( std::sqrt );
    template< typename ScalerType = ScalerType >
    const std::function< ScalerType( ScalerType ) > MathFunctions< ScalerType >::HyperbolicSineConstant =
        static_cast< ScalerType( * )( ScalerType ) >( std::sinh );
    template< typename ScalerType = ScalerType >
    const std::function< ScalerType( ScalerType ) > MathFunctions< ScalerType >::ArcTanConstant =
        static_cast< ScalerType( * )( ScalerType ) >( std::atan );
    template< typename ScalerType = ScalerType >
    const std::function< ScalerType( ScalerType ) > MathFunctions< ScalerType >::NaturalLogConstant =
        static_cast< ScalerType( * )( ScalerType ) >( std::log );
    template< typename ScalerType = ScalerType >
    const std::function< ScalerType( ScalerType, ScalerType ) > MathFunctions< ScalerType >::RaiseConstant =
        static_cast< ScalerType( * )( ScalerType, ScalerType ) >( std::pow );
}
