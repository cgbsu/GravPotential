#include "Derivatives.hpp"

#include <functional>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

namespace GravitationalLensing
{
    template< ScalarType BaseConstant, unsigned int CompileTimeRaiseerConstant >
    struct CompileTimeRaise
    {
        ScalarType constexpr operator()() {
            if constexpr( CompileTimeRaiseerConstant > 0 )
                return BaseConstant * CompileTimeRaise< BaseConstant, CompileTimeRaiseerConstant - 1 >()( );
            return 1;
        }
    };


    struct ConstantRaise
    {
        const ScalarType BaseConstant;
        const unsigned int PowerConstant;
        const ScalarType CurrentAccumulationConstant;
        ConstantRaise( const ScalarType BaseConstant_, const unsigned int PowerConstant_ ) :
            BaseConstant( BaseConstant_ ), PowerConstant( PowerConstant_ ), CurrentAccumulationConstant( 1.0 ) {}
        ConstantRaise( const ScalarType BaseConstant_, const unsigned int PowerConstant_, const ScalarType CurrentAccumulationConstant_ ) :
            BaseConstant( BaseConstant_ ), PowerConstant( PowerConstant_ ), CurrentAccumulationConstant( CurrentAccumulationConstant_ ) {}
        ScalarType constexpr operator()()
        {
            if( PowerConstant > 0 )
                return ConstantRaise( BaseConstant, PowerConstant - 1, CurrentAccumulationConstant * BaseConstant )( );
            return CurrentAccumulationConstant;
        }
    };

    constexpr unsigned int DefaultConstantSetConstant = 0;

    template< typename ScalarType = ScalarType, unsigned int ConstantSet = DefaultConstantSetConstant >
    struct Constants
    {
        /******************************************
        * Want to use a units library, there were *
        * complications when use one. *************
        ******************************************/
        constexpr static ScalarType HubbleNaughtConstant = 70.0 / 3.0856776e19; // //70.f; // kM / s Mpc
        constexpr static ScalarType OmegaSubMConstant = 0.3; // Unitless
        constexpr static ScalarType OmegaSubLambdaConstant = 1.0 - OmegaSubMConstant; // Unitless
        constexpr static ScalarType GravitationalConstant = 6.67408e-11; // m^3 / kg s^2 
        constexpr static ScalarType GravitationalScintificNotationDigitsConstant = 6.67408; // m^3 / kg s^2 
        constexpr static ScalarType GravitationalScintificNotationCompileTimeRaiserOf10Constant = 1e-11; // m^3 / kg s^2 
        constexpr static ScalarType PiConstant = ( ScalarType ) M_PI;
        constexpr static ScalarType KiloparsecInMetersConstant = 30856775814913700000.0;
        constexpr static ScalarType SolarMassInKilogramsConstant = 2e30; //kg
        constexpr static ScalarType LightYearToParsecConstant = 0.30660139; //Parsecs
    };

    template< typename ScalarType = ScalarType >
    struct MathFunctions
    {
        //Normally a cast would do the trick, that is the way that this is suppose to be done, but it wont compile with MSVC for some reason.//
        static const std::function< ScalarType( ScalarType ) > SquareRootConstant;
        static const std::function< ScalarType( ScalarType ) > HyperbolicSineConstant;
        static const std::function< ScalarType( ScalarType ) > ArcTanConstant;
        static const std::function< ScalarType( ScalarType ) > NaturalLogConstant;
        static const std::function< ScalarType( ScalarType ) > Log10Constant;
        static const std::function< ScalarType( ScalarType ) > RoundUpConstant;
        static const std::function< ScalarType( ScalarType, ScalarType ) > RaiseConstant;
        static const std::function< ScalarType( ScalarType, ScalarType ) > ModuloConstant;
    };

    //TODO: This may negativly impact performence, may want to redo this.//
    template< typename ScalarType = ScalarType >
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::SquareRootConstant =
        static_cast< ScalarType( * )( ScalarType ) >( std::sqrt );
    template< typename ScalarType = ScalarType >
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::HyperbolicSineConstant =
        static_cast< ScalarType( * )( ScalarType ) >( std::sinh );
    template< typename ScalarType = ScalarType >
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::ArcTanConstant =
        static_cast< ScalarType( * )( ScalarType ) >( std::atan );
    template< typename ScalarType = ScalarType >
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::NaturalLogConstant =
        static_cast< ScalarType( * )( ScalarType ) >( std::log );
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::Log10Constant =
        static_cast< ScalarType( * )( ScalarType ) >( std::log10 );
    const std::function< ScalarType( ScalarType ) > MathFunctions< ScalarType >::RoundUpConstant =
        static_cast< ScalarType( * )( ScalarType ) >( std::ceil );
    template< typename ScalarType = ScalarType >
    const std::function< ScalarType( ScalarType, ScalarType ) > MathFunctions< ScalarType >::RaiseConstant =
        static_cast< ScalarType( * )( ScalarType, ScalarType ) >( std::pow );
    const std::function< ScalarType( ScalarType, ScalarType ) > MathFunctions< ScalarType >::ModuloConstant =
        static_cast< ScalarType( * )( ScalarType, ScalarType ) >( std::fmod );
}
