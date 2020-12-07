#include "ConstantTimeMathOperations.hpp"

namespace GravitationalLensing
{
	ScalarType Magnitude( Vector3 vector );

    template< unsigned int ConstantSetConstant = DefaultConstantSetConstant >
    struct LensRadialAcceleration
    {
        static constexpr ScalarType ConstantInnerPart = ( static_cast< ScalarType >( 3 ) *
            ( Constants< ScalarType, ConstantSetConstant >::HubbleNaughtConstant / static_cast< ScalarType >( 2 ) ) );
        static constexpr ScalarType OmegaRatioConstant = ( Constants< ScalarType, ConstantSetConstant >::OmegaSubMConstant /
            Constants< ScalarType, ConstantSetConstant >::OmegaSubLambdaConstant );
        ScalarType operator()( ScalarType time )
        {
            return MathFunctions< ScalarType >::RaiseConstant( OmegaRatioConstant, static_cast< ScalarType >( 1 ) / static_cast< ScalarType >( 3 ) ) *
                MathFunctions< ScalarType >::RaiseConstant( MathFunctions< ScalarType >::HyperbolicSineConstant(
                    ConstantInnerPart * MathFunctions< ScalarType >::SquareRootConstant( time ) ),
                    static_cast< ScalarType >( 2 ) / static_cast< ScalarType >( 3 ) );
        }
    };

    template< unsigned int ConstantSetConstant = DefaultConstantSetConstant >
    struct GravitationalPotential
    {
        const ScalarType TidalRadiusConstant;
        const ScalarType ScaledRadiusConstant;
        const ScalarType MassConstant;
        const ScalarType Term1FoldConstant = ( 1.0 + ConstantRaise( TidalRadiusConstant, 2 )( ) );
        const ScalarType Term0Constant = Constants< ScalarType, ConstantSetConstant >::GravitationalScintificNotationDigitsConstant *
            ( MassConstant / ScaledRadiusConstant );
        const ScalarType Term1Constant = ConstantRaise( TidalRadiusConstant, 2 )( ) /
            ConstantRaise( Term1FoldConstant, 2 )( );
        const ScalarType Term2Constant = ( 1.0 / TidalRadiusConstant ) - TidalRadiusConstant;
        const ScalarType Term3Constant = ( ConstantRaise( TidalRadiusConstant, 2 )( ) - 1.0 );
        const ScalarType Term4Constant = ( Constants< ScalarType, ConstantSetConstant >::PiConstant * Term3Constant ) / ( 2.0 * TidalRadiusConstant );
        const ScalarType Term5CacheConstant;
        GravitationalPotential( const ScalarType TidalRadiusConstant_, const ScalarType ScaledRadiusConstant_, const ScalarType MassConstant_ ) :
            TidalRadiusConstant( TidalRadiusConstant_ ), ScaledRadiusConstant( ScaledRadiusConstant_ ), MassConstant( MassConstant_ ),
            Term5CacheConstant( 2.0 * MathFunctions< ScalarType >::NaturalLogConstant( TidalRadiusConstant_ ) ) {
            /*std::cout << "Tidal Radius " << TidalRadiusConstant << "\n" <<
                "Scaled Radius " << ScaledRadiusConstant << "\n" <<
                "Mass Constant " << MassConstant << "\n" <<
                "Term5 Cache " << Term5CacheConstant << "\n" <<
                "Term1 Fold " << Term1FoldConstant << "\n" <<
                "Term0 " << Term0Constant << "\n" <<
                "Term1 " << Term1Constant << "\n" <<
                "Term2 " << Term2Constant << "\n" <<
                "Term3 " << Term3Constant << "\n" <<
                "Term4 " << Term4Constant << "\n";*/
        }
        ScalarType operator()( const Vector3& position )
        {
            const ScalarType MagnitudeConstant = Magnitude( position );
            const ScalarType RuntimeTerm0Constant = MathFunctions< ScalarType >::ArcTanConstant( MagnitudeConstant / TidalRadiusConstant );
            const ScalarType RuntimeTerm1Constant = ( Term2Constant - ( ( 2.0 * TidalRadiusConstant ) / MagnitudeConstant ) );
            const ScalarType RuntimeTerm2Constant = MathFunctions< ScalarType >::NaturalLogConstant( ( 1.0 + MathFunctions< ScalarType >::RaiseConstant(
                MagnitudeConstant / TidalRadiusConstant, 2.0 ) ) / MathFunctions< ScalarType >::RaiseConstant( 1.0 + MagnitudeConstant, 2 ) );
            const ScalarType RumtimeTerm3Constant = ( ( Term3Constant / ( 2.0 * MagnitudeConstant ) ) - 1.0 );
            return Term0Constant * Term1Constant * ( ( RuntimeTerm0Constant * RuntimeTerm1Constant ) +
                ( RuntimeTerm2Constant * RumtimeTerm3Constant ) + Term4Constant - Term5CacheConstant ) *
                Constants< ScalarType >::GravitationalScintificNotationCompileTimeRaiserOf10Constant;
        }
    };
}
