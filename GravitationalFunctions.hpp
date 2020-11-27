#include "ConstantTimeMathOperations.hpp"

namespace GravitationalLensing
{
	ScalerType Magnitude( Vector3 vector );

    template< unsigned int ConstantSetConstant = DefaultConstantSetConstant >
    struct LensRadialAcceleration
    {
        static constexpr ScalerType ConstantInnerPart = ( static_cast< ScalerType >( 3 ) *
            ( Constants< ScalerType, ConstantSetConstant >::HubbleNaughtConstant / static_cast< ScalerType >( 2 ) ) );
        static constexpr ScalerType OmegaRatioConstant = ( Constants< ScalerType, ConstantSetConstant >::OmegaSubMConstant /
            Constants< ScalerType, ConstantSetConstant >::OmegaSubLambdaConstant );
        ScalerType operator()( ScalerType time )
        {
            return MathFunctions< ScalerType >::RaiseConstant( OmegaRatioConstant, static_cast< ScalerType >( 1 ) / static_cast< ScalerType >( 3 ) ) *
                MathFunctions< ScalerType >::RaiseConstant( MathFunctions< ScalerType >::HyperbolicSineConstant(
                    ConstantInnerPart * MathFunctions< ScalerType >::SquareRootConstant( time ) ),
                    static_cast< ScalerType >( 2 ) / static_cast< ScalerType >( 3 ) );
        }
    };

    template< unsigned int ConstantSetConstant = DefaultConstantSetConstant >
    struct GravitationalPotential
    {
        const ScalerType TidalRadiusConstant;
        const ScalerType ScaledRadiusConstant;
        const ScalerType MassConstant;
        const ScalerType Term1FoldConstant = ( 1.0 + ConstantRaise( TidalRadiusConstant, 2 )( ) );
        const ScalerType Term0Constant = Constants< ScalerType, ConstantSetConstant >::GravitationalScintificNotationDigitsConstant *
            ( MassConstant / ScaledRadiusConstant );
        const ScalerType Term1Constant = ConstantRaise( TidalRadiusConstant, 2 )( ) /
            ConstantRaise( Term1FoldConstant, 2 )( );
        const ScalerType Term2Constant = ( 1.0 / TidalRadiusConstant ) - TidalRadiusConstant;
        const ScalerType Term3Constant = ( ConstantRaise( TidalRadiusConstant, 2 )( ) - 1.0 );
        const ScalerType Term4Constant = ( Constants< ScalerType, ConstantSetConstant >::PiConstant * Term3Constant ) / ( 2.0 * TidalRadiusConstant );
        const ScalerType Term5CacheConstant;
        GravitationalPotential( const ScalerType TidalRadiusConstant_, const ScalerType ScaledRadiusConstant_, const ScalerType MassConstant_ ) :
            TidalRadiusConstant( TidalRadiusConstant_ ), ScaledRadiusConstant( ScaledRadiusConstant_ ), MassConstant( MassConstant_ ),
            Term5CacheConstant( 2.0 * MathFunctions< ScalerType >::NaturalLogConstant( TidalRadiusConstant_ ) ) {
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
        ScalerType operator()( const Vector3& position )
        {
            const ScalerType MagnitudeConstant = Magnitude( position );
            const ScalerType RuntimeTerm0Constant = MathFunctions< ScalerType >::ArcTanConstant( MagnitudeConstant / TidalRadiusConstant );
            const ScalerType RuntimeTerm1Constant = ( Term2Constant - ( ( 2.0 * TidalRadiusConstant ) / MagnitudeConstant ) );
            const ScalerType RuntimeTerm2Constant = MathFunctions< ScalerType >::NaturalLogConstant( ( 1.0 + MathFunctions< ScalerType >::RaiseConstant(
                MagnitudeConstant / TidalRadiusConstant, 2.0 ) ) / MathFunctions< ScalerType >::RaiseConstant( 1.0 + MagnitudeConstant, 2 ) );
            const ScalerType RumtimeTerm3Constant = ( ( Term3Constant / ( 2.0 * MagnitudeConstant ) ) - 1.0 );
            return Term0Constant * Term1Constant * ( ( RuntimeTerm0Constant * RuntimeTerm1Constant ) +
                ( RuntimeTerm2Constant * RumtimeTerm3Constant ) + Term4Constant - Term5CacheConstant ) *
                Constants< ScalerType >::GravitationalScintificNotationCompileTimeRaiserOf10Constant;
        }
    };
}
