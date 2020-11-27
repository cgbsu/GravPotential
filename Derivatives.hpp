#include "Vector.hpp"
namespace GravitationalLensing
{
    template< typename FunctionOfVectorType >
    ScalerType FirstDerivative( FunctionOfVectorType function, VectorComponent WithRespectToConstant, ScalerType IntervalConstant, const Vector3& vector )
    {
        ScalerType intervalArray[ 3 ] = { 0.0 };
        size_t WithRespectToSizeTypeConstant = static_cast< size_t >( WithRespectToConstant );
        intervalArray[ WithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalerType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ WithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalerType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        //Center Divided Difference Method.//
        //( f(x+h) - f(x-h) )/ 2h//
        return ( FirstValueConstant - SecondValueConstant ) / ( 2.0 * IntervalConstant );
    }

    template< typename FunctionOfVectorType >
    ScalerType SecondDerivative( FunctionOfVectorType function, VectorComponent WithRespectToConstant, ScalerType IntervalConstant, const Vector3& vector )
    {
        ScalerType intervalArray[ 3 ] = { 0.0 };
        size_t WithRespectToSizeTypeConstant = static_cast< size_t >( WithRespectToConstant );
        intervalArray[ WithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalerType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ WithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalerType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        //Center Divided Difference Method.//
        //( f(x+h) - 2f(x) + f(x-h) ) / h^2//
        return ( FirstValueConstant - ( 2.0 * function( vector ) ) + SecondValueConstant ) / ( IntervalConstant * IntervalConstant );
    }

    template< typename FunctionOfVector >
    ScalerType MixedDerivative( FunctionOfVector function, VectorComponent FirstWithRespectToConstant,
        VectorComponent SecondWithRespectToConstant, ScalerType IntervalConstant, const Vector3& vector )
    {
        ScalerType intervalArray[ 3 ] = { 0.0 };
        const size_t FirstWithRespectToSizeTypeConstant = static_cast< size_t >( FirstWithRespectToConstant );
        const size_t SecondWithRespectToSizeTypeConstant = static_cast< size_t >( SecondWithRespectToConstant );
        //Ugly, but if we are going to generalize it has to be done.//
        intervalArray[ FirstWithRespectToSizeTypeConstant ] = IntervalConstant;
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalerType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalerType ThirdValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ FirstWithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalerType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalerType FourthValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        //Note, could use two different "h"'s//
        // ( f(x+h,y+h) + f(x-h,y-h) - f(x+h,y-h) - f(x-h,y+h) ) / 4h ^ 2 //
        return ( FirstValueConstant + SecondValueConstant - ThirdValueConstant - FourthValueConstant ) /
            ( 4.0 * IntervalConstant * IntervalConstant );
    }

    template< typename FunctionOfVectorType >
    struct Deriver
    {
        template< VectorComponent FirstWithRespectToConstant, VectorComponent SecondWithRespectToConstant >
        struct SecondDeriver
        {
            ScalerType interval;
            Vector3 input;
            FunctionOfVectorType function;
            SecondDeriver() = default;
            SecondDeriver( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalerType Derive() {
                return MixedDerivative( function, FirstWithRespectToConstant, SecondWithRespectToConstant, interval, input );
            }
            ScalerType operator()() {
                return Derive();
            }
            SecondDeriver< FirstWithRespectToConstant, FirstWithRespectToConstant >& operator=(
                SecondDeriver< FirstWithRespectToConstant, FirstWithRespectToConstant >& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ )
            {
                function = function_;
                interval = interval_;
                input = input_;
            }
        };

        template< VectorComponent WithRespectToConstant >
        struct SecondDeriver< WithRespectToConstant, WithRespectToConstant >
        {
            ScalerType interval;
            Vector3 input;
            FunctionOfVectorType function;
            SecondDeriver() = default;
            SecondDeriver( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalerType Derive() {
                return SecondDerivative( function, WithRespectToConstant, interval, input );
            }
            ScalerType operator()() {
                return Derive();
            }
            SecondDeriver< WithRespectToConstant, WithRespectToConstant >& operator=(
                SecondDeriver< WithRespectToConstant, WithRespectToConstant >& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ )
            {
                function = function_;
                interval = interval_;
                input = input_;
            }
        };


        template< VectorComponent WithRespectToConstant >
        struct FirstDeriver
        {
            SecondDeriver< WithRespectToConstant, VectorComponent::X > x;
            SecondDeriver< WithRespectToConstant, VectorComponent::Y > y;
            SecondDeriver< WithRespectToConstant, VectorComponent::Z > z;
            ScalerType interval;
            Vector3 input;
            FunctionOfVectorType function;
            FirstDeriver() = default;
            FirstDeriver( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalerType Derive() {
                return FirstDerivative( function, WithRespectToConstant, interval, input );
            }
            ScalerType operator()() {
                return Derive();
            }
            FirstDeriver& operator=( FirstDeriver& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ )
            {
                function = function_;
                interval = interval_;
                input = input_;
                x.Assign( function, interval, input );
                y.Assign( function, interval, input );
                z.Assign( function, interval, input );
            }
        };

        ScalerType interval;
        Vector3 input;
        FunctionOfVectorType function;
        FirstDeriver< VectorComponent::X > x;
        FirstDeriver< VectorComponent::Y > y;
        FirstDeriver< VectorComponent::Z > z;
        explicit Deriver( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ ) {
            Derive( function_, interval_, input_ );
        }
        Deriver& Derive( Vector3 input_ ) {
            return Derive( function, interval, input_ );
        }
        Deriver& Derive( FunctionOfVectorType function_, ScalerType interval_, Vector3 input_ )
        {
            function = function_;
            interval = interval_;
            input = input_;
            x.Assign( function, interval, input );
            y.Assign( function, interval, input );
            z.Assign( function, interval, input );
            return *this;
        }
    };

    namespace Testing
    {
        ScalerType TestFunctionToDerive( const Vector3& v );

        ScalerType TestFunctionToCompareX( const Vector3& v );
        ScalerType TestFunctionToCompareY( const Vector3& v );
        ScalerType TestFunctionToCompareZ( const Vector3& v );

        ScalerType TestFunctionSecondDerivitiveCompare( const Vector3& v );
        ScalerType TestFunctionToMixedDerive( const Vector3& v );

        ScalerType TestFunctionMixedDerivitiveXY( const Vector3& v );
        ScalerType TestFunctionMixedDerivitiveYZ( const Vector3& v );
        ScalerType TestFunctionMixedDerivitiveXZ( const Vector3& v );

        ScalerType TestFunctionSecondDerivitiveXX( const Vector3& v );
        ScalerType TestFunctionSecondDerivitiveYY( const Vector3& v );
        ScalerType TestFunctionSecondDerivitiveZZ( const Vector3& v );

        ScalerType TestMixedFunctionFirstDerivitiveX( const Vector3& v );
        ScalerType TestMixedFunctionFirstDerivitiveY( const Vector3& v );
        ScalerType TestMixedFunctionFirstDerivitiveZ( const Vector3& v );

        void DemoDerivatives();
    }
}
