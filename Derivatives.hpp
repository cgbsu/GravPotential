#include "Vector.hpp"
namespace GravitationalLensing
{
    template< typename FunctionOfVectorType >
    ScalarType FirstDerivative( FunctionOfVectorType function, VectorComponent WithRespectToConstant, ScalarType IntervalConstant, const Vector3& vector )
    {
        ScalarType intervalArray[ 3 ] = { 0.0 };
        size_t WithRespectToSizeTypeConstant = static_cast< size_t >( WithRespectToConstant );
        intervalArray[ WithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalarType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ WithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalarType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        //Center Divided Difference Method.//
        //( f(x+h) - f(x-h) ) / 2h//
        return ( FirstValueConstant - SecondValueConstant ) / ( 2.0 * IntervalConstant );
    }

    template< typename FunctionOfVectorType >
    ScalarType SecondDerivative( FunctionOfVectorType function, VectorComponent WithRespectToConstant, ScalarType IntervalConstant, const Vector3& vector )
    {
        ScalarType intervalArray[ 3 ] = { 0.0 };
        size_t WithRespectToSizeTypeConstant = static_cast< size_t >( WithRespectToConstant );
        intervalArray[ WithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalarType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ WithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalarType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        //Center Divided Difference Method.//
        //( f(x+h) - 2f(x) + f(x-h) ) / h^2//
        return ( FirstValueConstant - ( 2.0 * function( vector ) ) + SecondValueConstant ) / ( IntervalConstant * IntervalConstant );
    }

    template< typename FunctionOfVector >
    ScalarType MixedDerivative( FunctionOfVector function, VectorComponent FirstWithRespectToConstant,
        VectorComponent SecondWithRespectToConstant, ScalarType IntervalConstant, const Vector3& vector )
    {
        ScalarType intervalArray[ 3 ] = { 0.0 };
        const size_t FirstWithRespectToSizeTypeConstant = static_cast< size_t >( FirstWithRespectToConstant );
        const size_t SecondWithRespectToSizeTypeConstant = static_cast< size_t >( SecondWithRespectToConstant );
        //Ugly, but if we are going to generalize it has to be done.//
        intervalArray[ FirstWithRespectToSizeTypeConstant ] = IntervalConstant;
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalarType FirstValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalarType ThirdValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ FirstWithRespectToSizeTypeConstant ] = -IntervalConstant;
        const ScalarType SecondValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
        intervalArray[ SecondWithRespectToSizeTypeConstant ] = IntervalConstant;
        const ScalarType FourthValueConstant = function( vector.ComponentWiseAdd( intervalArray ) );
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
            ScalarType interval;
            Vector3 input;
            FunctionOfVectorType function;
            SecondDeriver() = default;
            SecondDeriver( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalarType Derive() {
                return MixedDerivative( function, FirstWithRespectToConstant, SecondWithRespectToConstant, interval, input );
            }
            ScalarType operator()() {
                return Derive();
            }
            SecondDeriver< FirstWithRespectToConstant, FirstWithRespectToConstant >& operator=(
                SecondDeriver< FirstWithRespectToConstant, FirstWithRespectToConstant >& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ )
            {
                function = function_;
                interval = interval_;
                input = input_;
            }
        };

        template< VectorComponent WithRespectToConstant >
        struct SecondDeriver< WithRespectToConstant, WithRespectToConstant >
        {
            ScalarType interval;
            Vector3 input;
            FunctionOfVectorType function;
            SecondDeriver() = default;
            SecondDeriver( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalarType Derive() {
                return SecondDerivative( function, WithRespectToConstant, interval, input );
            }
            ScalarType operator()() {
                return Derive();
            }
            SecondDeriver< WithRespectToConstant, WithRespectToConstant >& operator=(
                SecondDeriver< WithRespectToConstant, WithRespectToConstant >& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ )
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
            ScalarType interval;
            Vector3 input;
            FunctionOfVectorType function;
            FirstDeriver() = default;
            FirstDeriver( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ ) :
                function( function_ ), interval( interval_ ), input( input_ ) {
            }
            ScalarType Derive() {
                return FirstDerivative( function, WithRespectToConstant, interval, input );
            }
            ScalarType operator()() {
                return Derive();
            }
            FirstDeriver& operator=( FirstDeriver& other ) {
                Assign( other.function, other.interval, other.input );
                return *this;
            }
            void Assign( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ )
            {
                function = function_;
                interval = interval_;
                input = input_;
                x.Assign( function, interval, input );
                y.Assign( function, interval, input );
                z.Assign( function, interval, input );
            }
        };

        ScalarType interval;
        Vector3 input;
        FunctionOfVectorType function;
        FirstDeriver< VectorComponent::X > x;
        FirstDeriver< VectorComponent::Y > y;
        FirstDeriver< VectorComponent::Z > z;
        explicit Deriver( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ ) {
            Derive( function_, interval_, input_ );
        }
        Deriver& Derive( Vector3 input_ ) {
            return Derive( function, interval, input_ );
        }
        Deriver& Derive( FunctionOfVectorType function_, ScalarType interval_, Vector3 input_ )
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
        ScalarType TestFunctionToDerive( const Vector3& v );

        ScalarType TestFunctionToCompareX( const Vector3& v );
        ScalarType TestFunctionToCompareY( const Vector3& v );
        ScalarType TestFunctionToCompareZ( const Vector3& v );

        ScalarType TestFunctionSecondDerivitiveCompare( const Vector3& v );
        ScalarType TestFunctionToMixedDerive( const Vector3& v );

        ScalarType TestFunctionMixedDerivitiveXY( const Vector3& v );
        ScalarType TestFunctionMixedDerivitiveYZ( const Vector3& v );
        ScalarType TestFunctionMixedDerivitiveXZ( const Vector3& v );

        ScalarType TestFunctionSecondDerivitiveXX( const Vector3& v );
        ScalarType TestFunctionSecondDerivitiveYY( const Vector3& v );
        ScalarType TestFunctionSecondDerivitiveZZ( const Vector3& v );

        ScalarType TestMixedFunctionFirstDerivitiveX( const Vector3& v );
        ScalarType TestMixedFunctionFirstDerivitiveY( const Vector3& v );
        ScalarType TestMixedFunctionFirstDerivitiveZ( const Vector3& v );

        void DemoDerivatives();
    }
}
