#include "Derivatives.hpp"

#include <iostream>

namespace GravitationalLensing
{
    namespace Testing
    {
        ScalarType TestFunctionToDerive( const Vector3& v ) {
            return ( v.x_ * v.x_ ) + ( v.y_ * v.y_ ) + ( v.z_ * v.z_ );
        }

        ScalarType TestFunctionToCompareX( const Vector3& v ) {
            return ( 2.0 * v.x_ );
        }

        ScalarType TestFunctionToCompareY( const Vector3& v ) {
            return ( 2.0 * v.y_ );
        }

        ScalarType TestFunctionToCompareZ( const Vector3& v ) {
            return ( 2.0 * v.z_ );
        }

        ScalarType TestFunctionSecondDerivitiveCompare( const Vector3& v ) {
            return 2.0;
        }

        ScalarType TestFunctionToMixedDerive( const Vector3& v ) {
            return v.x_ * v.x_ * v.y_ * v.y_ * v.z_ * v.z_;
        }

        ScalarType TestFunctionMixedDerivitiveXY( const Vector3& v ) {
            return  4.0 * v.x_ * v.y_ * v.z_ * v.z_;
        }

        ScalarType TestFunctionMixedDerivitiveYZ( const Vector3& v ) {
            return  4.0 * v.z_ * v.y_ * v.x_ * v.x_;
        }

        ScalarType TestFunctionMixedDerivitiveXZ( const Vector3& v ) {
            return  4.0 * v.x_ * v.z_ * v.y_ * v.y_;
        }

        ///////////////////////////////////////////////////////////

        ScalarType TestFunctionSecondDerivitiveXX( const Vector3& v ) {
            return  2.0 * v.y_ * v.y_ * v.z_ * v.z_;
        }
        ScalarType TestFunctionSecondDerivitiveYY( const Vector3& v ) {
            return  2.0 * v.x_ * v.x_ * v.z_ * v.z_;
        }
        ScalarType TestFunctionSecondDerivitiveZZ( const Vector3& v ) {
            return  2.0 * v.y_ * v.y_ * v.x_ * v.x_;
        }

        ////////////////////////////////////////////////////////////

        ScalarType TestMixedFunctionFirstDerivitiveX( const Vector3& v ) {
            return  2.0 * v.y_ * v.y_ * v.z_ * v.z_ * v.x_;
        }
        ScalarType TestMixedFunctionFirstDerivitiveY( const Vector3& v ) {
            return  2.0 * v.x_ * v.x_ * v.z_ * v.z_ * v.y_;
        }
        ScalarType TestMixedFunctionFirstDerivitiveZ( const Vector3& v ) {
            return  2.0 * v.y_ * v.y_ * v.x_ * v.x_ * v.z_;
        }

        void DemoDerivatives()
        {
            Vector3 v{ 23.0, 43.0, 64.0 };
            ScalarType theory = FirstDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToDerive, VectorComponent::X, .1, v );
            ScalarType actual = TestFunctionToCompareX( v );
            std::cout << "First Derivative Function X " << theory << "\n";
            std::cout << "First Derivative Hand Written Function X " << actual << "\n";
            theory = FirstDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToDerive, VectorComponent::Y, .1, v );
            actual = TestFunctionToCompareY( v );
            std::cout << "First Derivative Function Y " << theory << "\n";
            std::cout << "First Derivative Hand Written Function Y " << actual << "\n";
            theory = FirstDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToDerive, VectorComponent::Z, .1, v );
            actual = TestFunctionToCompareZ( v );
            std::cout << "First Derivative Function Z " << theory << "\n";
            std::cout << "First Derivative Hand Written Function Z " << actual << "\n";
            theory = SecondDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToDerive, VectorComponent::X, .1, v );
            actual = TestFunctionSecondDerivitiveCompare( v );
            std::cout << "Second Derivative Function " << theory << "\n";
            std::cout << "Second Derivative Hand Written Function " << actual << "\n";
            theory = MixedDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToMixedDerive, VectorComponent::X, VectorComponent::Y, .1, v );
            actual = TestFunctionMixedDerivitiveXY( v );
            std::cout << "Mixed Derivative XY: " << theory << "\n";
            std::cout << "Mixed Derivative XY Hand Writen: " << actual << "\n";
            theory = MixedDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToMixedDerive, VectorComponent::X, VectorComponent::Z, .1, v );
            actual = TestFunctionMixedDerivitiveXZ( v );
            std::cout << "Mixed Derivative XZ: " << theory << "\n";
            std::cout << "Mixed Derivative XZ Hand Writen: " << actual << "\n";
            theory = MixedDerivative< ScalarType( * )( const Vector3& ) >( &TestFunctionToMixedDerive, VectorComponent::Y, VectorComponent::Z, .1, v );
            actual = TestFunctionMixedDerivitiveYZ( v );
            std::cout << "Mixed Derivative YZ: " << theory << "\n";
            std::cout << "Mixed Derivative YZ Hand Writen: " << actual << "\n";
            std::cout << "!!!!!!!!!!!!!!TESTING DERIVER!!!!!!!!!!!!!!\n";
            Deriver< ScalarType( * )( const Vector3& ) > deriver0( &TestFunctionToDerive, .1, v );
            std::cout << "<First Derivatives>\n";
            theory = deriver0.x();
            actual = TestFunctionToCompareX( v );
            std::cout << "Deriver First Derivative Function X " << theory << "\n";
            std::cout << "First Derivative Hand Written Function X " << actual << "\n";
            theory = deriver0.y();
            actual = TestFunctionToCompareY( v );
            std::cout << "Deriver First Derivative Function Y " << theory << "\n";
            std::cout << "First Derivative Hand Written Function Y " << actual << "\n";
            theory = deriver0.z();
            actual = TestFunctionToCompareZ( v );
            std::cout << "Deriver First Derivative Function Z " << theory << "\n";
            std::cout << "First Derivative Hand Written Function Z " << actual << "\n";
            std::cout << "<Second Derivatives>\n";
            theory = deriver0.x.x();
            actual = TestFunctionSecondDerivitiveCompare( v );
            std::cout << "Deriver First Derivative Function XX " << theory << "\n";
            std::cout << "First Derivative Hand Written Function XX " << actual << "\n";
            theory = deriver0.y.y();
            actual = TestFunctionSecondDerivitiveCompare( v );
            std::cout << "Deriver First Derivative Function YY " << theory << "\n";
            std::cout << "First Derivative Hand Written Function YY " << actual << "\n";
            theory = deriver0.z.z();
            actual = TestFunctionSecondDerivitiveCompare( v );
            std::cout << "Deriver First Derivative Function ZZ " << theory << "\n";
            std::cout << "First Derivative Hand Written Function ZZ " << actual << "\n";
            std::cout << "Testing with input as mixed function\n";
            Deriver< ScalarType( * )( const Vector3& ) > deriver1( &TestFunctionToMixedDerive, .1, v );
            theory = deriver1.x();
            actual = TestMixedFunctionFirstDerivitiveX( v );
            std::cout << "Deriver First Derivative X " << theory << "\n";
            std::cout << "Hand Written First Derivative X " << actual << "\n";
            theory = deriver1.y();
            actual = TestMixedFunctionFirstDerivitiveY( v );
            std::cout << "Deriver First Derivative Y " << theory << "\n";
            std::cout << "Hand Written First Derivative Y " << actual << "\n";
            theory = deriver1.z();
            actual = TestMixedFunctionFirstDerivitiveZ( v );
            std::cout << "Deriver First Derivative Z " << theory << "\n";
            std::cout << "Hand Written First Derivative Z " << actual << "\n";
            theory = deriver1.x.x();
            actual = TestFunctionSecondDerivitiveXX( v );
            std::cout << "Deriver First Derivative XX " << theory << "\n";
            std::cout << "Hand Written First Derivative XX " << actual << "\n";
            theory = deriver1.y.y();
            actual = TestFunctionSecondDerivitiveYY( v );
            std::cout << "Deriver First Derivative YY " << theory << "\n";
            std::cout << "Hand Written First Derivative YY " << actual << "\n";
            theory = deriver1.z.z();
            actual = TestFunctionSecondDerivitiveZZ( v );
            std::cout << "Deriver First Derivative ZZ " << theory << "\n";
            std::cout << "Hand Written First Derivative ZZ " << actual << "\n";
            theory = deriver1.x.y();
            actual = TestFunctionMixedDerivitiveXY( v );
            std::cout << "Deriver First Derivative XY " << theory << "\n";
            std::cout << "Hand Written First Derivative XY " << actual << "\n";
            theory = deriver1.y.x();
            actual = TestFunctionMixedDerivitiveXY( v );
            std::cout << "Deriver First Derivative XY " << theory << "\n";
            std::cout << "Hand Written First Derivative XY " << actual << "\n";
            theory = deriver1.x.z();
            actual = TestFunctionMixedDerivitiveXZ( v );
            std::cout << "Deriver First Derivative XZ " << theory << "\n";
            std::cout << "Hand Written First Derivative XZ " << actual << "\n";
            theory = deriver1.z.x();
            actual = TestFunctionMixedDerivitiveXZ( v );
            std::cout << "Deriver First Derivative XZ " << theory << "\n";
            std::cout << "Hand Written First Derivative ZX " << actual << "\n";
            theory = deriver1.z.y();
            actual = TestFunctionMixedDerivitiveYZ( v );
            std::cout << "Deriver First Derivative ZY " << theory << "\n";
            std::cout << "Hand Written First Derivative ZY " << actual << "\n";
            theory = deriver1.y.z();
            actual = TestFunctionMixedDerivitiveYZ( v );
            std::cout << "Deriver First Derivative YZ " << theory << "\n";
            std::cout << "Hand Written First Derivative YZ " << actual << "\n";
        }
    }
}
