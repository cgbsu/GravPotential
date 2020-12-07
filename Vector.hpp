namespace GravitationalLensing
{
    using ScalarType = long double;

    enum class VectorComponent : size_t
    {
        X = 0,
        Y = 1,
        Z = 2
    };

    struct alignas( ScalarType ) Vector3
    {
        ScalarType x_, y_, z_;
        Vector3( ScalarType x = 0.0, ScalarType y = 0.0, ScalarType z = 0.0 ) : x_( x ), y_( y ), z_( z ) {
        }
        Vector3& SetComponents( ScalarType x, ScalarType y, ScalarType z )
        {
            x_ = x;
            y_ = y;
            z_ = z;
            return *this;
        }
        Vector3 ComponentWiseAdd( ScalarType xDifference, ScalarType yDifference, ScalarType zDifference ) const {
            return Vector3{ x_ + xDifference, y_ + yDifference, z_ + zDifference };
        }
        Vector3 ComponentWiseAdd( ScalarType difference[ 3 ] ) const {
            return Vector3{ x_ + difference[ 0 ], y_ + difference[ 1 ], z_ + difference[ 2 ] };
        }
        Vector3& ComponentWiseAddEquals( ScalarType xDifference, ScalarType yDifference, ScalarType zDifference )
        {
            x_ += xDifference;
            y_ += yDifference;
            z_ += zDifference;
            return *this;
        }
        Vector3 operator+( Vector3& other ) {
            return ComponentWiseAdd( other.x_, other.y_, other.z_ );
        }
        Vector3& operator+=( Vector3& other ) {
            return ComponentWiseAddEquals( other.x_, other.y_, other.z_ );
        }
        Vector3 operator-() {
            return Vector3( -x_, -y_, -z_ );
        }
        ScalarType operator[]( VectorComponent index )
        {
            switch( index )
            {
            case VectorComponent::X:
                return x_;
            case VectorComponent::Y:
                return y_;
            case VectorComponent::Z:
                return z_;
            default:
                break;
            }
            return 0;
        }
    };

    #define SPECIFY_COMPONENT_SPECILIZATION_MACRO( COMPONENT_PARAMETER, VECTOR_COMPONENT_PARAMETER ) \
        template<> \
        struct SpecifyComponent< VectorComponent:: VECTOR_COMPONENT_PARAMETER > \
        { \
            const ScalarType& ComponentConstant; \
            SpecifyComponent( const Vector3& vector ) : \
                    ComponentConstant( vector. COMPONENT_PARAMETER ) { \
            } \
            operator const ScalarType&() { \
                return ComponentConstant; \
            } \
        };

    template< VectorComponent ComponentSelectorConstant >
    struct SpecifyComponent {
        SpecifyComponent( const Vector3& vector ) = delete;
        operator const ScalarType& ( ) = delete;
    };

    SPECIFY_COMPONENT_SPECILIZATION_MACRO( x_, X )
    SPECIFY_COMPONENT_SPECILIZATION_MACRO( y_, Y )
    SPECIFY_COMPONENT_SPECILIZATION_MACRO( z_, Z )
}
