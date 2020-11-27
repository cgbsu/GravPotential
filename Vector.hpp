namespace GravitationalLensing
{
    using ScalerType = long double;

    enum class VectorComponent : size_t
    {
        X = 0,
        Y = 1,
        Z = 2
    };

    struct alignas( ScalerType ) Vector3
    {
        ScalerType x_, y_, z_;
        Vector3( ScalerType x = 0.0, ScalerType y = 0.0, ScalerType z = 0.0 ) : x_( x ), y_( y ), z_( z ) {
        }
        Vector3& SetComponents( ScalerType x, ScalerType y, ScalerType z )
        {
            x_ = x;
            y_ = y;
            z_ = z;
            return *this;
        }
        Vector3 ComponentWiseAdd( ScalerType xDifference, ScalerType yDifference, ScalerType zDifference ) const {
            return Vector3{ x_ + xDifference, y_ + yDifference, z_ + zDifference };
        }
        Vector3 ComponentWiseAdd( ScalerType difference[ 3 ] ) const {
            return Vector3{ x_ + difference[ 0 ], y_ + difference[ 1 ], z_ + difference[ 2 ] };
        }
        Vector3& ComponentWiseAddEquals( ScalerType xDifference, ScalerType yDifference, ScalerType zDifference )
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
        ScalerType operator[]( VectorComponent index )
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
            const ScalerType& ComponentConstant; \
            SpecifyComponent( const Vector3& vector ) : \
                    ComponentConstant( vector. COMPONENT_PARAMETER ) { \
            } \
            operator const ScalerType&() { \
                return ComponentConstant; \
            } \
        };

    template< VectorComponent ComponentSelectorConstant >
    struct SpecifyComponent {
        SpecifyComponent( const Vector3& vector ) = delete;
        operator const ScalerType& ( ) = delete;
    };

    SPECIFY_COMPONENT_SPECILIZATION_MACRO( x_, X )
    SPECIFY_COMPONENT_SPECILIZATION_MACRO( y_, Y )
    SPECIFY_COMPONENT_SPECILIZATION_MACRO( z_, Z )
}
