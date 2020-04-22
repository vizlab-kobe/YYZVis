#include "YinYangGrid.h"


namespace
{

template <typename ValueType>
inline void Bind(
    const YYZVis::YinYangVolumeObject* volume,
    const kvs::Vec3ui& base_index,
    kvs::Real32* grid_values,
    kvs::Vec3* grid_coords )
{
    const size_t dim0 = volume->dimR();
    const size_t dim1 = dim0 * volume->dimTheta();

    const size_t i = base_index[0];
    const size_t j = base_index[1];
    const size_t k = base_index[2];
    const size_t index0 = i + dim0 * j + dim1 * k;

    kvs::UInt32 index[8];
    index[0] = index0;
    index[1] = index[0] + 1;
    index[2] = index[1] + dim0;
    index[3] = index[0] + dim0;
    index[4] = index[0] + dim1;
    index[5] = index[1] + dim1;
    index[6] = index[2] + dim1;
    index[7] = index[3] + dim1;

    const ValueType* const values = volume->values().asValueArray<ValueType>().data();
    grid_values[0] = kvs::Real32( values[ index[0] ] );
    grid_values[1] = kvs::Real32( values[ index[1] ] );
    grid_values[2] = kvs::Real32( values[ index[2] ] );
    grid_values[3] = kvs::Real32( values[ index[3] ] );
    grid_values[4] = kvs::Real32( values[ index[4] ] );
    grid_values[5] = kvs::Real32( values[ index[5] ] );
    grid_values[6] = kvs::Real32( values[ index[6] ] );
    grid_values[7] = kvs::Real32( values[ index[7] ] );

    const kvs::Real32* const coords = volume->coords().data();
    grid_coords[0] = kvs::Vec3( coords + index[0] * 3 );
    grid_coords[1] = kvs::Vec3( coords + index[1] * 3 );
    grid_coords[2] = kvs::Vec3( coords + index[2] * 3 );
    grid_coords[3] = kvs::Vec3( coords + index[3] * 3 );
    grid_coords[4] = kvs::Vec3( coords + index[4] * 3 );
    grid_coords[5] = kvs::Vec3( coords + index[5] * 3 );
    grid_coords[6] = kvs::Vec3( coords + index[6] * 3 );
    grid_coords[7] = kvs::Vec3( coords + index[7] * 3 );
}

inline kvs::Real32 InterpolateValue( const kvs::Real32* values, const kvs::Real32* weights )
{
    return
        weights[0] * values[0] +
        weights[1] * values[1] +
        weights[2] * values[2] +
        weights[3] * values[3] +
        weights[4] * values[4] +
        weights[5] * values[5] +
        weights[6] * values[6] +
        weights[7] * values[7];
}

inline kvs::Vec3 InterpolateCoord( const kvs::Vec3* coords, const kvs::Real32* weights )
{
    return
        weights[0] * coords[0] +
        weights[1] * coords[1] +
        weights[2] * coords[2] +
        weights[3] * coords[3] +
        weights[4] * coords[4] +
        weights[5] * coords[5] +
        weights[6] * coords[6] +
        weights[7] * coords[7];
}

inline kvs::Real32 RandomNumber()
{
    // xorshift RGNs with period at least 2^128 - 1.
    static kvs::Real32 t24 = 1.0/16777216.0; /* 0.5**24 */
    static kvs::UInt32 x=123456789,y=362436069,z=521288629,w=88675123;
    kvs::UInt32 t;
    t=(x^(x<<11));
    x=y;y=z;z=w;
    w=(w^(w>>19))^(t^(t>>8));
    return t24*static_cast<kvs::Real32>(w>>8);
}

}

namespace YYZVis
{

YinYangGrid::YinYangGrid( const YYZVis::YinYangVolumeObject* volume ):
    m_base_index( 0, 0, 0 ),
    m_local_point( 0, 0, 0 ),
    m_reference_volume( volume )
{
    KVS_ASSERT( volume->coords().size() != 0 );

    std::memset( m_coords, 0, sizeof( kvs::Vec3 ) * 8 );
    std::memset( m_values, 0, sizeof( kvs::Real32 ) * 8 );
    std::memset( m_interpolation_functions, 0, sizeof( kvs::Real32 ) * 8 );
    std::memset( m_differential_functions, 0, sizeof( kvs::Real32 ) * 24 );
}

void YinYangGrid::bind( const kvs::Vec3ui& base_index )
{
    KVS_ASSERT( base_index.x() < m_reference_volume->dimR() - 1 );
    KVS_ASSERT( base_index.y() < m_reference_volume->dimTheta() - 1 );
    KVS_ASSERT( base_index.z() < m_reference_volume->dimPhi() - 1 );

    m_base_index = base_index;
    const YYZVis::YinYangVolumeObject* volume = m_reference_volume;
    switch ( volume->values().typeID() )
    {
    case kvs::Type::TypeInt8:   ::Bind<kvs::Int8>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeUInt8:  ::Bind<kvs::UInt8>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeInt16:  ::Bind<kvs::Int16>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeUInt16: ::Bind<kvs::UInt16>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeInt32:  ::Bind<kvs::Int32>( volume,  base_index, m_values, m_coords ); break;
    case kvs::Type::TypeUInt32: ::Bind<kvs::UInt32>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeInt64:  ::Bind<kvs::Int64>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeUInt64: ::Bind<kvs::UInt64>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeReal32: ::Bind<kvs::Real32>( volume, base_index, m_values, m_coords ); break;
    case kvs::Type::TypeReal64: ::Bind<kvs::Real64>( volume, base_index, m_values, m_coords ); break;
    default: break;
    }
}

void YinYangGrid::setLocalPoint( const kvs::Vec3& local ) const
{
    m_local_point = local;
    this->updateInterpolationFunctions( local );
    this->updateDifferentialFunctions( local );
}

void YinYangGrid::updateInterpolationFunctions( const kvs::Vec3& local ) const
{
    const float p = local.x();
    const float q = local.y();
    const float r = local.z();

    const float pq = p * q;
    const float qr = q * r;
    const float rp = r * p;
    const float pqr = pq * r;

    kvs::Real32* N = m_interpolation_functions;
    N[0] = 1.0f - p - q - r + pq + qr + rp - pqr;
    N[1] = p - pq - rp + pqr;
    N[2] = pq - pqr;
    N[3] = q - pq - qr + pqr;
    N[4] = r - rp - qr + pqr;
    N[5] = rp - pqr;
    N[6] = pqr;
    N[7] = qr - pqr;
}

void YinYangGrid::updateDifferentialFunctions( const kvs::Vec3& local ) const
{
    const float p = local.x();
    const float q = local.y();
    const float r = local.z();
    const float pq = p * q;
    const float qr = q * r;
    const float rp = r * p;

    const size_t nnodes = 8;
    kvs::Real32* const dN = m_differential_functions;
    kvs::Real32* const dNdp = dN;
    kvs::Real32* const dNdq = dNdp + nnodes;
    kvs::Real32* const dNdr = dNdq + nnodes;

    dNdp[0] =  - 1.0f + q +r - qr;
    dNdp[1] =  1.0f - q - r + qr;
    dNdp[2] =  q - qr;
    dNdp[3] =  - q + qr;
    dNdp[4] =  - r + qr;
    dNdp[5] =  r - qr;
    dNdp[6] =  qr;
    dNdp[7] =  - qr;

    dNdq[0] =  - 1.0f + p + r - rp;
    dNdq[1] =  - p + rp;
    dNdq[2] =  p - rp;
    dNdq[3] =  1.0f - p - r + rp;
    dNdq[4] =  - r + rp;
    dNdq[5] =  - rp;
    dNdq[6] =  rp;
    dNdq[7] =  r - rp;

    dNdr[0] =  - 1.0f + q + p - pq;
    dNdr[1] =  - p + pq;
    dNdr[2] =  - pq;
    dNdr[3] =  - q + pq;
    dNdr[4] =  1.0f - q - p + pq;
    dNdr[5] =  p - pq;
    dNdr[6] =  pq;
    dNdr[7] =  q - pq;
}

const kvs::Vec3 YinYangGrid::globalPoint() const
{
    const float* N = m_interpolation_functions;
    const kvs::Vec3* V = m_coords;
    return ::InterpolateCoord( V, N );
}

const kvs::Mat3 YinYangGrid::JacobiMatrix() const
{
    const kvs::UInt32 nnodes = 8;
    const float* dNdp = m_differential_functions;
    const float* dNdq = dNdp + nnodes;
    const float* dNdr = dNdq + nnodes;
    const kvs::Vec3* coords = m_coords;

    const kvs::Vec3 dx = ::InterpolateCoord( coords, dNdp );
    const kvs::Vec3 dy = ::InterpolateCoord( coords, dNdq );
    const kvs::Vec3 dz = ::InterpolateCoord( coords, dNdr );
    return kvs::Mat3( dx[0], dy[0], dz[0], dx[1], dy[1], dz[1], dx[2], dy[2], dz[2] );
}

const kvs::Vec3 YinYangGrid::center() const
{
    return (
        m_coords[0] +
        m_coords[1] +
        m_coords[2] +
        m_coords[3] +
        m_coords[4] +
        m_coords[5] +
        m_coords[6] +
        m_coords[7] ) / 8.0f;
}

const kvs::Real32 YinYangGrid::volume() const
{
    const size_t resolution = 3;
    const float sampling_length = 1.0f / (float)resolution;
    const float adjustment = sampling_length * 0.5f;

    kvs::Vec3 sampling_position( -adjustment, -adjustment, -adjustment );

    float sum_metric = 0;
    for ( size_t k = 0 ; k < resolution ; k++ )
    {
        sampling_position[ 2 ] +=  sampling_length;
        for( size_t j = 0 ; j < resolution ; j++ )
        {
            sampling_position[ 1 ] += sampling_length;
            for( size_t i = 0 ; i < resolution ; i++ )
            {
                sampling_position[ 0 ] += sampling_length;

                this->setLocalPoint( sampling_position );
                const kvs::Mat3 J = this->JacobiMatrix();
                const float metric_element = J.determinant();

                sum_metric += kvs::Math::Abs<float>( metric_element );
            }
            sampling_position[ 0 ] = -adjustment;
        }
        sampling_position[ 1 ] = -adjustment;
    }

    const float resolution3 = resolution * resolution * resolution;
    return sum_metric / resolution3;
}

const kvs::Real32 YinYangGrid::scalar() const
{
    const kvs::Real32* N = m_interpolation_functions;
    const kvs::Real32* S = m_values;
    return ::InterpolateValue( S, N );
}

const kvs::Vec3 YinYangGrid::gradientVector() const
{
    // Calculate a gradient vector in the local coordinate.
    const kvs::UInt32 nnodes = 8;
    const float* dNdp = m_differential_functions;
    const float* dNdq = m_differential_functions + nnodes;
    const float* dNdr = m_differential_functions + nnodes + nnodes;
    const kvs::Real32* S = m_values;

    const float dSdp = ::InterpolateValue( S, dNdp );
    const float dSdq = ::InterpolateValue( S, dNdq );
    const float dSdr = ::InterpolateValue( S, dNdr );
    const kvs::Vec3 g( dSdp, dSdq, dSdr );

    // Calculate a gradient vector in the global coordinate.
    const kvs::Mat3 J = this->JacobiMatrix();

//    float determinant = 0.0f;
//    const kvs::Vec3 G = 3.0f * J.inverted( &determinant ) * g;
//    return kvs::Math::IsZero( determinant ) ? kvs::Vec3::Zero() : G;
    return 3.0f * J.inverted() * g;
}

const kvs::Vec3 YinYangGrid::randomSampling() const
{
    const float p = ::RandomNumber();
    const float q = ::RandomNumber();
    const float r = ::RandomNumber();

    const kvs::Vec3 local( p, q, r );
    this->setLocalPoint( local );
    return this->globalPoint();
}

} // end of namespace YYZVis
