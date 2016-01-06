#include "YinYangVolumeObject.h"
#include <kvs/Endian>
#include <fstream>


namespace
{

kvs::ValueArray<kvs::Real32> CalculateCoords( const local::YinYangVolumeObject* object )
{
    const size_t dim_r = object->dimR(); // radius
    const size_t dim_theta = object->dimTheta(); // latitude
    const size_t dim_phi = object->dimPhi(); // longitude

    const local::YinYangVolumeObject::Range range_r = object->rangeR();
    const local::YinYangVolumeObject::Range range_theta = object->rangeTheta();
    const local::YinYangVolumeObject::Range range_phi = object->rangePhi();

    const size_t nnodes = object->numberOfNodes();
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();
    for ( size_t k = 0; k < dim_phi; k++ )
    {
        const float theta = range_theta.min + range_theta.d * k;
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );
        for ( size_t j = 0; j < dim_theta; j++ )
        {
            const float phi = range_phi.min + range_phi.d * j;
            const float sin_phi = std::sin( phi );
            const float cos_phi = std::cos( phi );
            for ( size_t i = 0; i < dim_r; i++ )
            {
                const float r = range_r.min + range_r.d * i;
                const float x = r * sin_phi * cos_theta;
                const float y = r * cos_phi;
                const float z = r * sin_phi * sin_theta;

                *(pcoords++) = ( object->gridType() == local::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( object->gridType() == local::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( object->gridType() == local::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    return coords;
}

kvs::ValueArray<kvs::UInt32> CalculateConnections( const local::YinYangVolumeObject* object )
{
    const size_t dim_r = object->dimR(); // radius
    const size_t dim_theta = object->dimTheta(); // latitude
    const size_t dim_phi = object->dimPhi(); // longitude

    const size_t nnodes = object->numberOfNodes();
    kvs::ValueArray<kvs::UInt32> connections( nnodes * 8 );
    kvs::UInt32* pconnections = connections.data();
    for ( size_t k = 0, index = 0; k < dim_phi - 1; k++, index += dim_r )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, index++ )
            {
                *(pconnections++) = index;
                *(pconnections++) = index + 1;
                *(pconnections++) = index + dim_r + 1;
                *(pconnections++) = index + dim_r;
                *(pconnections++) = index + ( dim_r * dim_theta );
                *(pconnections++) = index + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta );
            }
        }
    }

    return connections;
}

}

namespace local
{

kvs::StructuredVolumeObject* YinYangVolumeObject::ToStructuredVolumeObject( const local::YinYangVolumeObject* object )
{
    kvs::StructuredVolumeObject* volume = new kvs::StructuredVolumeObject();
    volume->setGridTypeToCurvilinear();
    volume->setVeclen( object->veclen() );
    volume->setResolution( kvs::Vec3ui( object->dimR(), object->dimTheta(), object->dimPhi() ) );
    volume->setCoords( ::CalculateCoords( object ) );
    volume->setValues( object->values() );
    volume->updateMinMaxValues();
    volume->updateMinMaxCoords();
    return volume;
}

kvs::UnstructuredVolumeObject* YinYangVolumeObject::ToUnstructuredVolumeObject( const local::YinYangVolumeObject* object )
{
    kvs::UnstructuredVolumeObject* volume = new kvs::UnstructuredVolumeObject();
    volume->setCellTypeToHexahedra();
    volume->setVeclen( object->veclen() );
    volume->setNumberOfNodes( object->numberOfNodes() );
    volume->setNumberOfCells( object->numberOfCells() );
    volume->setCoords( ::CalculateCoords( object ) );
    volume->setConnections( ::CalculateConnections( object ) );
    volume->setValues( object->values() );
    volume->updateMinMaxValues();
    volume->updateMinMaxCoords();
    return volume;
}

YinYangVolumeObject::YinYangVolumeObject():
    kvs::VolumeObjectBase(),
    m_grid_type( Yin ),
    m_dim_r( 0 ),
    m_dim_theta( 0 ),
    m_dim_phi( 0 )
{
    BaseClass::setVolumeType( Structured );

    m_range_r.min = 0.0f;
    m_range_r.max = 0.0f;
    m_range_r.d = 0.0f;

    m_range_theta.min = 0.0f;
    m_range_theta.max = 0.0f;
    m_range_theta.d = 0.0f;

    m_range_phi.min = 0.0f;
    m_range_phi.max = 0.0f;
    m_range_phi.d = 0.0f;
}

void YinYangVolumeObject::setDimR( const size_t dim_r, const float range_min, const float range_max )
{
    KVS_ASSERT( dim_r > 1 );

    m_dim_r = dim_r;

    m_range_r.max = range_max;
    m_range_r.min = range_min;
    m_range_r.d = ( m_range_r.max - m_range_r.min ) / ( m_dim_r - 1 );
}

void YinYangVolumeObject::setDimTheta( const size_t dim_theta, const size_t overwrap )
{
    KVS_ASSERT( dim_theta - overwrap > 1 );

    m_dim_theta = dim_theta;

    const float pi = 3.141593f;
    m_range_theta.max = pi - pi / 4.0f;
    m_range_theta.min = pi / 4.0f;
    m_range_theta.d = ( m_range_theta.max - m_range_theta.min ) / ( m_dim_theta - overwrap - 1 );
}

  void YinYangVolumeObject::setDimPhi( const size_t dim_phi, const size_t overwrap )
{
    KVS_ASSERT( dim_phi - overwrap > 1 );

    m_dim_phi = dim_phi;

    const float pi = 3.141593f;
    m_range_phi.max = ( 3 * pi ) / 4.0f;
    m_range_phi.min = -( 3 * pi ) / 4.0f;
    m_range_phi.d = ( m_range_phi.max - m_range_phi.min ) / ( m_dim_phi - overwrap - 1 );
}

size_t YinYangVolumeObject::numberOfNodes() const
{
    return m_dim_r * m_dim_theta * m_dim_phi;
}

size_t YinYangVolumeObject::numberOfCells() const
{
    return ( m_dim_r - 1 ) * ( m_dim_theta - 1 ) * ( m_dim_phi - 1 );
}

bool YinYangVolumeObject::readValues( const std::string& filename )
{
    const size_t nnodes = this->numberOfNodes();
    const size_t veclen = this->veclen();
    kvs::ValueArray<kvs::Real32> values( nnodes * veclen );

    std::ifstream ifs( filename.c_str(), std::ios_base::in | std::ios_base::binary );
    if ( !ifs ) { return false; }

    ifs.seekg( 4,std::ios::beg );
    ifs.read( (char*)values.data(), nnodes * veclen * sizeof( kvs::Real32 ) );
    kvs::Endian::Swap( values.data(), nnodes * veclen );

    this->setValues( kvs::AnyValueArray( values ) );

    return true;
}

} // end of namespace local
