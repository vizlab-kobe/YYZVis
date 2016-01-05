#include "ZhongVolumeObject.h"
#include <kvs/Endian>
#include <fstream>


namespace
{

kvs::ValueArray<kvs::Real32> CalculateCoords( const local::ZhongVolumeObject* object )
{
    const size_t dim = object->dim();
    const float r_min = object->rangeR().min;
    const float r_d = object->rangeR().d;

    //ix(= iy, iz),dix(= diy, diz)
    const float r_i = r_min + r_d * 2;
    const float ix_max = r_i;
    const float ix_min = -r_i;
    const float dix = ( ix_max - ix_min ) / ( dim - 1 );

    const size_t nnodes = object->numberOfNodes();
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();
    for ( size_t k = 0; k < dim; k++ )
    {
        const float z = ix_min + dix * k;
        for ( size_t j = 0; j < dim; j++ )
        {
            const float y = ix_min + dix * j;
            for ( size_t i = 0; i < dim; i++ )
            {
                const float x = ix_min + dix * i;
                *(pcoords++) = x;
                *(pcoords++) = y;
                *(pcoords++) = z;
            }
        }
    }

    return coords;
}

kvs::ValueArray<kvs::UInt32> CalculateConnections( const local::ZhongVolumeObject* object )
{
    const size_t dim = object->dim();

    const size_t nnodes = object->numberOfNodes();
    kvs::ValueArray<kvs::UInt32> connections( nnodes * 8 );
    kvs::UInt32* pconnections = connections.data();
    for ( size_t k = 0, index = 0; k < dim - 1; k++, index += dim )
    {
        for ( size_t j = 0; j < dim - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim - 1; i++, index++ )
            {
                *(pconnections++) = index;
                *(pconnections++) = index + 1;
                *(pconnections++) = index + ( dim * dim ) + 1;
                *(pconnections++) = index + ( dim * dim );
                *(pconnections++) = index + dim;
                *(pconnections++) = index + dim + 1;
                *(pconnections++) = index + dim + ( dim * dim ) + 1;
                *(pconnections++) = index + dim + ( dim * dim );
            }
        }
    }

    return connections;
}

}

namespace local
{

kvs::StructuredVolumeObject* ZhongVolumeObject::ToStructuredVolumeObject( const local::ZhongVolumeObject* object )
{
    kvs::StructuredVolumeObject* volume = new kvs::StructuredVolumeObject();
    volume->setGridTypeToUniform();
    volume->setVeclen( object->veclen() );
    volume->setResolution( kvs::Vec3ui( object->dim(), object->dim(), object->dim() ) );
    volume->setValues( object->values() );
    volume->updateMinMaxValues();
    volume->updateMinMaxCoords();
    return volume;
}

kvs::UnstructuredVolumeObject* ZhongVolumeObject::ToUnstructuredVolumeObject( const local::ZhongVolumeObject* object )
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

ZhongVolumeObject::ZhongVolumeObject():
    kvs::VolumeObjectBase(),
    m_dim( 0 ),
    m_dim_r( 0 )
{
    BaseClass::setVolumeType( Structured );

    m_range_r.min = 0.0f;
    m_range_r.max = 0.0f;
    m_range_r.d = 0.0f;
}

  void ZhongVolumeObject::setDimR( const size_t dim_r, const float range_min, const float range_max )
{
    KVS_ASSERT( dim_r > 1 );

    m_dim_r = dim_r;

    m_range_r.max = range_max;
    m_range_r.min = range_min;
    m_range_r.d = ( m_range_r.max - m_range_r.min ) / ( m_dim_r - 1 );
}

size_t ZhongVolumeObject::numberOfNodes() const
{
    return m_dim * m_dim * m_dim;
}

size_t ZhongVolumeObject::numberOfCells() const
{
    return ( m_dim - 1 ) * ( m_dim - 1 ) * ( m_dim - 1 );
}

bool ZhongVolumeObject::readValues( const std::string& filename )
{
    const size_t nnodes = this->numberOfNodes();
    const size_t veclen = this->veclen();
    kvs::ValueArray<kvs::Real32> values( nnodes * veclen );

    std::ifstream ifs( filename.c_str(), std::ios_base::in | std::ios_base::binary );
    if ( !ifs ) { return false; }

    ifs.seekg( 4,std::ios::beg );
    ifs.read( (char*)values.data(), nnodes * veclen * sizeof(kvs::Real32) );
    kvs::Endian::Swap( values.data(), nnodes * veclen );

    this->setValues( kvs::AnyValueArray( values ) );

    return true;
}

} // end of namespace local
