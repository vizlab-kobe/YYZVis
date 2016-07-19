#include "YinYangVolumeObject.h"
#include <kvs/Endian>
#include <fstream>


namespace
{

kvs::ValueArray<kvs::Real32> CalculateCoords( const YinYangVis::YinYangVolumeObject* object )
{
    const size_t dim_r = object->dimR(); // radius
    const size_t dim_theta = object->dimTheta(); // latitude
    const size_t dim_phi = object->dimPhi(); // longitude

    const YinYangVis::YinYangVolumeObject::Range range_r = object->rangeR();
    const YinYangVis::YinYangVolumeObject::Range range_theta = object->rangeTheta();
    const YinYangVis::YinYangVolumeObject::Range range_phi = object->rangePhi();

    const size_t nnodes = object->numberOfNodes();
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();
    for ( int k = 0; k < (int)dim_phi; k++ )
    {
      const float phi = range_phi.min + range_phi.d * ( k - 2 );
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );
        for ( int j = 0; j < (int)dim_theta; j++ )
        {
	    const float theta = range_theta.min + range_theta.d * ( j - 1 );
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            for ( int i = 0; i < (int)dim_r; i++ )
            {
                const float r = range_r.min + range_r.d * i;
                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( object->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( object->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( object->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    return coords;
}

kvs::ValueArray<kvs::UInt32> CalculateConnections( const YinYangVis::YinYangVolumeObject* object )
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
                *(pconnections++) = index + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + ( dim_r * dim_theta );
                *(pconnections++) = index + dim_r;
                *(pconnections++) = index + dim_r + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta );
            }
        }
    }

    return connections;
}

}

namespace YinYangVis
{

kvs::StructuredVolumeObject* YinYangVolumeObject::ToStructuredVolumeObject( const YinYangVis::YinYangVolumeObject* object )
{
    kvs::StructuredVolumeObject* volume = new kvs::StructuredVolumeObject();
    volume->setGridTypeToCurvilinear();
    volume->setVeclen( object->veclen() );
    volume->setResolution( kvs::Vec3ui( object->dimR(), object->dimTheta(), object->dimPhi() ) );
    volume->setCoords( object->coords().size() == 0 ? ::CalculateCoords( object ) : object->coords() );
    volume->setValues( object->values() );
    volume->setMinMaxValues( object->minValue(), object->maxValue() );
    volume->setMinMaxObjectCoords( object->minObjectCoord(), object->maxObjectCoord() );
    volume->setMinMaxExternalCoords( object->minExternalCoord(), object->maxExternalCoord() );
//    volume->updateMinMaxValues();
//    volume->updateMinMaxCoords();
    return volume;
}

kvs::UnstructuredVolumeObject* YinYangVolumeObject::ToUnstructuredVolumeObject( const YinYangVis::YinYangVolumeObject* object )
{
    kvs::UnstructuredVolumeObject* volume = new kvs::UnstructuredVolumeObject();
    volume->setCellTypeToHexahedra();
    volume->setVeclen( object->veclen() );
    volume->setNumberOfNodes( object->numberOfNodes() );
    volume->setNumberOfCells( object->numberOfCells() );
    volume->setCoords( object->coords().size() == 0 ? ::CalculateCoords( object ) : object->coords() );
    volume->setConnections( ::CalculateConnections( object ) );
    volume->setValues( object->values() );
    volume->setMinMaxValues( object->minValue(), object->maxValue() );
    volume->setMinMaxObjectCoords( object->minObjectCoord(), object->maxObjectCoord() );
    volume->setMinMaxExternalCoords( object->minExternalCoord(), object->maxExternalCoord() );
//    volume->updateMinMaxValues();
//    volume->updateMinMaxCoords();
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

void YinYangVolumeObject::calculateCoords()
{
    this->setCoords( ::CalculateCoords( this ) );
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

void YinYangVolumeObject::updateMinMaxCoords()
{
    kvs::Vec3 min_coord( 0.0f, 0.0f, 0.0f );
    kvs::Vec3 max_coord( 0.0f, 0.0f, 0.0f );

    const float* coord = this->coords().data();
    const float* const end = coord + this->coords().size();

    float x = *( coord++ );
    float y = *( coord++ );
    float z = *( coord++ );

    min_coord.set( x, y, z );
    max_coord.set( x, y, z );

    while ( coord < end )
    {
        x = *( coord++ );
        y = *( coord++ );
        z = *( coord++ );

        min_coord.x() = kvs::Math::Min( min_coord.x(), x );
        min_coord.y() = kvs::Math::Min( min_coord.y(), y );
        min_coord.z() = kvs::Math::Min( min_coord.z(), z );

        max_coord.x() = kvs::Math::Max( max_coord.x(), x );
        max_coord.y() = kvs::Math::Max( max_coord.y(), y );
        max_coord.z() = kvs::Math::Max( max_coord.z(), z );
    }

    this->setMinMaxObjectCoords( min_coord, max_coord );

    if ( !( this->hasMinMaxExternalCoords() ) )
    {
        this->setMinMaxExternalCoords(
            this->minObjectCoord(),
            this->maxObjectCoord() );
    }
}

} // end of namespace YinYangVis
