#include "YinYangVolumeObject.h"
#include <kvs/Endian>
#include <fstream>

namespace local
{

kvs::StructuredVolumeObject* YinYangVolumeObject::ToStructuredVolumeObject( const local::YinYangVolumeObject* object )
{
    kvs::StructuredVolumeObject* volume = new kvs::StructuredVolumeObject();
    volume->setGridTypeToCurvilinear();
    volume->setVeclen( object->veclen() );
    volume->setResolution( kvs::Vec3ui( object->dimR(), object->dimTheta(), object->dimPhi() ) );
    volume->setCoords( object->coords() );
    volume->setValues( object->values() );
    volume->updateMinMaxValues();
    volume->updateMinMaxCoords();
    return volume;
}

kvs::UnstructuredVolumeObject* YinYangVolumeObject::ToUnstructuredVolumeObject( const local::YinYangVolumeObject* object )
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

    kvs::UnstructuredVolumeObject* volume = new kvs::UnstructuredVolumeObject();
    volume->setCellTypeToHexahedra();
    volume->setVeclen( object->veclen() );
    volume->setNumberOfNodes( object->numberOfNodes() );
    volume->setNumberOfCells( object->numberOfCells() );
    volume->setCoords( object->coords() );
    volume->setConnections( connections );
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
    const float pi = 3.141593f;

    // r
    const float r_max = 1.0f;
    const float r_min = 0.35f;
    const float r_d = ( r_max - r_min ) / ( m_dim_r - 1 );

    // theta
    const float theta_max = pi - pi / 4.0f;
    const float theta_min = pi / 4.0f;
    const float theta_d = ( theta_max - theta_min ) / ( m_dim_theta - 3 );

    // phi
    const float phi_max = ( 3 * pi ) / 4.0f;
    const float phi_min = -( 3 * pi ) / 4.0f;
    const float phi_d = ( phi_max - phi_min ) / ( m_dim_phi - 5 );

    const size_t nnodes = this->numberOfNodes();
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();
    for ( size_t k = 0; k < m_dim_phi; k++ )
    {
        const float theta = theta_min + theta_d * k;
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );
        for ( size_t j = 0; j < m_dim_theta; j++ )
        {
            const float phi = phi_min + phi_d * j;
            const float sin_phi = std::sin( phi );
            const float cos_phi = std::cos( phi );
            for ( size_t i = 0; i < m_dim_r; i++ )
            {
                const float r = r_min + r_d * i;
                const float x = r * sin_phi * cos_theta;
                const float y = r * cos_phi;
                const float z = r * sin_phi * sin_theta;

                *(pcoords++) = x * ( m_grid_type == Yin ? 1 : -1 );
                *(pcoords++) = y;
                *(pcoords++) = z;
            }
        }
    }

    this->setCoords( coords );
}

bool YinYangVolumeObject::readValues( const std::string& filename )
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
