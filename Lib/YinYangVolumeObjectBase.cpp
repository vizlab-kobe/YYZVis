#include "YinYangVolumeObjectBase.h"
#include <kvs/Endian>
#include <fstream>


namespace
{

const std::string GridTypeName[2] = {
    "yin",
    "yang"
};

kvs::ValueArray<kvs::Real32> CalculateCoords( const YYZVis::YinYangVolumeObjectBase* object )
{
    const size_t dim_r = object->dimR(); // radius
    const size_t dim_theta = object->dimTheta(); // latitude
    const size_t dim_phi = object->dimPhi(); // longitude

    const YYZVis::YinYangVolumeObjectBase::Range range_r = object->rangeR();
    const YYZVis::YinYangVolumeObjectBase::Range range_theta = object->rangeTheta();
    const YYZVis::YinYangVolumeObjectBase::Range range_phi = object->rangePhi();

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

                *(pcoords++) = ( object->gridType() == YYZVis::YinYangVolumeObjectBase::Yin ) ? x : -x;
                *(pcoords++) = ( object->gridType() == YYZVis::YinYangVolumeObjectBase::Yin ) ? y : z;
                *(pcoords++) = ( object->gridType() == YYZVis::YinYangVolumeObjectBase::Yin ) ? z : y;
            }
        }
    }

    return coords;
}

kvs::ValueArray<kvs::UInt32> CalculateConnections( const YYZVis::YinYangVolumeObjectBase* object )
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

namespace YYZVis
{

kvs::StructuredVolumeObject* YinYangVolumeObjectBase::ToStructuredVolumeObject(
    const YYZVis::YinYangVolumeObjectBase* object )
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

kvs::UnstructuredVolumeObject* YinYangVolumeObjectBase::ToUnstructuredVolumeObject(
    const YYZVis::YinYangVolumeObjectBase* object )
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

YinYangVolumeObjectBase::YinYangVolumeObjectBase():
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

void YinYangVolumeObjectBase::shallowCopy( const YinYangVolumeObjectBase& object )
{
    BaseClass::shallowCopy( object );
    m_grid_type = object.m_grid_type;
    m_dim_r = object.m_dim_r;
    m_dim_theta = object.m_dim_theta;
    m_dim_phi = object.m_dim_phi;
    m_range_r = object.m_range_r;
    m_range_theta = object.m_range_theta;
    m_range_phi = object.m_range_phi;
}

void YinYangVolumeObjectBase::deepCopy( const YinYangVolumeObjectBase& object )
{
    BaseClass::deepCopy( object );
    m_grid_type = object.m_grid_type;
    m_dim_r = object.m_dim_r;
    m_dim_theta = object.m_dim_theta;
    m_dim_phi = object.m_dim_phi;
    m_range_r = object.m_range_r;
    m_range_theta = object.m_range_theta;
    m_range_phi = object.m_range_phi;
}

void YinYangVolumeObjectBase::print( std::ostream& os, const kvs::Indent& indent ) const
{
    if ( !this->hasMinMaxValues() ) this->updateMinMaxValues();
    os << indent << "Object type : " << "yin-yang volume object" << std::endl;
    BaseClass::print( os, indent );
    os << indent << "Grid type : " << ::GridTypeName[ this->gridType() ] << std::endl;
    os << indent << "Dimension : " << kvs::Vec3ui( this->dimR(), this->dimTheta(), this->dimPhi() ) << std::endl;
    os << indent << "Range (R): " << "[" << this->rangeR().min << ", " << this->rangeR().max << "]" << std::endl;
    os << indent << "Range (Theta): " << "[" << this->rangeTheta().min << ", " << this->rangeTheta().max << "]" << std::endl;
    os << indent << "Range (Phi): " << "[" << this->rangePhi().min << ", " << this->rangePhi().max << "]" << std::endl;
    os << indent << "Number of nodes : " << this->numberOfNodes() << std::endl;
    os << indent << "Number of cells : " << this->numberOfCells() << std::endl;
    os << indent << "Min. value : " << this->minValue() << std::endl;
    os << indent << "Max. value : " << this->maxValue() << std::endl;
}

void YinYangVolumeObjectBase::setDimR( const size_t dim_r, const float range_min, const float range_max )
{
    KVS_ASSERT( dim_r > 1 );

    m_dim_r = dim_r;

    m_range_r.max = range_max;
    m_range_r.min = range_min;
    m_range_r.d = ( m_range_r.max - m_range_r.min ) / ( m_dim_r - 1 );
}

void YinYangVolumeObjectBase::setDimTheta( const size_t dim_theta, const size_t overwrap )
{
    KVS_ASSERT( dim_theta - overwrap > 1 );

    m_dim_theta = dim_theta;

    const float pi = 3.141593f;
    m_range_theta.max = pi - pi / 4.0f;
    m_range_theta.min = pi / 4.0f;
    m_range_theta.d = ( m_range_theta.max - m_range_theta.min ) / ( m_dim_theta - overwrap - 1 );
}

void YinYangVolumeObjectBase::setDimPhi( const size_t dim_phi, const size_t overwrap )
{
    KVS_ASSERT( dim_phi - overwrap > 1 );

    m_dim_phi = dim_phi;

    const float pi = 3.141593f;
    m_range_phi.max = ( 3 * pi ) / 4.0f;
    m_range_phi.min = -( 3 * pi ) / 4.0f;
    m_range_phi.d = ( m_range_phi.max - m_range_phi.min ) / ( m_dim_phi - overwrap - 1 );
}

size_t YinYangVolumeObjectBase::numberOfNodes() const
{
    return m_dim_r * m_dim_theta * m_dim_phi;
}

size_t YinYangVolumeObjectBase::numberOfCells() const
{
    return ( m_dim_r - 1 ) * ( m_dim_theta - 1 ) * ( m_dim_phi - 1 );
}

void YinYangVolumeObjectBase::calculateCoords()
{
    this->setCoords( ::CalculateCoords( this ) );
}

bool YinYangVolumeObjectBase::readValues( const std::string& filename )
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

void YinYangVolumeObjectBase::updateMinMaxCoords()
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

} // end of namespace YYZVis
