#include "SlicePlane.h"


namespace YinYangVis
{

SlicePlane::SlicePlane():
    kvs::MapperBase(),
    kvs::PolygonObject()
{
}

SlicePlane::SlicePlane(
    const kvs::VolumeObjectBase* volume,
    const kvs::Vec4& coefficients,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->setPlane( coefficients );
    this->exec( volume );
}

SlicePlane::SlicePlane(
    const kvs::VolumeObjectBase* volume,
    const kvs::Vec3& point,
    const kvs::Vec3& normal,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->setPlane( point, normal );
    this->exec( volume );
}

void SlicePlane::setPlane( const kvs::Vector4f& coefficients )
{
    m_coefficients = coefficients;
}

void SlicePlane::setPlane( const kvs::Vector3f& point, const kvs::Vector3f& normal )
{
    m_coefficients = kvs::Vec4( normal, -point.dot( normal ) );
}

SlicePlane::SuperClass* SlicePlane::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is NULL.");
        return NULL;
    }

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is not volume dat.");
        return NULL;
    }

    // Check whether the volume can be processed or not.
    if ( volume->veclen() != 1 )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("The input volume is not a sclar field data.");
        return NULL;
    }

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );

    if ( YinYangVis::ZhongVolumeObject::DownCast( volume ) )
    {
        const YinYangVis::ZhongVolumeObject* zvolume = YinYangVis::ZhongVolumeObject::DownCast( volume );
        this->mapping( zvolume );
    }
    else if ( YinYangVis::YinYangVolumeObject::DownCast( volume ) )
    {
        const YinYangVis::YinYangVolumeObject* yvolume = YinYangVis::YinYangVolumeObject::DownCast( volume );
        this->mapping( yvolume );
    }

    return this;
}

void SlicePlane::mapping( const YinYangVis::ZhongVolumeObject* zvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;
    // std::vector<kvs::UInt8> colors;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // Calculate coords, normals and colors.
    // ...

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void SlicePlane::mapping( const YinYangVis::YinYangVolumeObject* yvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;
    // std::vector<kvs::UInt8> colors;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // Calculate coords, normals and colors.
    // ...

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

} // end of namespace YinYangVis
