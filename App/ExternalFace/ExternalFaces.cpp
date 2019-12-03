#include "ExternalFaces.h"


namespace YinYangVis
{

ExternalFaces::ExternalFaces():
    kvs::MapperBase(),
    kvs::PolygonObject()
{
}

ExternalFaces::ExternalFaces( const kvs::VolumeObjectBase* volume ):
    kvs::MapperBase(),
    kvs::PolygonObject()
{
    this->exec( volume );
}

ExternalFaces::ExternalFaces( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->exec( volume );
}

ExternalFaces::SuperClass* ExternalFaces::exec( const kvs::ObjectBase* object )
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

void ExternalFaces::mapping( const YinYangVis::ZhongVolumeObject* zvolume )
{
    // const size_t nfaces = NUMBER_OF_FACES;
    // const size_t nverts = NUMBER_OF_VERTICES;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // kvs::ValueArray<kvs::Real32> coords( nverts * 3 );
    // kvs::ValueArray<kvs::UInt8> colors( nverts * 3 );
    // kvs::ValueArray<kvs::Real32> normals( nfaces * 3 );

    // Calculate coords, colors, and normals.
    // ...

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

void ExternalFaces::mapping( const YinYangVis::YinYangVolumeObject* yvolume )
{
    // const size_t nfaces = NUMBER_OF_FACES;
    // const size_t nverts = NUMBER_OF_VERTICES;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // kvs::ValueArray<kvs::Real32> coords( nverts * 3 );
    // kvs::ValueArray<kvs::UInt8> colors( nverts * 3 );
    // kvs::ValueArray<kvs::Real32> normals( nfaces * 3 );

    // Calculate coords, colors, and normals.
    // ...

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

} // end of namespace YinYangVis
