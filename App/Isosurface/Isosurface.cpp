#include "Isosurface.h"


namespace YinYangVis
{

Isosurface::Isosurface():
    kvs::MapperBase(),
    kvs::PolygonObject(),
    m_isolevel( 0 ),
    m_duplication( true )
{
}

Isosurface::Isosurface(
    const kvs::VolumeObjectBase* volume,
    const double isolevel,
    const SuperClass::NormalType normal_type,
    const bool duplication,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject(),
    m_duplication( duplication )
{
    SuperClass::setNormalType( normal_type );
    this->setIsolevel( isolevel );
    this->exec( volume );
}

Isosurface::SuperClass* Isosurface::exec( const kvs::ObjectBase* object )
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

    // In the case of VertexNormal-type, the duplicated vertices are forcibly deleted.
    if ( SuperClass::normalType() == kvs::PolygonObject::VertexNormal )
    {
        m_duplication = false;
    }

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    const kvs::Real64 min_value = BaseClass::volume()->minValue();
    const kvs::Real64 max_value = BaseClass::volume()->maxValue();
    if ( kvs::Math::Equal( min_value, max_value ) ) { return this; }

    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
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

void Isosurface::mapping( const YinYangVis::ZhongVolumeObject* zvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
    // if ( m_duplication )
    // {
    //     ...
    // }
    // else
    // {
    //     ...
    // }

    // const kvs::ColorMap cmap = BaseClass::colorMap();
    // const kvs::RGBColor color = COLOR_OF_ISOSURFACES;

    // SuperClass::setCoords( coords );
    // SuperClass::setColor( color );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void Isosurface::mapping( const YinYangVis::YinYangVolumeObject* yvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
    // if ( m_duplication )
    // {
    //     ...
    // }
    // else
    // {
    //     ...
    // }

    // const kvs::ColorMap cmap = BaseClass::colorMap();
    // const kvs::RGBColor color = COLOR_OF_ISOSURFACES;

    // SuperClass::setCoords( coords );
    // SuperClass::setColor( color );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

} // end of namespace YinYangVis
