#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>


namespace local
{

/*===========================================================================*/
/**
 *  @brief  External face extraction class.
 */
/*===========================================================================*/
class ExternalFaces : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( local::ExternalFaces, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

public:
    ExternalFaces() {}
    ExternalFaces( const kvs::VolumeObjectBase* volume );
    ExternalFaces( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& tfunc );
    SuperClass* exec( const kvs::ObjectBase* object );

private:
    void mapping( const YYZVis::ZhongVolumeObject* zvolume );
    void mapping( const YYZVis::YinYangVolumeObject* yvolume );
    void calculate_coords( const YYZVis::ZhongVolumeObject* zvolume );
    void calculate_coords( const YYZVis::YinYangVolumeObject* yvolume );
    void calculate_colors( const YYZVis::YinYangVolumeObject* yvolume );
    void calculate_colors( const YYZVis::ZhongVolumeObject* zvolume );
};

} // end of namespace local
