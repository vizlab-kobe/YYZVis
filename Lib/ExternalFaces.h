#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include "YinYangVolumeObjectBase.h"
#include "ZhongVolumeObject.h"


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  External face extraction class.
 */
/*===========================================================================*/
class ExternalFaces : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( YYZVis::ExternalFaces, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

public:
    ExternalFaces() {}
    ExternalFaces( const kvs::VolumeObjectBase* volume );
    ExternalFaces( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& tfunc );
    SuperClass* exec( const kvs::ObjectBase* object );

private:
    void mapping( const YYZVis::ZhongVolumeObject* zvolume );
    void mapping( const YYZVis::YinYangVolumeObjectBase* yvolume );
    void calculate_coords( const YYZVis::YinYangVolumeObjectBase* yvolume );
    void calculate_colors( const YYZVis::YinYangVolumeObjectBase* yvolume );
    void calculate_coords( const YYZVis::ZhongVolumeObject* zvolume );
    void calculate_colors( const YYZVis::ZhongVolumeObject* zvolume );
};

} // end of namespace YYZVis
