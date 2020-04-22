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
    void calculate_yinyang_coords( const YYZVis::YinYangVolumeObject* yvolume );
    void calculate_normal(
        const float x0, const float y0, const float z0,
        const float x1, const float y1, const float z1,
        const float x2, const float y2, const float z2, kvs::Real32* normal, size_t index );
    void calculate_yinyang_colors( const YYZVis::YinYangVolumeObject* yvolume );
    void calculate_colors( const YYZVis::ZhongVolumeObject* zvolume );
    void GetColorIndices(
        kvs::AnyValueArray value,
        const kvs::Real64 min_value,
        const kvs::Real64 max_value,
        const size_t colormap_resolution,
        const kvs::UInt32 node_index[4],
        kvs::UInt32 (*color_index)[4]);
};

} // end of namespace local
