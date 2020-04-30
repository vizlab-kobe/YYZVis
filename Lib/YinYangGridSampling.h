#pragma once
#include <kvs/MapperBase>
#include <kvs/Camera>
#include <kvs/PointObject>
#include <kvs/VolumeObjectBase>
#include <kvs/Module>
#include "YinYangVolumeObjectBase.h"


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Partilce generation class for Yin-Yang grid.
 */
/*===========================================================================*/
class YinYangGridSampling : public kvs::MapperBase, public kvs::PointObject
{
    kvsModule( YYZVis::YinYangGridSampling, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PointObject );

private:
    const kvs::Camera* m_camera; ///< camera (reference)
    size_t m_subpixel_level; ///< subpixel level
    float m_sampling_step; ///< sampling step in the object coordinate
    float m_object_depth; ///< object depth

public:
    YinYangGridSampling(
        const kvs::VolumeObjectBase* volume,
        const size_t subpixel_level,
        const float sampling_step,
        const kvs::TransferFunction& transfer_function,
        const float object_depth = 0.0f );

    YinYangGridSampling(
        const kvs::Camera* camera,
        const kvs::VolumeObjectBase* volume,
        const size_t subpixel_level,
        const float sampling_step,
        const kvs::TransferFunction& transfer_function,
        const float object_depth = 0.0f );

    SuperClass* exec( const kvs::ObjectBase* object );

    size_t subpixelLevel() const { return m_subpixel_level; }
    float samplingStep() const { return m_sampling_step; }
    float objectDepth() const { return m_object_depth; }

    void attachCamera( const kvs::Camera* camera ) { m_camera = camera; }
    void setSubpixelLevel( const size_t subpixel_level ) { m_subpixel_level = subpixel_level; }
    void setSamplingStep( const float step ) { m_sampling_step = step; }
    void setObjectDepth( const float depth ) { m_object_depth = depth; }

private:
    void mapping_metro_yin( const YYZVis::YinYangVolumeObjectBase* volume );
    void mapping_metro_yang( const YYZVis::YinYangVolumeObjectBase* volume );
    void mapping_uniform( const YYZVis::YinYangVolumeObjectBase* volume );
};

} // end of namespace YYZVis
