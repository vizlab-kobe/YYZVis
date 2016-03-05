#pragma once
#include <kvs/MapperBase>
#include <kvs/Camera>
#include <kvs/PointObject>
#include <kvs/VolumeObjectBase>
#include <kvs/Module>
#include "ZhongVolumeObject.h"


namespace YinYangVis
{

class ZhongGridSampling : public kvs::MapperBase, public kvs::PointObject
{
    kvsModule( YinYangVis::ZhongGridSampling, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PointObject );

private:
    const kvs::Camera* m_camera; ///< camera (reference)
    size_t m_subpixel_level; ///< subpixel level
    float m_sampling_step; ///< sampling step in the object coordinate
    float m_object_depth; ///< object depth

public:
    ZhongGridSampling(
        const kvs::VolumeObjectBase* volume,
        const size_t subpixel_level,
        const float sampling_step,
        const kvs::TransferFunction& transfer_function,
        const float object_depth = 0.0f );

    ZhongGridSampling(
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
    void mapping( const YinYangVis::ZhongVolumeObject* volume );
};

} // end of namespace YinYangVis
