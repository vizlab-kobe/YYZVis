#pragma once
#include <kvs/MapperBase>
#include <kvs/Camera>
#include <kvs/PointObject>
#include <kvs/VolumeObjectBase>
#include <kvs/Module>
#include "YinYangVolumeObject.h"


namespace local
{

class YinYangGridSampling : public kvs::MapperBase, public kvs::PointObject
{
    kvsModule( local::YinYangGridSampling, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PointObject );

private:
    const kvs::Camera* m_camera; ///< camera (reference)
    size_t m_repeat_level; ///< repeat level
    float m_sampling_step; ///< sampling step in the object coordinate
    float m_object_depth; ///< object depth

public:
    YinYangGridSampling(
        const kvs::VolumeObjectBase* volume,
        const size_t repeat_level,
        const float sampling_step,
        const kvs::TransferFunction& transfer_function,
        const float object_depth = 0.0f );

    YinYangGridSampling(
        const kvs::Camera* camera,
        const kvs::VolumeObjectBase* volume,
        const size_t repeat_level,
        const float sampling_step,
        const kvs::TransferFunction& transfer_function,
        const float object_depth = 0.0f );

    SuperClass* exec( const kvs::ObjectBase* object );

    size_t repeatLevel() const { return m_repeat_level; }
    float samplingStep() const { return m_sampling_step; }
    float objectDepth() const { return m_object_depth; }

    void attachCamera( const kvs::Camera* camera ) { m_camera = camera; }
    void setRepeatLevel( const size_t repeat_level ) { m_repeat_level = repeat_level; }
    void setSamplingStep( const float step ) { m_sampling_step = step; }
    void setObjectDepth( const float depth ) { m_object_depth = depth; }

private:
    void mapping( const local::YinYangVolumeObject* volume );
};

} // end of namespace local
