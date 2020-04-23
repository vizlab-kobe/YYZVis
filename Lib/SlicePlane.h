#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include <kvs/Vector3>
#include <kvs/Vector4>
#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"


namespace YYZVis
{

class SlicePlane : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( YYZVis::SlicePlane, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

private:
    kvs::Vec4 m_coefficients; ///< coeficients of a slice plane

public:
    SlicePlane();
    SlicePlane(
        const kvs::VolumeObjectBase* volume,
        const kvs::Vec4& coefficients,
        const kvs::TransferFunction& transfer_function );
    SlicePlane(
        const kvs::VolumeObjectBase* volume,
        const kvs::Vec3& point,
        const kvs::Vec3& normal,
        const kvs::TransferFunction& transfer_function );

    void setPlane( const kvs::Vec4& coefficients );
    void setPlane( const kvs::Vec3& point, const kvs::Vec3& normal );
    SuperClass* exec( const kvs::ObjectBase* object );

private:
    void mapping( const YYZVis::ZhongVolumeObject* zvolume );
    void mapping( const YYZVis::YinYangVolumeObject* yvolume );
    void extract_plane( const YYZVis::YinYangVolumeObject* yvolume );
    void extract_plane( const YYZVis::ZhongVolumeObject* zvolume );
    size_t calculate_hexahedra_table_index( const size_t* local_index ) const;
    float substitute_plane_equation( const kvs::Vec3& vertex ) const;
    const kvs::Vec3 interpolate_vertex( const kvs::Vec3& vertex0, const kvs::Vec3& vertex1 ) const;
    double interpolate_value( const YYZVis::YinYangVolumeObject* yvolume, const size_t index0, const size_t index1 ) const;
    double interpolate_value( const YYZVis::ZhongVolumeObject* zvolume, const size_t index0, const size_t index1 ) const;
};

} // end of namespace YYZVis
