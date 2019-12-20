#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include <kvs/Vector3>
#include <kvs/Vector4>
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>


namespace YinYangVis
{

class SlicePlane : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( YinYangVis::SlicePlane, Mapper );
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
    void mapping( const YinYangVis::ZhongVolumeObject* zvolume );
    void mapping( const YinYangVis::YinYangVolumeObject* yvolume );
    void extract_yinyang_plane( const YinYangVis::YinYangVolumeObject* yvolume );
    void extract_zhong_plane( const YinYangVis::ZhongVolumeObject* zvolume );
    size_t calculate_hexahedra_table_index( const size_t* local_index ) const;
    float substitute_plane_equation( const kvs::Vector3f& vertex ) const;
    const kvs::Vector3f interpolate_vertex(
                                           const kvs::Vector3f& vertex0,
                                           const kvs::Vector3f& vertex1 ) const;
    double interpolate_yinyang_value(
			     const YinYangVis::YinYangVolumeObject* yvolume,
                             const size_t                         index0,
                             const size_t                         index1 ) const;
    double interpolate_zhong_value(
			     const YinYangVis::ZhongVolumeObject* zvolume,
                             const size_t                         index0,
                             const size_t                         index1 ) const;
};

} // end of namespace YinYangVis
