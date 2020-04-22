#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>


namespace YYZVis
{

class Isosurface : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( YYZVis::Isosurface, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

private:
    double m_isolevel; ///< isosurface level
    bool m_duplication; ///< duplication flag

public:
    Isosurface();
    Isosurface(
        const kvs::VolumeObjectBase* volume,
        const double isolevel,
        const SuperClass::NormalType normal_type,
        const bool duplication,
        const kvs::TransferFunction& transfer_function );

    void setIsolevel( const double isolevel ) { m_isolevel = isolevel; }
    SuperClass* exec( const kvs::ObjectBase* object );

private:
    void mapping( const YYZVis::ZhongVolumeObject* zvolume );
    void mapping( const YYZVis::YinYangVolumeObject* yvolume );
    void extract_yinyang_surfaces_with_duplication( const YYZVis::YinYangVolumeObject* yvolume );
    void extract_zhong_surfaces_with_duplication( const YYZVis::ZhongVolumeObject* zvolume );
    size_t calculate_table_index( const kvs::AnyValueArray values, const size_t* local_index ) const;

    const kvs::Vec3 interpolate_vertex( const kvs::AnyValueArray values,
					const kvs::Real32* coords,			                                                                                                      const int vertex0,
					const int vertex1 ) const;
    const kvs::RGBColor calculate_color();
};

} // end of namespace YYZVis
