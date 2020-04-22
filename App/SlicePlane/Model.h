#pragma once
#include "Input.h"
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <kvs/LineObject>
#include <kvs/PolygonObject>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/TransferFunction>


namespace local
{

/*===========================================================================*/
/**
 *  @brief  Model class.
 */
/*===========================================================================*/
class Model
{
public:
    typedef YYZVis::YinYangVolumeObject YinVolume;
    typedef YYZVis::YinYangVolumeObject YangVolume;
    typedef YYZVis::ZhongVolumeObject ZhongVolume;

private:
    const local::Input& m_input; ///< input information from commandline arguments
    YinVolume m_yin_volume; ///< yin volume data
    YangVolume m_yang_volume; ///< yang volume data
    ZhongVolume m_zhong_volume; ///< zhong volume data
    kvs::Vec3 m_plane_point; ///< point on the slice plane
    kvs::Vec3 m_plane_normal; ///< normal vector of the slice plane

public:
    Model( const local::Input& input );

    const local::Input& input() const { return m_input; }
    const YinVolume& constYinVolume() const { return m_yin_volume; }
    const YangVolume& constYangVolume() const { return m_yang_volume; }
    const ZhongVolume& constZhongVolume() const { return m_zhong_volume; }

    const kvs::Vec3& planePoint() const { return m_plane_point; }
    const kvs::Vec3& planeNormal() const { return m_plane_normal; }
    void setPlane( const kvs::Vec3& p, const kvs::Vec3& n ) { m_plane_point = p; m_plane_normal = n; }

    YinVolume* newYinVolume() const { return new YinVolume( m_yin_volume ); }
    YangVolume* newYangVolume() const { return new YangVolume( m_yang_volume ); }
    ZhongVolume* newZhongVolume() const { return new ZhongVolume( m_zhong_volume ); }

    kvs::LineObject* newYinMeshes( const size_t dim_edge = 10 ) const;
    kvs::LineObject* newYangMeshes( const size_t dim_edge = 10 ) const;

    kvs::LineObject* newYinEdges() const;
    kvs::LineObject* newYangEdges() const;
    kvs::LineObject* newZhongEdges() const;

    kvs::PolygonObject* newYinFaces() const;
    kvs::PolygonObject* newYangFaces() const;
    kvs::PolygonObject* newZhongFaces() const;
    kvs::PolygonObject* newFaces( const kvs::UnstructuredVolumeObject* volume ) const;

    kvs::PolygonObject* newYinSlice() const;
    kvs::PolygonObject* newYangSlice() const;
    kvs::PolygonObject* newZhongSlice() const;
    kvs::PolygonObject* newSlice( const kvs::UnstructuredVolumeObject* volume ) const;

private:
    void import_yin_volume();
    void import_yang_volume();
    void import_zhong_volume();
    void update_min_max_values();
    void update_min_max_coords();
};

} // end of namespace local
