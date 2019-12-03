#pragma once
#include "Input.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
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
    typedef YinYangVis::YinYangVolumeObject YinVolume;
    typedef YinYangVis::YinYangVolumeObject YangVolume;
    typedef YinYangVis::ZhongVolumeObject ZhongVolume;

private:
    local::Input& m_input; ///< input information from commandline arguments
    YinVolume m_yin_volume; ///< yin volume data
    YangVolume m_yang_volume; ///< yang volume data
    ZhongVolume m_zhong_volume; ///< zhong volume data
    float m_isovalue; ///< value for isosurface extraction

public:
    Model( local::Input& input );

    const local::Input& input() const { return m_input; }
    const YinVolume& constYinVolume() const { return m_yin_volume; }
    const YangVolume& constYangVolume() const { return m_yang_volume; }
    const ZhongVolume& constZhongVolume() const { return m_zhong_volume; }

    float isovalue() const { return m_isovalue; }
    void setIsovalue( const float value ) { m_isovalue = value; }

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

//    kvs::PolygonObject* newYinIsosurfaces() const;
//    kvs::PolygonObject* newYangIsosurfaces() const;
//    kvs::PolygonObject* newZhongIsosurfaces() const;
//    kvs::PolygonObject* newIsosurfaces( const kvs::UnstructuredVolumeObject* volume ) const;

private:
    void import_yin_volume();
    void import_yang_volume();
    void import_zhong_volume();
    void update_min_max_values();
    void update_min_max_coords();
};

} // end of namespace local
