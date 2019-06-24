#pragma once
#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/LineObject>


namespace YinYangVis
{

namespace Edge
{

kvs::LineObject* CreateLineMeshObject( const YinYangVis::YinYangVolumeObject* volume, const size_t dim_edge = 10 );
kvs::LineObject* CreateLineEdgeObject( const YinYangVis::YinYangVolumeObject* volume);
kvs::LineObject* CreateLineEdgeObject( const YinYangVis::ZhongVolumeObject* volume );

}

} // end of namespace YinYangVis
