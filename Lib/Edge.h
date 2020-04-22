#pragma once
#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/LineObject>


namespace YYZVis
{

namespace Edge
{

kvs::LineObject* CreateLineMeshObject( const YYZVis::YinYangVolumeObject* volume, const size_t dim_edge = 10 );
kvs::LineObject* CreateLineEdgeObject( const YYZVis::YinYangVolumeObject* volume);
kvs::LineObject* CreateLineEdgeObject( const YYZVis::ZhongVolumeObject* volume );

}

} // end of namespace YYZVis
