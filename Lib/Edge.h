#pragma once
#include "YinYangVolumeObjectBase.h"
#include "ZhongVolumeObject.h"
#include <kvs/LineObject>


namespace YYZVis
{

namespace Edge
{

kvs::LineObject* CreateLineMeshObject( const YYZVis::YinYangVolumeObjectBase* volume, const size_t dim_edge = 10 );
kvs::LineObject* CreateLineEdgeObject( const YYZVis::YinYangVolumeObjectBase* volume);
kvs::LineObject* CreateLineEdgeObject( const YYZVis::ZhongVolumeObject* volume );

}

} // end of namespace YYZVis
