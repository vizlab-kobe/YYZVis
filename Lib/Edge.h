#pragma once
#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/LineObject>

namespace YinYangVis
{

namespace Edge
{

    kvs::LineObject* CreateLineObject( const YinYangVis::YinYangVolumeObject* volume );
    kvs::LineObject* CreateLineObject( const YinYangVis::ZhongVolumeObject* volume );

}

} // end of namespace YinYangVis
