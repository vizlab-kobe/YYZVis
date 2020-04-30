#pragma once
#include "YinVolumeObject.h"
#include "YangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/Math>


namespace YYZVis
{

inline void UpdateMinMaxValues(
    YYZVis::YinVolumeObject* yin_volume,
    YYZVis::YangVolumeObject* yng_volume,
    YYZVis::ZhongVolumeObject* zng_volume )
{
    const kvs::Real32 min_value0 = yin_volume->minValue();
    const kvs::Real32 min_value1 = yng_volume->minValue();
    const kvs::Real32 min_value2 = zng_volume->minValue();
    const kvs::Real32 max_value0 = yin_volume->maxValue();
    const kvs::Real32 max_value1 = yng_volume->maxValue();
    const kvs::Real32 max_value2 = zng_volume->maxValue();
    const kvs::Real32 min_value = kvs::Math::Min( min_value0, min_value1, min_value2 );
    const kvs::Real32 max_value = kvs::Math::Max( max_value0, max_value1, max_value2 );
    yin_volume->setMinMaxValues( min_value, max_value );
    yng_volume->setMinMaxValues( min_value, max_value );
    zng_volume->setMinMaxValues( min_value, max_value );
}

} // end of namespace YYZVis
