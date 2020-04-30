#pragma once
#include "YinVolumeObject.h"
#include "YangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/Math>


namespace YYZVis
{

inline void UpdateMinMaxCoords(
    YYZVis::YinVolumeObject* yin_volume,
    YYZVis::YangVolumeObject* yng_volume,
    YYZVis::ZhongVolumeObject* zng_volume )
{
    const kvs::Vec3& min_coord0 = yin_volume->minObjectCoord();
    const kvs::Vec3& min_coord1 = yng_volume->minObjectCoord();
    const kvs::Vec3& min_coord2 = zng_volume->minObjectCoord();
    const kvs::Vec3& max_coord0 = yin_volume->maxObjectCoord();
    const kvs::Vec3& max_coord1 = yng_volume->maxObjectCoord();
    const kvs::Vec3& max_coord2 = zng_volume->maxObjectCoord();
    const kvs::Real32 min_x = kvs::Math::Min( min_coord0.x(), min_coord1.x(), min_coord2.x() );
    const kvs::Real32 min_y = kvs::Math::Min( min_coord0.y(), min_coord1.y(), min_coord2.y() );
    const kvs::Real32 min_z = kvs::Math::Min( min_coord0.z(), min_coord1.z(), min_coord2.z() );
    const kvs::Real32 max_x = kvs::Math::Min( max_coord0.x(), max_coord1.x(), max_coord2.x() );
    const kvs::Real32 max_y = kvs::Math::Min( max_coord0.y(), max_coord1.y(), max_coord2.y() );
    const kvs::Real32 max_z = kvs::Math::Min( max_coord0.z(), max_coord1.z(), max_coord2.z() );
    const kvs::Vec3 min_coord( min_x, min_y, min_z );
    const kvs::Vec3 max_coord( max_x, max_y, max_z );
    yin_volume->setMinMaxObjectCoords( min_coord, max_coord );
    yin_volume->setMinMaxExternalCoords( min_coord, max_coord );
    yng_volume->setMinMaxObjectCoords( min_coord, max_coord );
    yng_volume->setMinMaxExternalCoords( min_coord, max_coord );
    zng_volume->setMinMaxObjectCoords( min_coord, max_coord );
    zng_volume->setMinMaxExternalCoords( min_coord, max_coord );
}

} // end of namespace YYZVis
