#include "Convert.h"
#include "sgrid.h"

namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinYangVolumeObject& yin_volume,
    YYZVis::YinYangVolumeObject& yang_volume,
    YYZVis::ZhongVolumeObject& zhong_volume )
{
    // Vector length
    const int veclen = 1;
    int dim_r, dim_t, dim_p;

    local::Sgrid sgrid( yin_volume, yang_volume, zhong_volume );
 
    dim_r = sgrid.sgrid__size.nr;
    dim_t = sgrid.sgrid__size.nt;
    dim_p = sgrid.sgrid__size.np;

    // Grid resolution (rx, ry, rz)
    const kvs::Vec3u resolution( dim_r, dim_t, dim_p );

    // Values (physical quantity)
    kvs::ValueArray<kvs::Real32> values( sgrid.sgrid__values );
    
    // Coordinate value
    kvs::ValueArray<kvs::Real32> coords( sgrid.sgrid__coords );

    kvs::StructuredVolumeObject object;
    object.setGridTypeToCurvilinear();
    object.setVeclen( veclen );
    object.setResolution( resolution );
    object.setValues( values );
    object.setCoords( coords );

    return object;
}

    
} // end of namespace local
