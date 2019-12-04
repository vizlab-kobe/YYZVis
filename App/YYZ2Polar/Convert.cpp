#include "Convert.h"
#include "sgrid.h"


namespace local
{

kvs::StructuredVolumeObject Convert(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume )
{
    local::Sgrid sgrid( yin_volume, yang_volume, zhong_volume );

    // Vector length
    const int veclen = 1;

    // Grid resolution (rx, ry, rz)
    int dim_r = sgrid.sgrid__size.nr;
    int dim_t = sgrid.sgrid__size.nt;
    int dim_p = sgrid.sgrid__size.np;
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
    object.updateMinMaxCoords();

    return object;
}

} // end of namespace local
