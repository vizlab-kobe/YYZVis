#include "Convert.h"
#include "sgrid.h"

namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinYangVolumeObject& yin_volume,
    YYZVis::YinYangVolumeObject& yang_volume,
    YYZVis::ZhongVolumeObject& zhong_volume )
{
    local::Sgrid sgrid( yin_volume, yang_volume, zhong_volume );

    const size_t dim_r = sgrid.sgrid__size.nr;
    const size_t dim_t = sgrid.sgrid__size.nt;
    const size_t dim_p = sgrid.sgrid__size.np;
    const kvs::Vec3u resolution( dim_r, dim_t, dim_p );
    const kvs::ValueArray<kvs::Real32> values( sgrid.sgrid__values );
    const kvs::ValueArray<kvs::Real32> coords( sgrid.sgrid__coords );
    const int veclen = 1;

    kvs::StructuredVolumeObject object;
    object.setGridTypeToCurvilinear();
    object.setVeclen( veclen );
    object.setResolution( resolution );
    object.setValues( values );
    object.setCoords( coords );
    return object;
}

} // end of namespace local
