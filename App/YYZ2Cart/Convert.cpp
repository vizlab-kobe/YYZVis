#include "Convert.h"
#include "cgrid.h"

namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinYangVolumeObject& yin_volume,
    YYZVis::YinYangVolumeObject& yang_volume,
    YYZVis::ZhongVolumeObject& zhong_volume )
{
    local::Cgrid cgrid( yin_volume, yang_volume, zhong_volume );

    const size_t dim_x = cgrid.cgrid__size.nx;
    const size_t dim_y = cgrid.cgrid__size.ny;
    const size_t dim_z = cgrid.cgrid__size.nz;
    const kvs::Vec3u resolution( dim_x, dim_y, dim_z );
    const kvs::ValueArray<kvs::Real32> values( cgrid.cgrid__values );
    const kvs::ValueArray<kvs::Real32> coords( cgrid.cgrid__coords );
    const int veclen = 1;

    kvs::StructuredVolumeObject object;
    object.setGridTypeToUniform();
    object.setVeclen( veclen );
    object.setResolution( resolution );
    object.setValues( values );
    object.setCoords( coords );
    return object;
}
} // end of namespace local
