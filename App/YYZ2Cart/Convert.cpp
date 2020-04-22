#include "Convert.h"
#include "cgrid.h"

namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinYangVolumeObject& yin_volume,
    YYZVis::YinYangVolumeObject& yang_volume,
    YYZVis::ZhongVolumeObject& zhong_volume )
{
    // Vector length
    const int veclen = 1;
    int dim_x, dim_y, dim_z;

    local::Cgrid cgrid( yin_volume, yang_volume, zhong_volume );
    
    dim_x = cgrid.cgrid__size.nx;
    dim_y = cgrid.cgrid__size.ny;
    dim_z = cgrid.cgrid__size.nz;
    
    // Grid resolution (rx, ry, rz)
    const kvs::Vec3u resolution( dim_x, dim_y, dim_z );

    // Values (physical quantity)
    kvs::ValueArray<kvs::Real32> values( cgrid.cgrid__values );

    // Coordinate value
    kvs::ValueArray<kvs::Real32> coords( cgrid.cgrid__coords );

    kvs::StructuredVolumeObject object;
    object.setGridTypeToUniform();
    object.setVeclen( veclen );
    object.setResolution( resolution );
    object.setValues( values );
    object.setCoords( coords );

    return object;
}
} // end of namespace local
