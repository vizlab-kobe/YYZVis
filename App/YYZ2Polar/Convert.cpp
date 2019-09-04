#include "Convert.h"
#include "Sgrid.h"

namespace local
{

kvs::StructuredVolumeObject Convert(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume )
{
    
    local::Sgrid sgrid( yin_volume );
    
    // Vector length
    //const int veclen = 1;

    // Grid resolution (rx, ry, rz)
    //const kvs::Vec3u resolution( rx, ry, rz );

    // Values (physical quantity)
    //kvs::ValueArray<kvs::Real32> values( rx, * ry * rz );

    // Coordinate value
    //kvs::ValueArray<kvs::Real32> coords( rx, * ry * rz * 3 );

    kvs::StructuredVolumeObject object;
    //object.setGridTypeToCurvilinear();
    //object.setVeclen( veclen );
    //object.setResolution( resolution );
    //object.setValues( values );
    return object;
}

    
} // end of namespace local
