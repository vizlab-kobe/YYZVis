#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <kvs/StructuredVolumeObject>


namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinYangVolumeObject& yin_volume,
    YYZVis::YinYangVolumeObject& yang_volume,
    YYZVis::ZhongVolumeObject& zhong_volume );

} // end of namespace local
