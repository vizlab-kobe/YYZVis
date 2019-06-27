#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
#include <kvs/StructuredVolumeObject>


namespace local
{

kvs::StructuredVolumeObject Convert(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume );

} // end of namespace local
