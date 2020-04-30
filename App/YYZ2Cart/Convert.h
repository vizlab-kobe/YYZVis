#include <YYZVis/Lib/YinVolumeObject.h>
#include <YYZVis/Lib/YangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <kvs/StructuredVolumeObject>


namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinVolumeObject& yin_volume,
    YYZVis::YangVolumeObject& yng_volume,
    YYZVis::ZhongVolumeObject& zng_volume );

} // end of namespace local
