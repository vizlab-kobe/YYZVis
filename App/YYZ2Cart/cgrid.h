#include <kvs/ValueArray>
#include <vector>
#include <kvs/AnyValueArray>
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>

namespace local
{
  class Cgrid
  {
  public:
    struct Cgrid__size{
      int nx, ny, nz;
    };
    Cgrid__size cgrid__size;

  public:
    kvs::ValueArray<float> cgrid__x, cgrid__y, cgrid__z;
    std::vector<float> cgrid__coords;
    std::vector<float> cgrid__values;

    float cgrid__dx, cgrid__dy, cgrid__dz,
          cgrid__x_min, cgrid__x_max,
          cgrid__y_min, cgrid__y_max,
          cgrid__z_min, cgrid__z_max;
          
  public:
    Cgrid( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume  );
    void cgrid__make();
    void cgrid__output();
    void set_minmax();
    void set_xyz();
    void set_metric();

    void mapping__localize( const YinYangVis::YinYangVolumeObject& yin_volume,  const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume );
    void iFind_zhong(float x, float y, float z, int index, const YinYangVis::ZhongVolumeObject& z_object );
    void ogrid__find_near_corner( float rad, float theta, float phi, int index, const YinYangVis::YinYangVolumeObject& object );
    int  igrid__find_nearleft( char axis, float c, const YinYangVis::ZhongVolumeObject& object );
    void ogrid_to_cgrid_localize( int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index, const YinYangVis::YinYangVolumeObject& object );
    void igrid_to_cgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index, const YinYangVis::ZhongVolumeObject& object );
  };
}  // end of namespace local

