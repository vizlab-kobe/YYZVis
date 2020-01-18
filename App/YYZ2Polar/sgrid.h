#include <kvs/ValueArray>
#include <vector>
#include <kvs/AnyValueArray>
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>

namespace local
{
  class Sgrid
  {
  public:
    struct Sgrid__size{
      int nr, nt, np;
    };
    Sgrid__size sgrid__size;

  public:
    float yoy_min_value=0,yoy_max_value=0;
     float zhong_min_value=0,zhong_max_value=0;
    kvs::ValueArray<float> sgrid__rad, sgrid__theta, sgrid__phi;
    std::vector<float> sgrid__coords;
    std::vector<float> sgrid__values;

    float sgrid__drad, sgrid__dtht, sgrid__dphi,
          sgrid__rad_min, sgrid__rad_max,
          sgrid__tht_min, sgrid__tht_max,
          sgrid__phi_min, sgrid__phi_max;
          
  public:
    Sgrid( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume  );
    void sgrid__make( const YinYangVis::YinYangVolumeObject& yoy_object );
    void sgrid__output();
    void set_minmax( const YinYangVis::YinYangVolumeObject& yoy_object );
    void set_rtp( const YinYangVis::YinYangVolumeObject& yoy_object );
    void set_metric( );

    void mapping__localize( const YinYangVis::YinYangVolumeObject& yin_volume,  const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume );
    void iFind( float rad, float tht, float phi, int index, const YinYangVis::YinYangVolumeObject& yoy_object);
    void iFind_zhong(float rad, float tht, float phi, int index, const YinYangVis::ZhongVolumeObject& z_object );
    void sgrid__rtp2xyz ( float rad, float tht, float phi, float cart[3] );
    void sgrid__xyz2rtp(float cart[3], float polar[3], const YinYangVis::YinYangVolumeObject& object);
    void ogrid__find_near_corner( float rad, float theta, float phi, int index, const YinYangVis::YinYangVolumeObject& object );
    int  igrid__find_nearleft( char axis, float c, const YinYangVis::ZhongVolumeObject& object );
    void ogrid_to_sgrid_localize( int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index, const YinYangVis::YinYangVolumeObject& object );
    void igrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index, const YinYangVis::ZhongVolumeObject& object );
  };
}  // end of namespace local

