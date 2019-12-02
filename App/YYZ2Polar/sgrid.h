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

    struct Ogrid__size{
      int nr, nt, np;
    };
    Ogrid__size ogrid__size;
     
  public:
    
    kvs::ValueArray<float> sgrid__rad, sgrid__theta, sgrid__phi;
    std::vector<float> sgrid__coords;
    std::vector<float> sgrid__values;
    kvs::ValueArray<float> ogrid__rad, ogrid__theta, ogrid__phi;
    kvs::AnyValueArray  ogrid__values;
    kvs::ValueArray<float> igrid__x, igrid__y, igrid__z; 
    kvs::AnyValueArray  igrid__values;

    float sgrid__drad, sgrid__dtht, sgrid__dphi,
      sgrid__rad_min, sgrid__rad_max,
      sgrid__tht_min, sgrid__tht_max,
      sgrid__phi_min, sgrid__phi_max,
      SGRID__THETA_FROM, SGRID__THETA_TO,
      SGRID__PHI_FROM, SGRID__PHI_TO,
      
      ogrid__drad, ogrid__dtht, ogrid__dphi,
              
      OGRID__THETA_FROM, OGRID__THETA_TO,
      OGRID__PHI_FROM, OGRID__PHI_TO,

      igrid__dx, igrid__dy, igrid__dz,
      igrid__dim;
    
    YinYangVis::YinYangVolumeObject::Range ogrid__range_r;
    YinYangVis::YinYangVolumeObject::Range ogrid__range_theta;
    YinYangVis::YinYangVolumeObject::Range ogrid__range_phi;
    
    
  public:
    Sgrid( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume  );
    void sgrid__make();
    void sgrid__output();
    void sgrid__rtp2xyz();
    void sgrid__xyz2rtp();
    void set_minmax();
    void set_rtp();
    void set_metric();

    void ogrid__make( const YinYangVis::YinYangVolumeObject& yoy_object );
    void set_o_minmax( const YinYangVis::YinYangVolumeObject& yoy_object );
    void set_o_nrtp( const YinYangVis::YinYangVolumeObject& yoy_object );
    void set_o_metric( const YinYangVis::YinYangVolumeObject& yoy_object );

    void igrid__make( const YinYangVis::ZhongVolumeObject& z_object ); 

    void mapping__localize();
    void iFind( float rad, float tht, float phi, int index );
    void sgrid__rtp2xyz ( float rad, float tht, float phi, float cart[3] );
    void sgrid__xyz2rtp(float cart[3], float polar[3]);
    void ogrid__find_near_corner( float rad, float theta, float phi, int index );
    int  igrid__find_nearleft( char axis, float c );
    void ogrid_to_sgrid_localize( int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index );
    void igrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index );
  };
}  // end of namespace local

