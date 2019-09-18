#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <YinYangVis/Lib/YinYangVolumeObject.h>

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
    kvs::ValueArray<float> sgrid__cosphi, sgrid__costht, sgrid__sinphi, sgrid__sintht;
    kvs::ValueArray<float> sgrid__rad, sgrid__theta, sgrid__phi;

    kvs::ValueArray<float> sgrid__value;

    kvs::ValueArray<float> ogrid__rad, ogrid__theta, ogrid__phi;

    kvs::AnyValueArray  ogrid__value;

    
    float sgrid__drad, sgrid__dtht, sgrid__dphi,
      sgrid__rad_min, sgrid__rad_max,
      sgrid__tht_min, sgrid__tht_max,
      sgrid__phi_min, sgrid__phi_max,
      SGRID__THETA_FROM, SGRID__THETA_TO,
      SGRID__PHI_FROM, SGRID__PHI_TO,
      
      ogrid__drad, ogrid__dtht, ogrid__dphi,
      /*ogrid__rad_max, ogrid__rad_min,  
      ogrid__tht_max, ogrid__tht_min,          
      ogrid__phi_max, ogrid__phi_min,      */               
      OGRID__THETA_FROM, OGRID__THETA_TO,
      OGRID__PHI_FROM, OGRID__PHI_TO;
    
    YinYangVis::YinYangVolumeObject::Range ogrid__range_r;
    YinYangVis::YinYangVolumeObject::Range ogrid__range_theta;
    YinYangVis::YinYangVolumeObject::Range ogrid__range_phi;
    
    
  public:
    Sgrid( const YinYangVis::YinYangVolumeObject& object);
    void sgrid__localize();
    void sgrid__make();
    void sgrid__output();
    void sgrid__rtp2xyz();
    void sgrid__xyz2rtp();
    void set_minmax();
    void set_nrtp();
    void set_drtp();
    void set_metric();

    void ogrid__make( const YinYangVis::YinYangVolumeObject& object );
    void set_o_minmax( const YinYangVis::YinYangVolumeObject& object );
    void set_o_nrtp( const YinYangVis::YinYangVolumeObject& object );
    void set_o_metric( const YinYangVis::YinYangVolumeObject& object );

    void mapping__localize();
    void iFind(float rad, float tht, float phi, int i, int j, int k );
    void sgrid__rtp2xyz ( float rad, float tht, float phi, float cart[3] );
    void sgrid__xyz2rtp(float cart[3], float polar[3]);
    void ogrid__find_near_corner(float rad,float theta,float phi);
    void ogrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1);
  };
}  // end of namespace local

