#include <kvs/ValueArray>
#include <YinYangVis/Lib/YinYangVolumeObject.h>

namespace local
{
class Sgrid
{
public:
    struct Sgrid__size{
        int nr, nt, np;
    };
    
public:
    kvs::ValueArray<float> sgrid__cosphi, sgrid__costht, sgrid__sinphi, sgrid__sintht;
    kvs::ValueArray<float> sgrid__rad, sgrid__theta, sgrid__phi;
    
    float sgrid__drad, sgrid__dtht, sgrid__dphi,
          sgrid__rad_min, sgrid__rad_max,
          sgrid__tht_min, sgrid__tht_max,
          sgrid__phi_min, sgrid__phi_max,
          SGRID__THETA_FROM, SGRID__THETA_TO,
          SGRID__PHI_FROM, SGRID__PHI_TO,
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
    void set_ogrid( const YinYangVis::YinYangVolumeObject& object );
};
}  // end of namespace local

