#include "sgrid.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <vector>
#include <math.h>

namespace local
{
  Sgrid::Sgrid( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume )
  {
    this->ogrid__make( yin_volume );
    this->igrid__make( zhong_volume ); 
    this->sgrid__make();
    this->sgrid__localize();
    this->mapping__localize();
    this->ogrid__make( yang_volume );
    this->mapping__localize();
  }

  void Sgrid::ogrid__make( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    this->set_o_minmax( yoy_object );
    this->set_o_nrtp( yoy_object );
    this->set_o_metric( yoy_object );
    ogrid__values = yoy_object.values();
  }

   void Sgrid::igrid__make( const YinYangVis::ZhongVolumeObject& z_object  )
   {
     float igrid__xmin, igrid__xmax;
     float igrid__ymin, igrid__ymax;
     float igrid__zmin, igrid__zmax;
     
     igrid__dim = z_object.dim();

     igrid__x.allocate(igrid__dim);
     igrid__y.allocate(igrid__dim);
     igrid__z.allocate(igrid__dim);
     
     igrid__xmin = -z_object.rangeR().max;
     igrid__xmax = z_object.rangeR().max;
     igrid__ymin = -z_object.rangeR().max;
     igrid__ymax = z_object.rangeR().max;
     igrid__zmin = -z_object.rangeR().max;
     igrid__zmax = z_object.rangeR().max;
   
     igrid__dx = ( igrid__xmax*2 ) / (igrid__dim-1);
     igrid__dy = ( igrid__xmax*2 ) / (igrid__dim-1);
     igrid__dz = ( igrid__xmax*2 ) / (igrid__dim-1);

     igrid__values = z_object.values();
     
for ( int i = 0; i < igrid__dim; i++ )
      {
	igrid__x[i] = igrid__xmin + igrid__dx * i;
	igrid__y[i] = igrid__ymin + igrid__dy * i;
	igrid__z[i] = igrid__zmin + igrid__dz * i;
      }
  
   }

  void Sgrid::set_o_minmax( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    OGRID__THETA_FROM = M_PI/4;
    OGRID__THETA_TO   = M_PI - M_PI/4;
    OGRID__PHI_FROM   = -(3*M_PI)/4.0;
    OGRID__PHI_TO     =  (3*M_PI)/4.0;// この辺全部いらない
        
    ogrid__range_r = yoy_object.rangeR();
    ogrid__range_theta = yoy_object.rangeTheta();
    ogrid__range_phi = yoy_object.rangePhi();
  }

  void Sgrid::set_o_nrtp( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    ogrid__size.nr = yoy_object.dimR();
    ogrid__size.nt = yoy_object.dimTheta();
    ogrid__size.np = yoy_object.dimPhi();
  
      ogrid__rad.allocate(ogrid__size.nr);
      ogrid__theta.allocate(ogrid__size.nt);
      ogrid__phi.allocate(ogrid__size.np);
      //  ogrid__value.allocate( ogrid__size.nr * ogrid__size.nt * ogrid__size.np );
  }

  void Sgrid::set_o_metric( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    int i, j, k;
      
    for ( i = 0; i < ogrid__size.nr; i++ )
      {
	ogrid__rad[i] = yoy_object.rangeR().min + ogrid__range_r.d * i;
      }
     
    for ( j = 0; j < ogrid__size.nt; j++ )
      {
	ogrid__theta[j] = yoy_object.rangeTheta().min + ogrid__range_theta.d * j; 
      }
     
    for ( k = 0; k < ogrid__size.np; k++ )
      {
	ogrid__phi[k] = yoy_object.rangePhi().min + ogrid__range_phi.d * k;
      }
  }
  
  void Sgrid::set_minmax()
  {
    SGRID__THETA_FROM = 0;
    SGRID__THETA_TO   = M_PI;
    SGRID__PHI_FROM   = -M_PI;
    SGRID__PHI_TO     =  M_PI;
        
    sgrid__rad_min = ogrid__range_r.d / 2;
    sgrid__rad_max = ogrid__range_r.max + ogrid__range_r.d;
        
    // Changed by Yoshiki Ueda on 2017.06.23
    //sgrid__tht_min = SGRID__THETA_FROM + ogrid__dtht / 2;
    //sgrid__tht_max = SGRID__THETA_TO - ogrid__dtht / 2;
    sgrid__tht_min = SGRID__THETA_FROM + ogrid__range_theta.d / 10000;
    sgrid__tht_max = SGRID__THETA_TO - ogrid__range_theta.d / 10000;
        
    sgrid__phi_min = SGRID__PHI_FROM + ogrid__range_phi.d / 2;
    // Changed by Yoshiki Ueda on 2017.06.23
    // sgrid__phi_max = SGRID__PHI_TO - ogrid__dphi / 2;
    sgrid__phi_max = SGRID__PHI_TO + ogrid__range_phi.d / 2;

  }
    
  void Sgrid::set_nrtp()
  {
    float drad, dtht, dphi;
    float rad_culling, tht_culling, phi_culling;
        
    rad_culling = 2.5; //namelist__double('DIM1_culling')
    tht_culling = 2.5; //namelist__double('DIM2_culling')
    phi_culling = 2.5; //namelist__double('DIM3_culling')
        
    /* if(rad_culling < 1.0_DP ||
       tht_culling < 1.0_DP ||
       phi_culling < 1.0_DP) then
       call ut__fatal("< sgrid > Error : Culling values must be bigger than 1.0")
       end if*/
        
    drad = rad_culling * ogrid__range_r.d;
    dtht = tht_culling * ogrid__range_theta.d;
    dphi = phi_culling * ogrid__range_phi.d;
        
    sgrid__size.nr = ( sgrid__rad_max - sgrid__rad_min ) / drad + 1;
    sgrid__size.nt = ( sgrid__tht_max - sgrid__tht_min ) / dtht + 1;
    sgrid__size.np = ( sgrid__phi_max - sgrid__phi_min ) / dphi + 1;
        
    if ( sgrid__size.nt%2 == 1 )
      {
	sgrid__size.nt = sgrid__size.nt + 1;
      }

    if ( sgrid__size.np%2 == 1 )
      {
	sgrid__size.np = sgrid__size.np + 1;
      }

    sgrid__rad.allocate(sgrid__size.nr);
    sgrid__theta.allocate(sgrid__size.nt);
    sgrid__phi.allocate(sgrid__size.np);
    sgrid__sintht.allocate(sgrid__size.nt);
    sgrid__sinphi.allocate(sgrid__size.np);
    sgrid__costht.allocate(sgrid__size.nt);
    sgrid__cosphi.allocate(sgrid__size.np);
   
    sgrid__values.reserve(sgrid__size.nr*sgrid__size.nt*sgrid__size.np);
    sgrid__coords.reserve(3*sgrid__size.nr*sgrid__size.nt*sgrid__size.np);
    
    for(int i = 0; i<sgrid__size.nr*sgrid__size.nt*sgrid__size.np; i++)
      {
	sgrid__values.push_back(0);
      }
    //    sgrid__value.allocate(sgrid__size.nr * sgrid__size.nt * sgrid__size.np);
  }
    
  void Sgrid::set_drtp()
  {
    sgrid__drad = ( sgrid__rad_max - sgrid__rad_min ) / (sgrid__size.nr - 1);
    sgrid__dtht = ( sgrid__tht_max - sgrid__tht_min ) / (sgrid__size.nt - 1);
    sgrid__dphi = ( sgrid__phi_max - sgrid__phi_min ) / (sgrid__size.np - 1);

  }
    
  void Sgrid::set_metric()
  {
    int i, j, k;
    
    for ( i = 0; i < sgrid__size.nr; i++ )
      {
	sgrid__rad[i] = sgrid__rad_min + sgrid__drad * i;
      }
    for ( j = 0; j < sgrid__size.nt; j++ )
      {
	sgrid__theta[j] = sgrid__tht_min + sgrid__dtht * j;
      }
    for ( k = 0; k < sgrid__size.np; k++ )
      {
	sgrid__phi[k] = sgrid__phi_min + sgrid__dphi * k;
      }

     for ( i = 0; i < sgrid__size.nr; i++ )
       {
	 for ( j = 0; j < sgrid__size.nt; j++ )
	   {
	     for ( k = 0; k < sgrid__size.np; k++ )
	       {
		 sgrid__coords.push_back(sgrid__rad[i]);
		 sgrid__coords.push_back(sgrid__rad[j]);
		 sgrid__coords.push_back(sgrid__rad[k]);
	       }
	   }

       }
  }
    
  void Sgrid::sgrid__make()
  {
    this->set_minmax();
    this->set_nrtp();
    this->set_drtp();
    this->set_metric();
  }
  
  void Sgrid::sgrid__localize()
  {
    int j, k;
    float tht, phi;
        
    for ( j = 0; j < sgrid__size.nt; j++)
      {
	tht = sgrid__theta[j];
	sgrid__sintht[j] = sin(tht);
	sgrid__costht[j] = cos(tht);
      }
        
    for ( k = 0; k < sgrid__size.np; k++)
      {
	phi = sgrid__phi[k];
	sgrid__sinphi[k] = sin(phi);
	sgrid__cosphi[k] = cos(phi);
      }
	
  }

  void Sgrid::mapping__localize()
  {
    int i, j, k;
    float rad, tht, phi;
    int index = 0;
    for( k = 0; k < sgrid__size.np; k++ )
      {
        phi = sgrid__phi[k];
	for( j = 0; j < sgrid__size.nt; j++ )
	  {
	    tht = sgrid__theta[j];
	    for( i = 0; i < sgrid__size.nr; i++ )
	      {
		rad = sgrid__rad[i];
	        this->iFind( rad, tht, phi, index );
		index++;
	      }
	  }
      }
  }

  void Sgrid::iFind(float rad, float tht, float phi, int index )
  {
    float cart[3];     //{ x, y, z }
    int i1, j1, k1;

    float polar[3];    //{ r, t, p }
    float wx1, wy1, wz1;


    this->sgrid__rtp2xyz( rad, tht, phi, cart );
   
    if ( cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2] >= ogrid__range_r.max )
      {
	return;
	  }
    else if ( cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2] <= ogrid__range_r.min )
      {
	// Zhong
	i1 =  igrid__find_nearleft('x', cart[0]);
	j1 =  igrid__find_nearleft('y', cart[1]);
	k1 =  igrid__find_nearleft('z', cart[2]);

	wx1 = ( igrid__x[i1 + 1] - cart[0] ) / igrid__dx;
	wy1 = ( igrid__y[j1 + 1] - cart[1] ) / igrid__dy;
	wz1 = ( igrid__z[k1 + 1] - cart[2] ) / igrid__dz;

	this->igrid_to_sgrid_localize(i1, j1, k1, wx1, wy1, wz1, index);      
	return;
	  }
    
    this->sgrid__xyz2rtp( cart, polar );
   
    this->ogrid__find_near_corner( polar[0], polar[1], polar[2], index);
      }
 
  void Sgrid::sgrid__rtp2xyz( float rad, float tht, float phi, float cart[3] )
  {
    cart[0] = rad*sin(tht)*cos(phi);
    cart[1] = rad*sin(tht)*sin(phi);
    cart[2] = rad*cos(tht);
  }

  void Sgrid::sgrid__xyz2rtp( float cart[3],float polar[3] )
  {
    float x, y, z;
    float r, t, p;

    std::string flag;

    x = cart[0];
    y = cart[1];
    z = cart[2];
  
    r = sqrt(  x*x + y*y + z*z );
 
    t = acos( z/r );

    flag = "YANG";   //default

    if ( t >= OGRID__THETA_FROM && t <= OGRID__THETA_TO )
      {
	p = atan2( y, x );
	if ( p >= OGRID__PHI_FROM && p <= OGRID__PHI_TO )
	  {
	    flag = "YIN";
	  }
      }

    if ( flag ==  "YIN" )
      {
	float x_ = -x;
	float y_ = z;
	float z_ = y;
	t = acos( z_/r );
	p = atan2( y_, x_ );
      
      }
    
    polar[0] = r;
    polar[1] = t;
    polar[2] = p;

  }

  int Sgrid::igrid__find_nearleft( char axis, float c )
  {
    int i;
    
    switch (axis){
    case 'x':
       for ( i = 0; i < igrid__dim; i++ )
      {
	if ( c >= igrid__x[i] && c <= igrid__x[i+1] ) 
	  {
	    return i;
	  }
      }
       break;
    case 'y':
        for ( i = 0; i < igrid__dim; i++ )
      {
	if ( c >= igrid__y[i] && c <= igrid__y[i+1] ) 
	  {
	    return i;
	    break;
	  }
      }
      break;
    case 'z':
        for ( i = 0; i < igrid__dim; i++ )
      {
	if ( c >= igrid__z[i] && c <= igrid__z[i+1] ) 
	  {
	    return i;
	    break;
	  }
      }
      break;
    }
  }

  void Sgrid::ogrid__find_near_corner( float rad,float theta,float phi,int index  )
  {                            
    int i1, j1, k1;                         
    float wr1, wt1, wp1;          
    int i, j, k; 
    int nr, nt, np;
      //real(DP), parameter :: ERROR_MARK = -2007.0605_DP

    nr = ogrid__size.nr;
    nt = ogrid__size.nt;
    np = ogrid__size.np;

    for ( i = 0; i < nr-2; i++ )
      {
	if ( rad >= ogrid__rad[i] && rad <= ogrid__rad[i+1] ) 
	  {
	    i1 = i;
	    wr1 = (ogrid__rad[i+1]-rad)/ogrid__range_r.d;
	    break;
	  }
      }

    for ( j = 0; j < nt-2; j++ )
      {
	if ( theta >= ogrid__theta[j] && theta <= ogrid__theta[j+1] ) 
	  {
	    j1 = j;
	    wt1 = (ogrid__theta[j+1]-theta)/ogrid__range_theta.d;
	    break;
	  }
      }

    for ( k = 0; k < np-2; k++ )
      {
	if ( phi >= ogrid__phi[k] && phi <= ogrid__phi[k+1] ) 
	  {
	    k1 = k;
	    wp1 = (ogrid__phi[k+1]-phi)/ogrid__range_phi.d;
	    break;
	  }
      }
   
    this-> ogrid_to_sgrid_localize( i1, j1, k1, wr1, wt1, wp1, rad, theta, phi,index );
   
      }
  void Sgrid::ogrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index )
   {
     float wr2, wt2, wp2;
     float w[8];
     float v[8];
     float value;
     size_t a = k1*ogrid__size.nr*ogrid__size.nt + j1*ogrid__size.nr + i1 ;
     
     wr2 = 1 - wr1;
     wt2 = 1 - wt1;
     wp2 = 1 - wp1;

     v[0] = ogrid__values.at<float>( a );   
     v[1] = ogrid__values.at<float>( a + 1 );
     v[2] = ogrid__values.at<float>( a + ogrid__size.nr );
     v[3] = ogrid__values.at<float>( a + ogrid__size.nr + 1 );
     v[4] = ogrid__values.at<float>( a + (ogrid__size.nr*ogrid__size.nt) );
     v[5] = ogrid__values.at<float>( a + (ogrid__size.nr*ogrid__size.nt) + 1 );
     v[6] = ogrid__values.at<float>( a + ogrid__size.nr + (ogrid__size.nr*ogrid__size.nt) );
     v[7] = ogrid__values.at<float>( a + ogrid__size.nr + (ogrid__size.nr*ogrid__size.nt) + 1 );

     
     w[0] = wr1 * wt1 * wp1;
     w[1] = wr2 * wt1 * wp1;
     w[2] = wr1 * wt2 * wp1;
     w[3] = wr2 * wt2 * wp1;
     w[4] = wr1 * wt1 * wp2;
     w[5] = wr2 * wt1 * wp2;
     w[6] = wr1 * wt2 * wp2;
     w[7] = wr2 * wt2 * wp2;

     value = v[0] * w[0] + v[1] * w[1]
           + v[2] * w[2] + v[3] * w[3]
           + v[4] * w[4] + v[5] * w[5]
           + v[6] * w[6] + v[7] * w[7];

     sgrid__values[index] = value;
   
   }

  void Sgrid::igrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index )
   {
     float wr2, wt2, wp2;
     float w[8];
     float v[8];
     float value;
     size_t a = k1*igrid__dim*igrid__dim +j1*igrid__dim + i1 ;
     
     wr2 = 1 - wr1;
     wt2 = 1 - wt1;
     wp2 = 1 - wp1;

     v[0] = igrid__values.at<float>( a );   
     v[1] = igrid__values.at<float>( a + 1 );
     v[2] = igrid__values.at<float>( a + igrid__dim );
     v[3] = igrid__values.at<float>( a + igrid__dim + 1 );
     v[4] = igrid__values.at<float>( a + (igrid__dim*igrid__dim) );
     v[5] = igrid__values.at<float>( a + (igrid__dim*igrid__dim) + 1 );
     v[6] = igrid__values.at<float>( a + igrid__dim + (igrid__dim*igrid__dim) );
     v[7] = igrid__values.at<float>( a + igrid__dim + (igrid__dim*igrid__dim) + 1 );

     
     w[0] = wr1 * wt1 * wp1;
     w[1] = wr2 * wt1 * wp1;
     w[2] = wr1 * wt2 * wp1;
     w[3] = wr2 * wt2 * wp1;
     w[4] = wr1 * wt1 * wp2;
     w[5] = wr2 * wt1 * wp2;
     w[6] = wr1 * wt2 * wp2;
     w[7] = wr2 * wt2 * wp2;

     value = v[0] * w[0] + v[1] * w[1]
           + v[2] * w[2] + v[3] * w[3]
           + v[4] * w[4] + v[5] * w[5]
           + v[6] * w[6] + v[7] * w[7];

     sgrid__values[index] = value;
   
   }
}  // end of namespace local
