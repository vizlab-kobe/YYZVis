#include "sgrid.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <kvs/ValueArray>
//#include <bits/stdc++.h>
#include <math.h>

namespace local
{
  Sgrid::Sgrid( const YinYangVis::YinYangVolumeObject& object )
  {
    this->ogrid__make( object );
    this->sgrid__make();
    this->sgrid__localize();
    this->mapping__localize();
  }

  void Sgrid::ogrid__make( const YinYangVis::YinYangVolumeObject& object )
  {
    this->set_o_minmax( object );
    this->set_o_nrtp( object );
    //  this->set_o_drtp();
    this->set_o_metric();
    
  }
  void Sgrid::set_o_minmax( const YinYangVis::YinYangVolumeObject& object )
  {
    OGRID__THETA_FROM = M_PI/4;
    OGRID__THETA_TO   = M_PI - M_PI/4;
    OGRID__PHI_FROM   = -(3*M_PI)/4.0;
    OGRID__PHI_TO     =  (3*M_PI)/4.0;
        
    ogrid__range_r = object.rangeR();
    ogrid__range_theta = object.rangeTheta();
    ogrid__range_phi = object.rangePhi();
    
  }

  void Sgrid::set_o_nrtp( const YinYangVis::YinYangVolumeObject& object )
  {
    ogrid__size.nr = object.dimR();
    ogrid__size.nt = object.dimTheta();
    ogrid__size.np = object.dimPhi();
  }

  void Sgrid::set_o_metric()
  {
    int i, j, k;

    for ( i = 0; i < ogrid__size.nr; i++ )
      {
	ogrid__rad[i] = ogrid__range_r.min + ogrid__range_r.d * i;
      }
    for ( j = 0; j < ogrid__size.nt; j++ )
      {
	ogrid__theta[j] = OGRID__THETA_FROM + ogrid__range_theta.d * j; 
      }
    for ( k = 0; k < ogrid__size.np; k++ )
      {
	ogrid__phi[k] = OGRID__PHI_FROM + ogrid__range_phi.d * k;
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
    /*  call ut__message('<sgrid> sgrid__size%nr,nt,np = ', &
	sgrid__size%nr,              &
	sgrid__size%nt,              &
	sgrid__size%np)*/
  }
    
  void Sgrid::set_drtp()
  {
    sgrid__drad = ( sgrid__rad_max - sgrid__rad_min ) / sgrid__size.nr - 1;
    sgrid__dtht = ( sgrid__tht_max - sgrid__tht_min ) / sgrid__size.nt - 1;
    sgrid__dphi = ( sgrid__phi_max - sgrid__phi_min ) / sgrid__size.np - 1;
  }
    
  void Sgrid::set_metric()
  {
    int i, j, k;
    float tht, phi;
    //   size_t nr = sgrid__size.nr;
    //kvs::ValueArray<float> sgrid__rad( nr );
    //kvs::ValueArray<float> sgrid__theta( sgrid__size.nt );
	
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

    for( k = 0; k < sgrid__size.np; k++ )
      {
        phi = sgrid__phi[k];
	for( j = 0; j < sgrid__size.nt; j++ )
	  {
	    tht = sgrid__theta[j];
	    for( i = 0; i < sgrid__size.nr; i++ )
	      {
		rad = sgrid__rad[i];
	        this->iFind( rad, tht, phi, i, j, k );
	      }
	  }
      }
  }

  void Sgrid::iFind(float rad, float tht, float phi, int i, int j, int k )
  {
    float cart[3];     //{ x, y, z }
    int i1, j1, k1;

    float polar[3];    //{ r, t, p }
    float wr1, wt1, wp1;
    float wx1, wy1, wz1;


    this->sgrid__rtp2xyz( rad, tht, phi, cart );
    
    if ( cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2] >= ogrid__range_r.max )
      {
	return;
	  }
    else if ( cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2] <= ogrid__range_r.min )
      {
	// Zhong
	/*	i1 =  igrid__find_nearleft_int('x', x)
	  j1 =  igrid__find_nearleft_int('y', y)
	  k1 =  igrid__find_nearleft_int('z', z)

	  wx1 = ( igrid__pos_x(i1 + 1) - x ) / igrid__dx
	  wy1 = ( igrid__pos_y(j1 + 1) - y ) / igrid__dy
	  wz1 = ( igrid__pos_z(k1 + 1) - z ) / igrid__dz

	  this->igrid_to_sgrid_localize(i,  j,  k, i1, j1, k1, wx1,wy1,wz1);
	*/
	return;
	  }
    
    this->sgrid__xyz2rtp( cart, polar );

    //this->ogrid__find_near_corner(r,t,p,i1,j1,k1,wr1,wt1,wp1);
    // this->ogrid_to_sgrid_localize(i, j, k,i1, j1, k1,wr1,wt1,wp1);

      }
 
  void Sgrid::sgrid__rtp2xyz( float rad, float tht, float phi, float cart[3] )
  {
    cart[0] = rad*sin(tht)*cos(phi);
    cart[1] = rad*sin(tht)*sin(phi);
    cart[2] = rad*cos(tht);

  }

  void Sgrid::sgrid__xyz2rtp( float cart[3],float polar[3] )
  {
    float rad;
    float x, y, z;
    float r, t, p;

    std::string flag;

    x = cart[0];
    y = cart[1];
    z = cart[2];

    rad = sqrt( x*x + y*y + z*z );

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

  void Sgrid::ogrid__find_near_corner(float rad,float theta,float phi,int i1,int j1,int k1,float wr1,float wt1,float wp1)
  {                            
    //   integer,  intent(out) :: i1, j1, k1                                  
    //real(DP), intent(out) :: wr1, wt1, wp1             
     
    int nr, nt, np, i, j, k;
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
      }
      
}  // end of namespace local
