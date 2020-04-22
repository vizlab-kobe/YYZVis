#include "sgrid.h"
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <vector>
#include <math.h>

namespace local
{
  Sgrid::Sgrid( const YYZVis::YinYangVolumeObject& yin_volume, const YYZVis::YinYangVolumeObject& yang_volume, const YYZVis::ZhongVolumeObject& zhong_volume )
  {
    this->sgrid__make( yin_volume );
    this->mapping__localize( yin_volume, yang_volume, zhong_volume );
  }

      
  void Sgrid::sgrid__make( const YYZVis::YinYangVolumeObject& yoy_object )
  {
    this->set_minmax( yoy_object );
    this->set_rtp( yoy_object );
    this->set_metric();
  }
 
  void Sgrid::set_minmax( const YYZVis::YinYangVolumeObject& yoy_object )
  {
    sgrid__rad_min = yoy_object.rangeR().d / 2;
    sgrid__rad_max = yoy_object.rangeR().max + yoy_object.rangeR().d;
 
    sgrid__tht_min = yoy_object.rangeTheta().d / 2;
    sgrid__tht_max = M_PI - yoy_object.rangeTheta().d / 2;

    //sgrid__phi_min = yoy_object.rangeTheta().d / 2;
    //sgrid__phi_max = 2*M_PI - yoy_object.rangeTheta().d / 2;
    
    sgrid__phi_min = -M_PI + yoy_object.rangeTheta().d / 2;
    sgrid__phi_max = M_PI - yoy_object.rangeTheta().d / 2;
  }
    
  void Sgrid::set_rtp( const YYZVis::YinYangVolumeObject& yoy_object )
  {
    float rad_culling, tht_culling, phi_culling;
       
    rad_culling = 2.5;
    tht_culling = 2.5;
    phi_culling = 2.5;

    sgrid__drad = rad_culling * yoy_object.rangeR().d;
    sgrid__dtht = tht_culling * yoy_object.rangeTheta().d;
    sgrid__dphi = phi_culling * yoy_object.rangePhi().d;

    sgrid__size.nr = ( sgrid__rad_max - sgrid__rad_min ) / sgrid__drad + 1;
    sgrid__size.nt = ( sgrid__tht_max - sgrid__tht_min ) / sgrid__dtht + 1;
    sgrid__size.np = ( sgrid__phi_max - sgrid__phi_min ) / sgrid__dphi + 1;

    sgrid__rad.allocate(sgrid__size.nr);
    sgrid__theta.allocate(sgrid__size.nt);
    sgrid__phi.allocate(sgrid__size.np);
   
    sgrid__values.reserve(sgrid__size.nr*sgrid__size.nt*sgrid__size.np);
    sgrid__coords.reserve(3*sgrid__size.nr*sgrid__size.nt*sgrid__size.np);
    
    for(int i = 0; i<sgrid__size.nr*sgrid__size.nt*sgrid__size.np; i++)
      {
	sgrid__values.push_back(0);
      }

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
		 sgrid__coords.push_back(sgrid__theta[j]);
		 sgrid__coords.push_back(sgrid__phi[k]);
	       }
	   }
       }
  }

  
  void Sgrid::mapping__localize( const YYZVis::YinYangVolumeObject& yin_volume, const YYZVis::YinYangVolumeObject& yang_volume, const YYZVis::ZhongVolumeObject& zhong_volume )
  {
    int i, j, k;
    float cart[3];     //{ x, y, z }
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
		this->sgrid__rtp2xyz( rad, tht, phi, cart );

		 if( rad > yin_volume.rangeR().max ) continue;
                    
		 if( rad <= yin_volume.rangeR().min )
		    {
		      // ZHONG
	              this->iFind_zhong( cart[0], cart[1], cart[2], index, zhong_volume );
	              continue;
		    }
		 if( tht >= yin_volume.rangeTheta().min && tht <= yin_volume.rangeTheta().max &&
		     phi >= yin_volume.rangePhi().min && phi <= yin_volume.rangePhi().max )
	            {
	              // YIN
		      this->iFind( cart[0], cart[1], cart[2], index, yin_volume );
	            }
		else
		    {
		      // YANG
		      this->iFind( cart[0], cart[1], cart[2], index, yang_volume );
		    }
		index++;
	      }	 
	  } 
      }
  }

  void Sgrid::iFind(float x, float y, float z, int index, const YYZVis::YinYangVolumeObject& yoy_object )
  {
    float polar[3];    //{ r, t, p }
    
    this->sgrid__xyz2rtp( x, y, z, polar, yoy_object );
    this->ogrid__find_near_corner( polar[0], polar[1], polar[2], index, yoy_object);
      }

  void Sgrid::iFind_zhong(float x, float y, float z, int index, const YYZVis::ZhongVolumeObject& z_object )
  {
    int i1,j1,k1;
    float wx1,wy1,wz1;
    float x_next,y_next,z_next;
    float dx;
   
	i1 =  igrid__find_nearleft('x', x, z_object);
	j1 =  igrid__find_nearleft('y', y, z_object);
	k1 =  igrid__find_nearleft('z', z, z_object);

        x_next =  -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(i1+1);
	y_next =  -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(j1+1);
	z_next =  -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(k1+1);
	dx = (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1);
	
	wx1 = ( x_next - x ) / dx;
	wy1 = ( y_next - y ) / dx;   //dx = dy = dz
	wz1 = ( z_next - z ) / dx;   //dx = dy = dz

  	this->igrid_to_sgrid_localize(i1, j1, k1, wx1, wy1, wz1, index, z_object);      

      }

  void Sgrid::sgrid__rtp2xyz( float rad, float tht, float phi, float cart[3] )
  {
    cart[0] = rad*sin(tht)*cos(phi);
    cart[1] = rad*sin(tht)*sin(phi);
    cart[2] = rad*cos(tht);
  }

  void Sgrid::sgrid__xyz2rtp( float x, float y, float z, float polar[3], const YYZVis::YinYangVolumeObject& object )
  {
    float r, t, p;
   
    r = sqrt(  x*x + y*y + z*z );
 
    if ( object.gridType() == YYZVis::YinYangVolumeObject::Yin )
      {
	//YIN
        t = acos( z/r );
	p = atan2( y, x );
      }
    if( object.gridType() == YYZVis::YinYangVolumeObject::Yang)
      {
	//YANG
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

  int Sgrid::igrid__find_nearleft( char axis, float c, const YYZVis::ZhongVolumeObject& object )
  {
    int i;
    float x, y, z;
    float x_next, y_next, z_next;
    
    switch (axis){
    case 'x':
      for ( i = 0; i < object.dim(); i++ )
      {
	x =  -( object.rangeR().min + object.rangeR().d * 2 ) + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1) * i;
	  x_next = x + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1);
	if ( c >= x && c <= x_next ) 
	  {
	    return i;
	  }
      }
       break;
    case 'y':
      for ( i = 0; i < object.dim(); i++ )
      {
	y =  -( object.rangeR().min + object.rangeR().d * 2 ) + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1) * i;
	y_next = y + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1);
	if ( c >= y && c <= y_next ) 
	  {
	    return i;
	    break;
	  }
      }
      break;
    case 'z':
      for ( i = 0; i < object.dim(); i++ )
      {
	z =  -( object.rangeR().min + object.rangeR().d * 2 ) + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1) * i;
	z_next = z + (object.rangeR().min+object.rangeR().d*2)*2/(object.dim()-1);
	if ( c >= z && c <= z_next ) 
	  {
	    return i;
	    break;
	  }
      }
      break;
    }
    return 0;
  }

  void Sgrid::ogrid__find_near_corner( float rad, float theta, float phi, int index, const YYZVis::YinYangVolumeObject& object  )
  {                            
    int i1=0, j1=0, k1=0;                         
    float wr1=0, wt1=0, wp1=0;          
    int i, j, k;
    float ogrid_rad, ogrid_theta, ogrid_phi;
    float ogrid_rad_next, ogrid_theta_next, ogrid_phi_next;

    for ( i = 0; i < object.dimR() - 2; i++ )
      {
	ogrid_rad = object.rangeR().min + object.rangeR().d * i;
	ogrid_rad_next = ogrid_rad + object.rangeR().d;
	if ( rad >= ogrid_rad && rad <= ogrid_rad_next ) 
	  {
	    i1 = i;
	    wr1 = ( ogrid_rad_next - rad ) / object.rangeR().d;
	    break;
	  }
      }

    for ( j = 0; j < object.dimTheta() - 2; j++ )
      {
	ogrid_theta = object.rangeTheta().min + object.rangeTheta().d * j;
	ogrid_theta_next = ogrid_theta + object.rangeTheta().d;
	if ( theta >= ogrid_theta && theta <= ogrid_theta_next ) 
	  {
	    j1 = j;
	    wt1 = ( ogrid_theta_next - theta ) / object.rangeTheta().d;
	    break;
	  } 
      }

    for ( k = 0; k < object.dimPhi() - 2; k++ )
      {
	ogrid_phi = object.rangePhi().min + object.rangePhi().d * k;
	ogrid_phi_next = ogrid_phi + object.rangePhi().d;
	if ( phi >= ogrid_phi && phi <= ogrid_phi_next ) 
	  {
	    k1 = k;
	    wp1 = ( ogrid_phi_next - phi ) / object.rangePhi().d;
	    break;
	  }
      }
   
    this-> ogrid_to_sgrid_localize( i1, j1, k1, wr1, wt1, wp1, rad, theta, phi,index, object );
      }
  
  void Sgrid::ogrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index, const YYZVis::YinYangVolumeObject& object )
   {
     float wr2, wt2, wp2;
     float w[8];
     float v[8];
     float value;
     size_t a = k1*object.dimR()*object.dimTheta() + j1*object.dimR() + i1;

     wr2 = 1 - wr1;
     wt2 = 1 - wt1;
     wp2 = 1 - wp1;

     v[0] = object.values().at<float>( a );   
     v[1] = object.values().at<float>( a + 1 );
     v[2] = object.values().at<float>( a + object.dimR() );
     v[3] = object.values().at<float>( a + object.dimR() + 1 );
     v[4] = object.values().at<float>( a + (object.dimR()*object.dimTheta()) );
     v[5] = object.values().at<float>( a + (object.dimR()*object.dimTheta()) + 1 );
     v[6] = object.values().at<float>( a + object.dimR() + (object.dimR()*object.dimTheta()) );
     v[7] = object.values().at<float>( a + object.dimR() + (object.dimR()*object.dimTheta()) + 1 );

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

  void Sgrid::igrid_to_sgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index, const YYZVis::ZhongVolumeObject& object)
   {
     float wr2, wt2, wp2;
     float w[8];
     float v[8];
     float value;
     size_t a = k1*object.dim()*object.dim() +j1*object.dim() + i1 ;
     
     wr2 = 1 - wr1;
     wt2 = 1 - wt1;
     wp2 = 1 - wp1;

     v[0] = object.values().at<float>( a );   
     v[1] = object.values().at<float>( a + 1 );
     v[2] = object.values().at<float>( a + object.dim() );
     v[3] = object.values().at<float>( a + object.dim() + 1 );
     v[4] = object.values().at<float>( a + (object.dim()*object.dim()) );
     v[5] = object.values().at<float>( a + (object.dim()*object.dim()) + 1 );
     v[6] = object.values().at<float>( a + object.dim() + (object.dim()*object.dim()) );
     v[7] = object.values().at<float>( a + object.dim() + (object.dim()*object.dim()) + 1 );

     
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
