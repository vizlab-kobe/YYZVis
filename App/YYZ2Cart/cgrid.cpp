#include "cgrid.h"
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <vector>
#include <math.h>

namespace local
{
  Cgrid::Cgrid( const YYZVis::YinYangVolumeObject& yin_volume, const YYZVis::YinYangVolumeObject& yang_volume, const YYZVis::ZhongVolumeObject& zhong_volume )
  {
    this->cgrid__make( yin_volume ); 
    this->mapping__localize( yin_volume, yang_volume, zhong_volume );
  }

  void Cgrid::cgrid__make( const YYZVis::YinYangVolumeObject& yin_volume )
  {
    this->set_minmax( yin_volume );
    this->set_xyz();
    this->set_metric();
  }  
  
  void Cgrid::set_minmax( const YYZVis::YinYangVolumeObject& yin_volume )
  {
    
    cgrid__x_min = -yin_volume.rangeR().max;
    cgrid__x_max =  yin_volume.rangeR().max;
    cgrid__y_min = -yin_volume.rangeR().max;
    cgrid__y_max =  yin_volume.rangeR().max;
    cgrid__z_min = -yin_volume.rangeR().max;
    cgrid__z_max =  yin_volume.rangeR().max;

  }
    
  void Cgrid::set_xyz()
  {   
    cgrid__size.nx = 200;
    cgrid__size.ny = 200;
    cgrid__size.nz = 200;
        
    cgrid__dx = (cgrid__x_max - cgrid__x_min) / ( cgrid__size.nx - 1);
    cgrid__dy = (cgrid__y_max - cgrid__y_min) / ( cgrid__size.ny - 1);
    cgrid__dz = (cgrid__z_max - cgrid__z_min) / ( cgrid__size.nz - 1);
    
    cgrid__x.allocate(cgrid__size.nx);
    cgrid__y.allocate(cgrid__size.ny);
    cgrid__z.allocate(cgrid__size.nz);
   
    cgrid__values.reserve(cgrid__size.nx*cgrid__size.ny*cgrid__size.nz);
    cgrid__coords.reserve(3*cgrid__size.nx*cgrid__size.ny*cgrid__size.nz);
    
    for(int i = 0; i<cgrid__size.nx*cgrid__size.ny*cgrid__size.nz; i++)
      {
	cgrid__values.push_back(0);
      }

  }
    
  void Cgrid::set_metric()
  {
    int i, j, k;

    for ( i = 0; i < cgrid__size.nx; i++ )
      {
	cgrid__x[i] = cgrid__x_min + cgrid__dx * i;
      }
    for ( j = 0; j < cgrid__size.ny; j++ )
      {
	cgrid__y[j] = cgrid__y_min + cgrid__dy * j;
      }
    for ( k = 0; k < cgrid__size.nz; k++ )
      {
	cgrid__z[k] = cgrid__z_min + cgrid__dz * k;
      }

    for ( i = 0; i < cgrid__size.nx; i++ )
      {
	for ( j = 0; j < cgrid__size.ny; j++ )
	  {
	    for ( k = 0; k < cgrid__size.nz; k++ )
	      {
		cgrid__coords.push_back(cgrid__x[i]);
		cgrid__coords.push_back(cgrid__y[j]);
		cgrid__coords.push_back(cgrid__z[k]);
	      }
	  }

      }
  }

void Cgrid::mapping__localize( const YYZVis::YinYangVolumeObject& yin_volume, const YYZVis::YinYangVolumeObject& yang_volume, const YYZVis::ZhongVolumeObject& zhong_volume )
  {
    int i, j, k;
    float x, y, z;
    float r, t, p;
    int index = 0;

    for( k = 0; k < cgrid__size.nz; k++ )
      {
        x = cgrid__x[k];
	for( j = 0; j < cgrid__size.ny; j++ )
	  {
	    y = cgrid__y[j];
	    for( i = 0; i < cgrid__size.nx; i++ )
	      {
		z = cgrid__x[i];
		
		r = sqrt( x*x + y*y + z*z );
		t = acos( z/r );
		p = atan2( y, x );

		if ( r > yin_volume.rangeR().max ) continue;
		if ( r <= yin_volume.rangeR().min )
		  {
		     // ZHONG
		     this->iFind_zhong( x, y, z, index, zhong_volume );
		     continue;
		  }
		if ( t >= yin_volume.rangeTheta().min && t <= yin_volume.rangeTheta().max &&
		     p >= yin_volume.rangePhi().min && p <= yin_volume.rangePhi().max )
	            {
	              // YIN
		      this->ogrid__find_near_corner( r, t, p, index, yin_volume );
	            }
		else
		    {
		      // YANG
		      float x_ = -x;
		      float y_ = z;
		      float z_ = y;
		      t = acos( z_/r );
		      p = atan2( y_, x_ );
		      this->ogrid__find_near_corner( r, t, p, index, yang_volume );
		    }
		index++;
	      }	 
	  } 
      }
  }
  
  void Cgrid::iFind_zhong(float x, float y, float z, int index, const YYZVis::ZhongVolumeObject& z_object )
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

    this->igrid_to_cgrid_localize(i1, j1, k1, wx1, wy1, wz1, index, z_object);      
    return;
   }

  int Cgrid::igrid__find_nearleft( char axis, float c, const YYZVis::ZhongVolumeObject& object )
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
    return 0; // ???
  }

  void Cgrid::ogrid__find_near_corner( float rad, float theta, float phi, int index, const YYZVis::YinYangVolumeObject& object  )
  {                            
    int i1, j1, k1=0;                         
    float wr1, wt1, wp1=0;          
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
   
    this-> ogrid_to_cgrid_localize( i1, j1, k1, wr1, wt1, wp1, rad, theta, phi,index, object );
   
      }
  void Cgrid::ogrid_to_cgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index, const YYZVis::YinYangVolumeObject& object )
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
 
     cgrid__values[index] = value;
   }

  void Cgrid::igrid_to_cgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index, const YYZVis::ZhongVolumeObject& object)
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

     cgrid__values[index] = value;
   
   }
}  // end of namespace local
