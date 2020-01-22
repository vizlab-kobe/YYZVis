#include "cgrid.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <vector>
#include <math.h>

namespace local
{
  Cgrid::Cgrid( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume )
  {
    this->cgrid__make(); 
    this->mapping__localize( yin_volume, yang_volume, zhong_volume );
  }

  void Cgrid::cgrid__make()
  {
    this->set_minmax();
    this->set_xyz();
    this->set_metric();
  }  
  
  void Cgrid::set_minmax()
  {
    cgrid__x_min = -1.0;
    cgrid__x_max =  1.0;
    cgrid__y_min = -1.0;
    cgrid__y_max =  1.0;
    cgrid__z_min = -1.0;
    cgrid__z_max =  1.0;
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

void Cgrid::mapping__localize( const YinYangVis::YinYangVolumeObject& yin_volume, const YinYangVis::YinYangVolumeObject& yang_volume, const YinYangVis::ZhongVolumeObject& zhong_volume )
  {
    int i, j, k;
    float x, y, z;
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
		if ( sqrt( x*x + y*y + z*z ) <= yin_volume.rangeR().min )
		  {
		this->iFind_zhong( x, y, z, index, zhong_volume );
		continue;
		  }
		float p = atan2( y, x );
		if (  p >= yin_volume.rangePhi().min && p <= yin_volume.rangePhi().max )
	            {
	              // YIN
		      this->iFind( x, y, z, index, yin_volume );
	            }
		else
		    {
		      // YANG
		      this->iFind( x, y, z, index, yang_volume );
		    }
		index++;
	      }	 
	  } 
      }
  }

  void Cgrid::iFind(float x, float y, float z, int index, const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    float polar[3];    //{ r, t, p }
   
    if ( sqrt ( x*x + y*y + z*z ) >= yoy_object.rangeR().max )
      {
	return;
      }
    else if ( sqrt ( x*x + y*y + z*z ) <= yoy_object.rangeR().min )
      {
	return;
      }

    this->cgrid__xyz2rtp( x, y, x, polar, yoy_object );

    this->ogrid__find_near_corner( polar[0], polar[1], polar[2], index, yoy_object);
      }
  
  void Cgrid::iFind_zhong(float x, float y, float z, int index, const YinYangVis::ZhongVolumeObject& z_object )
  {
    int i1,j1,k1;
    float wx1,wy1,wz1;
    
    i1 =  igrid__find_nearleft('x', x, z_object);
    j1 =  igrid__find_nearleft('y', y, z_object);
    k1 =  igrid__find_nearleft('z', z, z_object);

    wx1 = ( -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(i1+1) - x ) / (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1);
    wy1 = (  -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(j1+1) - y ) / (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1);
    wz1 = (  -(z_object.rangeR().min+z_object.rangeR().d*2) + (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1)*(k1+1) - z ) / (z_object.rangeR().min+z_object.rangeR().d*2)*2/(z_object.dim()-1);
	
    /*	wx1 = ( igrid__x[i1 + 1] - cart[0] ) / igrid__dx;
	wy1 = ( igrid__y[j1 + 1] - cart[1] ) / igrid__dy;
	wz1 = ( igrid__z[k1 + 1] - cart[2] ) / igrid__dz;*/

    this->igrid_to_cgrid_localize(i1, j1, k1, wx1, wy1, wz1, index, z_object);      
    return;
   }

  void Cgrid::cgrid__xyz2rtp( float x, float y, float z, float polar[3], const YinYangVis::YinYangVolumeObject& object )
  {
    float r, t, p;

    r = sqrt(  x*x + y*y + z*z );
 
    t = acos( z/r );

    if ( t >= object.rangeTheta().min && t <= object.rangeTheta().max )
      {
	p = atan2( y, x );
      }
    if( object.gridType() == YinYangVis::YinYangVolumeObject::Yang)
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

  int Cgrid::igrid__find_nearleft( char axis, float c, const YinYangVis::ZhongVolumeObject& object )
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
  }

  void Cgrid::ogrid__find_near_corner( float rad, float theta, float phi, int index, const YinYangVis::YinYangVolumeObject& object  )
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
  void Cgrid::ogrid_to_cgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, float rad, float tht, float phi, int index, const YinYangVis::YinYangVolumeObject& object )
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

  void Cgrid::igrid_to_cgrid_localize(int i1, int j1, int k1, float wr1, float wt1, float wp1, int index, const YinYangVis::ZhongVolumeObject& object)
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
