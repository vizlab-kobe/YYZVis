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
    this->ogrid__make( yin_volume );
    this->igrid__make( zhong_volume ); 
    this->cgrid__make(); 
    this->mapping__localize();
    this->ogrid__make( yang_volume );
    this->mapping__localize();
 
  }

  void Cgrid::ogrid__make( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    this->set_o_minmax( yoy_object );
    this->set_o_nrtp( yoy_object );
    this->set_o_metric( yoy_object );
    ogrid__values = yoy_object.values();
  }

  void Cgrid::igrid__make( const YinYangVis::ZhongVolumeObject& z_object  )
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

  void Cgrid::set_o_minmax( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    OGRID__THETA_FROM = yoy_object.rangeTheta().min;
    OGRID__THETA_TO = yoy_object.rangeTheta().max;
    OGRID__PHI_FROM = yoy_object.rangePhi().min;
    OGRID__PHI_TO = yoy_object.rangePhi().max;
       
    ogrid__range_r = yoy_object.rangeR();
    ogrid__range_theta = yoy_object.rangeTheta();
    ogrid__range_phi = yoy_object.rangePhi();
  }

  void Cgrid::set_o_nrtp( const YinYangVis::YinYangVolumeObject& yoy_object )
  {
    ogrid__size.nr = yoy_object.dimR();
    ogrid__size.nt = yoy_object.dimTheta();
    ogrid__size.np = yoy_object.dimPhi();
    
    ogrid__rad.allocate(ogrid__size.nr);
    ogrid__theta.allocate(ogrid__size.nt);
    ogrid__phi.allocate(ogrid__size.np);
  }

  void Cgrid::set_o_metric( const YinYangVis::YinYangVolumeObject& yoy_object )
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
  
  void Cgrid::set_minmax()
  {
    CGRID__X_FROM = -1.0;
    CGRID__X_TO   = 1.0;
    CGRID__Y_FROM   = -1.0;
    CGRID__Y_TO     =  1.0;
    CGRID__Z_FROM   = -1.0;
    CGRID__Z_TO     =  1.0;
  }
    
  void Cgrid::set_xyz()
  {
    float dx,dy,dz;
    float culling;
        
    culling = 2.5;

    cgrid__size.nx = 350;
    cgrid__size.ny = 350;
    cgrid__size.nz = 350;
        
    cgrid__dx = culling * (CGRID__X_TO - CGROD__X_TO) / cgrid__size_nx;
    cgrid__dy = culling * (CGRID__Y_TO - CGROD__Y_TO) / cgrid__size_ny;
    cgrid__dz = culling * (CGRID__Z_TO - CGROD__Z_TO) / cgrid__size_nz;

    cgrid__size.nx = ( CGRID__X_TO - CGRID__X_FROM ) / dx + 1;
    cgrid__size.ny = ( CGRID__Y_TO - CGRID__Y_FROM ) / dy + 1;
    cgrid__size.nz = ( CGRID__Z_TO - CGRID__Z_FROM ) / dz + 1;
    
    cgrid__x.allocate(cgrid__size.nx);
    cgrid__y.allocate(cgrid__size.ny);
    cgrid__z.allocate(cgrid__size.nz);
   
    cgrid__values.reserve(cgrid__size.nx*cgrid__size.ny*cgrid__size.nz);
    cgrid__coords.reserve(3*scgrid__size.nx*cgrid__size.ny*cgrid__size.nz);
    
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
    
  void Cgrid::cgrid__make()
  {
    this->set_minmax();
    this->set_rtp();
    this->set_metric();
  }
  

  void Sgrid::mapping__localize()
  {
    int i, j, k;
    float x, y, z;
    int index = 0;
    
    for( k = 0; k < cgrid__size.nx; k++ )
      {
        z = cgrid__z[k];
	for( j = 0; j < sgrid__size.nt; j++ )
	  {
	    y = cgrid__y[j];
	    for( i = 0; i < sgrid__size.nr; i++ )
	      {
		x = cgrid__x[i];
	        this->iFind( x, y, z, index );
		index++;
	      }
	  }
      }
  }

  void Sgrid::iFind(float x, float y, float z, int index )
  {
    //   float cart[3];     //{ x, y, z }
    int i1, j1, k1;

    //  float polar[3];    //{ r, t, p }
    float wx1, wy1, wz1;
    
   
    if ( x*x + y*y + z*z >= ogrid__range_r.max )
      {
	return;
      }
    else if ( x*x + y*y + z*z <= ogrid__range_r.min )
      {
	// Zhong
	i1 =  igrid__find_nearleft('x', x);
	j1 =  igrid__find_nearleft('y', y);
	k1 =  igrid__find_nearleft('z', z);

	wx1 = ( igrid__x[i1 + 1] - x ) / igrid__dx;
	wy1 = ( igrid__y[j1 + 1] - y ) / igrid__dy;
	wz1 = ( igrid__z[k1 + 1] - z ) / igrid__dz;

	this->igrid_to_sgrid_localize(i1, j1, k1, wx1, wy1, wz1, index);      
	return;
      }
    
    this->ogrid__find_near_corner( polar[0], polar[1], polar[2], index);
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

  void Sgrid::ogrid__find_near_corner( float x,float y,float z,int index  )
  {                            
    int i1, j1, k1;                         
    float wx1, wy1, wz1;          
    int i, j, k; 

    for ( i = 0; i < ogrid__size.nr-2; i++ )
      {
	if ( x >= ogrid__xrad[i] && x <= ogrid__rad[i+1] ) 
	  {
	    i1 = i;
	    wx1 = (ogrid__rad[i+1]-rad)/ogrid__range_r.d;
	    break;
	  }
      }

    for ( j = 0; j < ogrid__size.nt-2; j++ )
      {
	if ( y >= ogrid__theta[j] && y <= ogrid__theta[j+1] ) 
	  {
	    j1 = j;
	    wy1 = (ogrid__theta[j+1]-y)/ogrid__range_theta.d;
	    break;
	  } 
      }

    for ( k = 0; k < ogrid__size.np-2; k++ )
      {
	if ( phi >= ogrid__phi[k] && phi <= ogrid__phi[k+1] ) 
	  {
	    k1 = k;
	    wz1 = (ogrid__phi[k+1]-phi)/ogrid__range_phi.d;
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
