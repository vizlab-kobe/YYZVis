#include "ExternalFaces.h"


namespace YinYangVis
{

ExternalFaces::ExternalFaces():
    kvs::MapperBase(),
    kvs::PolygonObject()
{
}

ExternalFaces::ExternalFaces( const kvs::VolumeObjectBase* volume ):
    kvs::MapperBase(),
    kvs::PolygonObject()
{
    this->exec( volume );
}

ExternalFaces::ExternalFaces( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->exec( volume );
}

ExternalFaces::SuperClass* ExternalFaces::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is NULL.");
        return NULL;
    }

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is not volume dat.");
        return NULL;
    }

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );

    if ( YinYangVis::ZhongVolumeObject::DownCast( volume ) )
    {
        const YinYangVis::ZhongVolumeObject* zvolume = YinYangVis::ZhongVolumeObject::DownCast( volume );
        this->mapping( zvolume );
    }
    else if ( YinYangVis::YinYangVolumeObject::DownCast( volume ) )
    {
        const YinYangVis::YinYangVolumeObject* yvolume = YinYangVis::YinYangVolumeObject::DownCast( volume );
        this->mapping( yvolume );
    }
    
    return this;
}

void ExternalFaces::mapping( const YinYangVis::ZhongVolumeObject* zvolume )
{
    // Calculate coords, colors, and normals.
     this->calculate_zhong_coords( zvolume );
     this->calculate_zhong_colors( zvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

void ExternalFaces::mapping( const YinYangVis::YinYangVolumeObject* yvolume )
{

    // Calculate coords, colors, and normals.
      this->calculate_yinyang_coords( yvolume );
      this->calculate_yinyang_colors( yvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

void ExternalFaces::calculate_zhong_coords( const YinYangVis::ZhongVolumeObject* zvolume )  
{
  const size_t dim = zvolume->dim();
  float x_min,x_max,
        y_min,y_max,
        z_min,z_max;
  float dx,dy,dz;

  kvs::ValueArray<float> xcoords, ycoords, zcoords; 

  const size_t nfaces = ( dim - 1 ) * ( dim - 1 ) * 6 * 2;
  const size_t nverts = nfaces * 3;
  kvs::ValueArray<kvs::Real32> coords( 3 * nverts);
  kvs::Real32* coord = coords.data();
  kvs::ValueArray<kvs::Real32> normals( 3 * nfaces );
  kvs::Real32* normal = normals.data();

  x_min = -0.35;
  x_max = 0.35;
  y_min = -0.35;
  y_max = 0.35;
  z_min = -0.35;
  z_max = 0.35;
  std ::cout << zvolume->rangeR().max;
  dx = ( x_max*2 ) / ( dim - 1 );
  dy = ( y_max*2 ) / ( dim - 1 );
  dz = ( z_max*2 ) / ( dim - 1 );

  xcoords.allocate( dim );
  ycoords.allocate( dim );
  zcoords.allocate( dim );

  for ( size_t i = 0; i < dim; i++ )
    {
      xcoords[i] = x_min + dx * i;
      ycoords[i] = y_min + dy * i;
      zcoords[i] = z_min + dz * i;
    }

  // XY (Z=-0.35) plane.
  {
    const float z = z_min;
    const kvs::Vector3f n( 0.0f, 0.0f, -1.0f );
    for ( size_t j = 0; j < dim - 1; j++ )
      {
	const float y = y_min + dy * j;
	for ( size_t i = 0; i < dim - 1; i++ )
	  {
	    const float x = x_min + dx * i;
	    // v3
	    *( coord++ ) = x;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z;
	    // v2
	    *( coord++ ) = x + dx;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z;
	    // v1
	    *( coord++ ) = x + dx;
	    *( coord++ ) = y;
	    *( coord++ ) = z;

	    // v1
	    *( coord++ ) = x + dx;
	    *( coord++ ) = y;
	    *( coord++ ) = z;
	    // v0
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z;
	    // v3
	    *( coord++ ) = x;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z;

	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	  }
      }
  }

  // XY (Z=0.35) plane.
  {
    const float z = z_max;
    const kvs::Vector3f n( 0.0f, 0.0f, 1.0f );
    for ( size_t j = 0; j < dim - 1; j++ )
      {
	const float y = y_min + dy * j;
	for ( size_t i = 0; i < dim - 1; i++ )
	  {
	    const float x = x_min + dx * i;
	// v0
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z;
	// v1
	*( coord++ ) = x + dx;
	*( coord++ ) = y;
	*( coord++ ) = z;
	// v2
	*( coord++ ) = x + dx;
	*( coord++ ) = y + dy;
	*( coord++ ) = z;

	// v2
	*( coord++ ) = x + dx;
	*( coord++ ) = y + dy;
	*( coord++ ) = z;
	// v3
	*( coord++ ) = x;
	*( coord++ ) = y + dy;
	*( coord++ ) = z;
	// v0
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z;

	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
      }
  }
}

  // YZ (X=-0.35) plane.
  {
    const float x = x_min;
    const kvs::Vector3f n( -1.0f, 0.0f, 0.0f );
    for ( size_t j = 0; j < dim - 1; j++ )
      {
	const float y = y_min + dy * j;
	for ( size_t k = 0; k < dim - 1; k++ )
	  {
	    const float z = z_min + dz * k;
	    // v0
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z;
	    // v1
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z + dz;
	    // v2
	    *( coord++ ) = x;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z + dz;

	    // v2
	    *( coord++ ) = x;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z + dz;
	    // v3
	    *( coord++ ) = x;
	    *( coord++ ) = y + dy;
	    *( coord++ ) = z;
	    // v0
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z;

	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	  }
      }
  }

  // YZ (X=0.35) plane.
  {
    const float x = x_max;
    const kvs::Vector3f n( 1.0f, 0.0f, 0.0f );
    for ( size_t j = 0; j < dim - 1; j++ )
      {
	const float y = y_min + dy * j;
    for ( size_t k = 0; k < dim - 1; k++ )
      {
	const float z = z_min + dz * k;
	// v3
	*( coord++ ) = x;
	*( coord++ ) = y + dy;
	*( coord++ ) = z;
	// v2
	*( coord++ ) = x;
	*( coord++ ) = y + dy;
	*( coord++ ) = z + dz;
	// v1
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z + dz;

	// v1
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z + dz;
	// v0
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z;
	// v3
	*( coord++ ) = x;
	*( coord++ ) = y + dy;
	*( coord++ ) = z;

	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
      }
  }
}

// XZ (Y=-0.35) plane.
  {
    const float y = y_min;
    const kvs::Vector3f n( 0.0f, -1.0f, 0.0f );
    for ( size_t k = 0; k < dim - 1; k++ )
      {
	const float z = z_min + dz * k;
	for ( size_t i = 0; i < dim - 1; i++ )
      {
	const float x = x_min + dx * i;
	// v0
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z;
	// v1
	*( coord++ ) = x + dz;
	*( coord++ ) = y;
	*( coord++ ) = z;
	// v2
	*( coord++ ) = x + dx;
	*( coord++ ) = y;
	*( coord++ ) = z + dz;

	// v2
	*( coord++ ) = x + dx;
	*( coord++ ) = y;
	*( coord++ ) = z + dz;
	// v3
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z + dz;
	// v0
	*( coord++ ) = x;
	*( coord++ ) = y;
	*( coord++ ) = z;

	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
	// n0
	*( normal++ ) = n.x();
	*( normal++ ) = n.y();
	*( normal++ ) = n.z();
      }
  }
}

// XZ (Y=0.35) plane.
  {
    const float y = y_max;
      const kvs::Vector3f n( 0.0f, 1.0f, 0.0f );
    for ( size_t k = 0; k < dim - 1; k++ )
      {
	const float z = z_min + dz * k;
	for ( size_t i = 0; i < dim - 1; i++ )
	  {
	    const float x = x_min + dx * i;
	    // v3
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z + dz;
	      // v2
            *( coord++ ) = x + dx;
	    *( coord++ ) = y;
	    *( coord++ ) = z + dz;
	    // v1
	    *( coord++ ) = x + dx;
	    *( coord++ ) = y;
	    *( coord++ ) = z;

	    // v1
	    *( coord++ ) = x + dx;
	    *( coord++ ) = y;
	    *( coord++ ) = z;
	    // v0
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z;
	    // v3
	    *( coord++ ) = x;
	    *( coord++ ) = y;
	    *( coord++ ) = z + dz;

	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	    // n0
	    *( normal++ ) = n.x();
	    *( normal++ ) = n.y();
	    *( normal++ ) = n.z();
	  }
      }
  }

		    	
SuperClass::setCoords( coords );
SuperClass::setNormals( normals );    
  
}

void ExternalFaces::calculate_yinyang_coords( const YinYangVis::YinYangVolumeObject* yvolume )
{
     size_t dim_r = 0, dim_theta = 0, dim_phi = 0;

     float theta_from,theta_to,phi_from,phi_to,r_from,r_to;
     float r_d,theta_d,phi_d;
     size_t index = 0;

     r_from = yvolume->rangeR().min;
     r_to =  yvolume->rangeR().max;
     theta_from =  yvolume->rangeTheta().min;
     theta_to =  yvolume->rangeTheta().max;
     phi_from =  yvolume->rangePhi().min;
     phi_to =  yvolume->rangePhi().max;

     r_d = yvolume->rangeR().d;
     theta_d = yvolume->rangeTheta().d;
     phi_d = yvolume->rangePhi().d;

     dim_r = yvolume->dimR();
     dim_theta= yvolume->dimTheta();
     dim_phi = yvolume->dimPhi();

    const size_t nfaces = ((dim_r - 1) * (dim_theta - 1) * 2
                        + (dim_r - 1) * (dim_phi - 1) * 2
			+ (dim_theta - 1) * (dim_phi - 1) * 2) * 2;
    const size_t nverts = nfaces * 3;

    kvs::ValueArray<kvs::Real32> coords( 3 * nverts );
    kvs::Real32* coord = coords.data();
    kvs::ValueArray<kvs::Real32> normals( 3 * nfaces );
    kvs::Real32* normal = normals.data();
    
    // phi = 0	
    {
      const float phi = phi_from;
      const float sin_phi = std::sin( phi );
      const float cos_phi = std::cos( phi );
      for ( size_t j = 0; j < dim_theta - 1; j++ )
	{
	  const float theta = theta_from + theta_d * j;
	  const float theta_next = theta + theta_d;
	  const float sin_theta = std::sin( theta );
	  const float cos_theta = std::cos( theta );
	  const float sin_theta_next = std::sin( theta_next );
	  const float cos_theta_next = std::cos( theta_next );
	  for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
	      const float r = r_from + r_d * i;
	      const float r_next = r + r_d;
	      
	      // v3
	      const float x3 = r * sin_theta_next * cos_phi;
	      const float y3 = r * sin_theta_next * sin_phi;
	      const float z3 = r * cos_theta_next;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
	      // v2
	      const float x2 = r_next * sin_theta_next * cos_phi;
	      const float y2 = r_next * sin_theta_next * sin_phi;
	      const float z2 = r_next * cos_theta_next;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
	      // v1
	      const float x1 = r_next * sin_theta * cos_phi;
	      const float y1 = r_next * sin_theta * sin_phi;
	      const float z1 = r_next * cos_theta;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;

	      // v1
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
	      // v0
	      const float x0 = r * sin_theta * cos_phi;
	      const float y0 = r * sin_theta * sin_phi;
	      const float z0 = r * cos_theta;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
	      // v3
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;

	      if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
		{
	          this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
	          index++;
	          this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
		  index++;
		}
	      else
		{
		  this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
	          index++;
	          this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
		  index++;
		}
	    }
	}

    }
    
    // phi = dim_phi - 1
    {
      const float phi = phi_to;  //grid_size.z() * ( dim_phi - 1 );
      const float sin_phi = std::sin( phi );
      const float cos_phi = std::cos( phi );
      for ( size_t j = 0; j < dim_theta - 1; j++ )
	{
	  const float theta = theta_from + theta_d * j;
	  const float theta_next = theta + theta_d;
	  const float sin_theta = std::sin( theta );
	  const float cos_theta = std::cos( theta );
	  const float sin_theta_next = std::sin( theta_next );
	  const float cos_theta_next = std::cos( theta_next );
	  for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
	      const float r = r_from + r_d * i;
	      const float r_next = r + r_d;
	      
	      // v0
	      const float x0 = r * sin_theta * cos_phi;
	      const float y0 = r * sin_theta * sin_phi;
	      const float z0 = r * cos_theta;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
	      // v1
	      const float x1 = r_next * sin_theta * cos_phi;
	      const float y1 = r_next * sin_theta * sin_phi;
	      const float z1 = r_next * cos_theta;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
	      // v2
	      const float x2 = r_next * sin_theta_next * cos_phi;
	      const float y2 = r_next * sin_theta_next * sin_phi;
	      const float z2 = r_next * cos_theta_next;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;

	      // v2
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
	      // v3
	      const float x3 = r * sin_theta_next * cos_phi;
	      const float y3 = r * sin_theta_next * sin_phi;
	      const float z3 = r * cos_theta_next;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
	      // v0
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	      *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;

	       if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
		 {
		   this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
		   index++;
		   this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
		   index++;
		 }
	       else
		 {
		   this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
		   index++;
		   this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
		   index++;
		 }
	    }
	}
     
    // r = 0
    const float r = r_from;
    for  ( size_t j = 0; j < dim_theta - 1; j++ ) 
    {
         const float theta = theta_from + theta_d * j;
	 const float theta_next = theta + theta_d;
	 const float sin_theta = std::sin( theta );
	 const float cos_theta = std::cos( theta );
	 const float sin_theta_next = std::sin( theta_next );
	 const float cos_theta_next = std::cos( theta_next );
         for ( size_t k = 0; k < dim_phi - 1; k++ )
	 {
	   const float phi = phi_from + phi_d * k;
	   const float phi_next = phi + phi_d;
		    
	   const float sin_phi = std::sin( phi );
	   const float cos_phi = std::cos( phi );
	   const float sin_phi_next = std::sin( phi_next );
	   const float cos_phi_next = std::cos( phi_next );
	   
           // v0
	   const float x0 = r * sin_theta * cos_phi;
	   const float y0 = r * sin_theta * sin_phi;
	   const float z0 = r * cos_theta;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
	   // v1
	   const float x1 = r * sin_theta * cos_phi_next;
	   const float y1 = r * sin_theta * sin_phi_next;
	   const float z1 = r * cos_theta;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
	   // v2
	   const float x2 = r * sin_theta_next * cos_phi_next;
	   const float y2 = r * sin_theta_next * sin_phi_next;
	   const float z2 = r * cos_theta_next;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;

	   // v2
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
	   // v3
	   const float x3 = r * sin_theta_next * cos_phi;
	   const float y3 = r * sin_theta_next * sin_phi;
	   const float z3 = r * cos_theta_next;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
	   // v0
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	   *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;

	   if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
	     {
	       this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
	       index++;
	       this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
	       index++;
	     }
	   else
	     {
	       this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
	       index++;
	       this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
	       index++;
	     }
	 }
    }
    }

    
    // r = dim_r - 1
    {
         const float r = r_to; //grid_size.x() * ( dim_r - 1 );
	 for ( size_t j = 0; j < dim_theta - 1; j++ )
	 {
	   const float theta = theta_from + theta_d * j;
	   const float theta_next = theta + theta_d;
	   const float sin_theta = std::sin( theta );
	   const float cos_theta = std::cos( theta );
	   const float sin_theta_next = std::sin( theta_next );
	   const float cos_theta_next = std::cos( theta_next );
	   for ( size_t k = 0; k < dim_phi - 1; k++ )
	     {
	       const float phi = phi_from + phi_d * k;
	       const float phi_next = phi + phi_d;
		    
	       const float sin_phi = std::sin( phi );
	       const float cos_phi = std::cos( phi );
	       const float sin_phi_next = std::sin( phi_next );
	       const float cos_phi_next = std::cos( phi_next );
	     
	        // v3
	       const float x3 = r * sin_theta_next * cos_phi;
	       const float y3 = r * sin_theta_next * sin_phi;
	       const float z3 = r * cos_theta_next;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
	       // v2
	       const float x2 = r * sin_theta_next * cos_phi_next;
	       const float y2 = r * sin_theta_next * sin_phi_next;
	       const float z2 = r * cos_theta_next;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
	       // v1
	       const float x1 = r * sin_theta * cos_phi_next;
	       const float y1 = r * sin_theta * sin_phi_next;
	       const float z1 = r * cos_theta;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;

	       // v1
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
	       // v0
	       const float x0 = r * sin_theta * cos_phi;
	       const float y0 = r * sin_theta * sin_phi;
	       const float z0 = r * cos_theta;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
	       // v3
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
	       *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;

	       if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
		 {
		   this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
		   index++;
		   this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
		   index++;
		 }
	       else
		 {
		   this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
		   index++;
		   this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
		   index++;
		 }
	     }
	 }
    }

    
    // theta = 0
    {
        const float theta = theta_from;
	const float sin_theta = std::sin( theta );
	const float cos_theta = std::cos( theta );
        for ( size_t k = 0; k < dim_phi - 1; k++ )
	{
	    const float phi = phi_from + phi_d * k;
	    const float phi_next = phi + phi_d;
	    for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
	        const float r = r_from + r_d * i;
		const float r_next = r + r_d;
		    
		const float sin_phi = std::sin( phi );
		const float cos_phi = std::cos( phi );
		const float sin_phi_next = std::sin( phi_next );
		const float cos_phi_next = std::cos( phi_next );
		
		// v0
		const float x0 = r * sin_theta * cos_phi;
		const float y0 = r * sin_theta * sin_phi;
		const float z0 = r * cos_theta;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
		// v1
		const float x1 = r_next * sin_theta * cos_phi;
		const float y1 = r_next * sin_theta * sin_phi;
		const float z1 = r_next * cos_theta;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
		// v2
		const float x2 = r_next * sin_theta * cos_phi_next;
		const float y2 = r_next * sin_theta * sin_phi_next;
		const float z2 = r_next * cos_theta;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;

		// v2
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
		// v3
		const float x3 = r * sin_theta * cos_phi_next;
		const float y3 = r * sin_theta * sin_phi_next;
		const float z3 = r * cos_theta;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
		// v0
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
		*( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;

		if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
		  {
		    this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
		    index++;
		    this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
		    index++;
		  }
		else
		  {
		    this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
		    index++;
		    this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
		    index++;
		  }
	     }
	}
    }

	// theta = dim_theta - 1
	{
	  const float theta = theta_to; //grid_size.y() * ( dim_theta - 1 );
	  const float sin_theta = std::sin( theta );
	  const float cos_theta = std::cos( theta );
	    for ( size_t k = 0; k < dim_phi - 1; k++ )
	    {
	        const float phi = phi_from + phi_d * k;
		const float phi_next = phi + phi_d;
		for ( size_t i = 0; i < dim_r - 1; i++ )
		{
		    const float r = r_from + r_d * i;
		    const float r_next = r + r_d;
		    
		    const float sin_phi = std::sin( phi );
		    const float cos_phi = std::cos( phi );
		    const float sin_phi_next = std::sin( phi_next );
		    const float cos_phi_next = std::cos( phi_next );
		    
		    // v3
		    const float x3 = r * sin_theta * cos_phi_next;
                    const float y3 = r * sin_theta * sin_phi_next;
                    const float z3 = r * cos_theta;
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;
		    // v2
		    const float x2 = r_next * sin_theta * cos_phi_next;
                    const float y2 = r_next * sin_theta * sin_phi_next;
                    const float z2 = r_next * cos_theta;
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z2 : y2;
		    // v1
		    const float x1 = r_next * sin_theta * cos_phi;
                    const float y1 = r_next * sin_theta * sin_phi;
                    const float z1 = r_next * cos_theta;
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;

		    // v1
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z1 : y1;
		    // v0
		    const float x0 = r * sin_theta * cos_phi;
                    const float y0 = r * sin_theta * sin_phi;
                    const float z0 = r * cos_theta;
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z0 : y0;
		    // v3
		    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                    *( coord++) = ( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z3 : y3;

		    if( yvolume->gridType() == YinYangVis::YinYangVolumeObject::Yin )
		      {
			this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
			index++;
			this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
			index++;
		      }
		    else
		      {
			this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
			index++;
			this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
			index++;
		      }
		}
	    }
	}

	SuperClass::setCoords( coords );
      	SuperClass::setNormals( normals );
}

void ExternalFaces::calculate_yinyang_colors( const YinYangVis::YinYangVolumeObject* yvolume
				    )
{
   size_t dim_r = 0, dim_theta = 0, dim_phi = 0;
   const kvs::Real64 min_value = yvolume->minValue();
   const kvs::Real64 max_value = yvolume->maxValue();
   kvs::AnyValueArray value =  yvolume->values();

   kvs::UInt32 node_index[4];
   kvs::UInt32 color_level[4];
   const kvs::ColorMap cmap( BaseClass::colorMap() );

   dim_r = yvolume->dimR();
   dim_theta= yvolume->dimTheta();
   dim_phi = yvolume->dimPhi();

   float theta_from,theta_to,phi_from,phi_to,r_from,r_to;
   float r_d,theta_d,phi_d;

   r_from = yvolume->rangeR().min;
   r_to =  yvolume->rangeR().max;
   theta_from =  yvolume->rangeTheta().min;
   theta_to =  yvolume->rangeTheta().max;
   phi_from =  yvolume->rangePhi().min;
   phi_to =  yvolume->rangePhi().max;

   r_d = yvolume->rangeR().d;
   theta_d = yvolume->rangeTheta().d;
   phi_d = yvolume->rangePhi().d;

   const size_t nfaces = ((dim_r - 1) * (dim_theta - 1) * 2
			  + (dim_r - 1) * (dim_phi - 1) * 2
			  + (dim_theta - 1) * (dim_phi - 1) * 2) * 2;
   const size_t nverts = nfaces * 3;

   kvs::ValueArray<kvs::UInt8> colors( 3 * nverts );
   kvs::UInt8* color = colors.data();

   // phi = 0 
    {
        const size_t k = 0;
	const size_t offset0 = k * dim_r;
        for ( size_t j = 0, offset = offset0 ; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset++ )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(),node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }

    // phi = dim_phi - 1
    {
        const size_t k = dim_phi - 1;
        const size_t offset0 = k * dim_r * dim_theta;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // r = 0
    {
        const size_t i = 0;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1 ; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t k = 0; k < dim_phi - 1; k++, offset += dim_r * dim_theta )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + ( dim_r * dim_theta );
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // r = dim_r - 1
    {
        const size_t i = dim_r - 1;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
	  for ( size_t k = 0; k < dim_phi - 1; k++, offset += ( dim_r * dim_theta ) )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + ( dim_r * dim_theta );
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }

    // theta = 0
    {
        const size_t j = 0;
        const size_t offset0 = j * dim_r;
        for ( size_t k = 0, offset = offset0; k < dim_phi - 1; k++, offset =  offset0 + k * ( dim_r * dim_theta ) )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + ( dim_r * dim_theta );
                node_index[3] = node_index[0] + ( dim_r * dim_theta );
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // theta = dim_theta - 1
    {
        const size_t j = dim_theta;
        const size_t offset0 = j * dim_r;
        for ( size_t k = 0, offset = offset0; k < dim_phi - 1; k++, offset = offset0 + k * ( dim_r * dim_theta ) )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + ( dim_r * dim_theta );
                node_index[3] = node_index[0] + ( dim_r * dim_theta );
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }
    SuperClass::setColors( colors );
    }

  void ExternalFaces::calculate_zhong_colors( const YinYangVis::ZhongVolumeObject* zvolume )
  {
    // Parameters of the volume data.
    const size_t dim = zvolume->dim();
    const kvs::Real64 min_value = zvolume->minValue();
    const kvs::Real64 max_value = zvolume->maxValue();
    kvs::AnyValueArray value =  zvolume->values();
    const size_t nnodes_per_line = dim;
    const size_t nnodes_per_slice = dim * dim;

    kvs::UInt32 node_index[4];
    kvs::UInt32 color_level[4];
 
    const size_t nfaces = (dim - 1) * (dim - 1) * 6 * 2;
    const size_t nverts = nfaces * 3;
    kvs::ValueArray<kvs::UInt8> colors( 3 * nverts );
    kvs::UInt8* color = colors.data();	       
    const kvs::ColorMap cmap( BaseClass::colorMap() );

      // XY (Z=-0.35) plane.
      {
        const size_t k = 0;
        const size_t offset0 = k * nnodes_per_line;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
	  {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
	      }
	  }
      }

      // XY (Z=0.35) plane.
      {
        const size_t k = dim - 1;
        const size_t offset0 = k * nnodes_per_slice;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
	  {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
	      }
	  }
      }

      // YZ (X=-0.35) plane.
      {
        const size_t i = 0;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
	  {
            for ( size_t k = 0; k < dim - 1; k++, offset += nnodes_per_slice )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + nnodes_per_slice;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
	      }
	  }
      }

      // YZ (X=0.35) plane.
      {
        const size_t i = dim - 1;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
	  {
            for ( size_t k = 0; k < dim - 1; k++, offset += nnodes_per_slice )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + nnodes_per_slice;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
	      }
	  }
      }

      // XZ (Y=-0.35) plane.
      {
        const size_t j = 0;
        const size_t offset0 = j * nnodes_per_line;
        for ( size_t k = 0, offset = offset0; k < dim - 1; k++, offset =  offset0 + k * nnodes_per_slice )
	  {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_slice;
                node_index[3] = node_index[0] + nnodes_per_slice;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
	      }
	  }
      }

      // XZ (Y=0.35) plane.
      {
        const size_t j = dim - 1;
        const size_t offset0 = j * nnodes_per_line;
        for ( size_t k = 0, offset = offset0; k < dim - 1; k++, offset = offset0 + k * nnodes_per_slice )
	  {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
	      {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_slice;
                node_index[3] = node_index[0] + nnodes_per_slice;
                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
	      }
	  }
      }

    SuperClass::setColors( colors );
}

void ExternalFaces::calculate_normal( const float x0, const float y0, const float z0,
					const float x1, const float y1, const float z1,
				      const float x2, const float y2, const float z2, kvs::Real32* normal, size_t index )
{
    float x4=0,y4=0,z4=0;
    float x5=0,y5=0,z5=0;
    float x6=0,y6=0,z6=0;
    float er=0;

    x4 = x1 - x0;
    y4 = y1 - y0;
    z4 = z1 - z0;

    x5 = x2 - x0;
    y5 = y2 - y0;
    z5 = z2 - z0;

    x6 = y4 * z5 - y5 * z4;
    y6 = z4 * x5 - x4 * z5;
    z6 = x4 * y5 - x5 * y4;

    er = sqrt ( x6 * x6 + y6 * y6 + z6 * z6 );
    x6 = x6 / er;
    y6 = y6 / er;
    z6 = z6 / er;

    normal[ 3 * index ] = x6;
    normal[ 3 * index + 1 ] = y6;
    normal[ 3 * index + 2 ] = z6;
}

void ExternalFaces::GetColorIndices(  kvs::AnyValueArray value,
				      const kvs::Real64 min_value,
				      const kvs::Real64 max_value,
				      const size_t colormap_resolution,
				      const kvs::UInt32 node_index[4],
				      kvs::UInt32 (*color_index)[4])
  {
   
    
    const kvs::Real64 normalize =
      static_cast<kvs::Real64>( colormap_resolution - 1 ) / ( max_value - min_value );

    // Scalar data.
    // if ( veclen == 1 )
    {
      for ( size_t i = 0; i < 4; i++ )
        {
	  (*color_index)[i] = kvs::UInt32( normalize * ( kvs::Real64( value.at<double>( node_index[i] ) - min_value ) ));
					   
	}
    }
    /*   // Vector data.
	 else
	 {
	 // In case of the vector component, the magnitude value is calculated.
	 kvs::Real64 magnitude[N]; memset( magnitude, 0, sizeof( kvs::Real64 ) * N );
	 for ( size_t i = 0; i < veclen; ++i )
	 {
	 for ( size_t j = 0; j < N; j++ )
	 {
	 magnitude[j] += kvs::Math::Square( kvs::Real64( value[ veclen * node_index[j] + i ] ) );
	 }
	 }

	 for ( size_t i = 0; i < N; i++ )
	 {
	 magnitude[i] = std::sqrt( magnitude[i] );
	 (*color_index)[i] = kvs::UInt32( normalize * ( magnitude[i] - min_value ) );
	 }
	 }*/
  }
  
} // end of namespace YinYangVis
