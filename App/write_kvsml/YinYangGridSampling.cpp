#include "YinYangGridSampling.h"
#include "YinYangVolumeObject.h"
#include "YinYangGrid.h"
#include "DensityMap.h"
#include <kvs/Timer>
#include <string>
#include <math.h>

namespace
{

inline kvs::Real32 RandomNumber()
{
    static kvs::UInt32 x=123456789,y=362436069,z=521288629,w=88675123;
    kvs::UInt32 t;
    t=(x^(x<<11));
    x=y;y=z;z=w;
    w=(w^(w>>19))^(t^(t>>8));
    return w * ( 1.0f / 4294967296.0f );
}

}


namespace
{

struct Particle
{
    kvs::Vec3 coord; ///< coordinate value
    kvs::Vec3 normal; ///< normal vector
    kvs::Real32 scalar; ///< scalar value
};

class Particles
{
private:
    std::vector<kvs::Real32> m_coords; ///< coorinate value array
    std::vector<kvs::Real32> m_normals; ///< normal vector array
    std::vector<kvs::UInt8> m_colors; ///< color value array

public:
    Particles() {}
    kvs::ValueArray<kvs::Real32> coords() const { return kvs::ValueArray<kvs::Real32>( m_coords ); }
    kvs::ValueArray<kvs::Real32> normals() const { return kvs::ValueArray<kvs::Real32>( m_normals ); }
    kvs::ValueArray<kvs::UInt8> colors() const { return kvs::ValueArray<kvs::UInt8>( m_colors ); }
    void push( const Particle& particle, const kvs::ColorMap& color_map )
    {
        const kvs::RGBColor color = color_map.at( particle.scalar );
        m_coords.push_back( particle.coord.x() );
        m_coords.push_back( particle.coord.y() );
        m_coords.push_back( particle.coord.z() );
        m_normals.push_back( particle.normal.x() );
        m_normals.push_back( particle.normal.y() );
        m_normals.push_back( particle.normal.z() );
        m_colors.push_back( color.r() );
        m_colors.push_back( color.g() );
        m_colors.push_back( color.b() );
    }
};

class Sampler
{
private:
    local::YinYangGrid* m_grid;
    local::DensityMap* m_density_map;
    Particles m_particles; ///< particles
    Particle m_current; ///< current sampled point
    Particle m_trial; ///< trial point
  //kvs::Real32 m_overlap_weight; //overlap weight of point( j, k ) : ( j = theta, k = phi )

public:
    Sampler(
        local::YinYangGrid* grid,
        local::DensityMap* density_map ):
        m_grid( grid ),
        m_density_map( density_map ) {}

    const Particles& particles() const { return m_particles; }

    void bind( const kvs::Vec3ui& base_index )
    {
        m_grid->bind( base_index );
    }

  size_t numberOfParticles( const local::YinYangVolumeObject* object )
    {
      const kvs::Real32 scalar = this->averaged_scalar();
	kvs::Real32 density;
	density =  m_density_map->at( scalar );
        const kvs::Real32 volume = m_grid->volume();
        return this->number_of_particles( density, volume );
    }
  
  size_t numberOfParticles( const local::YinYangVolumeObject* object, size_t overlap_flag )
  {
    const kvs::Real32 scalar = this->averaged_scalar();
    kvs::Real32 density;
    kvs::Real32 overlap_weight;
    if( overlap_flag == 0 )
      {
	overlap_weight = 1.0f;
      }
    else if( overlap_flag == 255 )
      {
	overlap_weight = 0.0f;
      }
    else
      {
	//std::cout << "checkflag = " << overlap_flag << std::endl;
	//overlap_weight = 0.0f;
	//overlap_weight = number_of_flags( checkOverlapFlag( object ) ) / 8.0f;
	overlap_weight = calc_yinyang_overlap_weight( object );
	//std::cout << "overlap_weight = " << overlap_weight << std::endl;
      }
    density = m_density_map->at( scalar ) * overlap_weight;
    
    const kvs::Real32 volume = m_grid->volume();
    return this->number_of_particles( density, volume );
  }
  
  
  kvs::Vec3 getGridCoords( size_t index)
    {
      return m_grid->coord( index );
    }

  kvs::Real32 sample()
    {
        m_current.coord = m_grid->randomSampling();
        m_current.normal = -m_grid->gradientVector();
        m_current.scalar = m_grid->scalar();
        return m_density_map->at( m_current.scalar );
    }

    kvs::Real32 sample( const kvs::Vec3& coord, const kvs::Vec3& normal, const kvs::Real32 scalar )
    {
        m_current.coord = coord;
        m_current.normal = normal;
        m_current.scalar = scalar;
        return m_density_map->at( m_current.scalar );
    }

  kvs::Real32 sample( const size_t max_loops )
    {
      kvs::Real32 density = this->sample();
      if ( kvs::Math::IsZero( density ) )
        {
	  for ( size_t i = 0; i < max_loops; i++ )
            {
	      density = this->sample();
	      if ( !kvs::Math::IsZero( density ) ) { break; }
            }
        }

        return density;
    }

  kvs::Real32 trySample()
    {
        m_trial.coord = m_grid->randomSampling();
        m_trial.normal = -m_grid->gradientVector();
        m_trial.scalar = m_grid->scalar();
        return m_density_map->at( m_current.scalar );
    }

  kvs::Real32 sampleOverlap( const local::YinYangVolumeObject* object, size_t& sample_stoper )
    {
      size_t counter = 0;
      while( 1 )
	{
	  m_current.coord = m_grid->randomSampling();
	  counter++;
	  if( judge_yinyang_overlap( object, m_current.coord ) == 0 ) break;
	  else if( counter > 500000 )
	    {
	      sample_stoper = 1;
	      //std::cout << "sample stop." << std::endl;
	      break;
	    }
	}
        m_current.normal = -m_grid->gradientVector();
        m_current.scalar = m_grid->scalar();
        return m_density_map->at( m_current.scalar );
    }

  kvs::Real32 sampleOverlap( const size_t max_loops, const local::YinYangVolumeObject* object, size_t& sample_stoper )
    {
      kvs::Real32 density = this->sampleOverlap( object, sample_stoper );
      if ( kvs::Math::IsZero( density ) )
        {
	  for ( size_t i = 0; i < max_loops; i++ )
            {
	      density = this->sampleOverlap( object, sample_stoper );
	      if ( !kvs::Math::IsZero( density ) ) { break; }
            }
        }

        return density;
    }

  kvs::Real32 trySampleOverlap( const local::YinYangVolumeObject* object, size_t& sample_stoper )
    {
      size_t counter = 0;
      while( 1 )
	{
	  m_trial.coord = m_grid->randomSampling();
	  counter++;
	  if( judge_yinyang_overlap( object, m_trial.coord ) == 0 ) break;
	  else if( counter > 500000 )
	    {
	      sample_stoper = 1;
	      break;
	    }
	}
        m_trial.normal = -m_grid->gradientVector();
        m_trial.scalar = m_grid->scalar();
        return m_density_map->at( m_current.scalar );
    }

    void accept( const kvs::ColorMap& cmap )
    {
        m_particles.push( m_current, cmap );
    }

    void acceptTrial( const kvs::ColorMap& cmap )
    {
        m_particles.push( m_trial, cmap );
        m_current.coord = m_trial.coord;
        m_current.normal = m_trial.normal;
        m_current.scalar = m_trial.scalar;
    }

    
  size_t checkOverlapFlag( const local::YinYangVolumeObject* object ) const
  {
    kvs::Real32 tht_max;
    kvs::Real32 tht_min;
    kvs::Real32 phi_max;
    kvs::Real32 phi_min;
    const size_t dim_theta = object->dimTheta();
    const size_t dim_phi = object->dimPhi();

    RangeYY range_theta = { object->rangeTheta().min, object->rangeTheta().max, object->rangeTheta().d };
    RangeYY range_phi = { object->rangePhi().min, object->rangePhi().max, object->rangePhi().d };

    tht_min = range_theta.min + range_theta.d * ( 0 - 1 );
    tht_max = range_theta.min + range_theta.d * ( ( dim_theta - 1 ) - 1 );
    phi_min = range_phi.min + range_phi.d * ( 0 - 2 );
    phi_max = range_phi.min + range_phi.d * ( ( dim_phi - 1 ) - 2 );

    kvs::Real32 x, y, z;
    kvs::Real32 r, tht, phi;

    size_t flag = 0;
    size_t cflag[8];
    for( int i = 0; i < 8; i++)
      {
	x = m_grid->coord(i).x();
	y = m_grid->coord(i).y();
	z = m_grid->coord(i).z();
	r = sqrt( x*x + y*y + z*z );
	tht = acos( z / r );
	phi = atan2( y, x );

	if( tht_min < tht && tht < tht_max && phi_min < phi && phi < phi_max )
	  {
	    cflag[i] = 1;
	  }
	else
	  {
	    cflag[i] = 0;
	  }
      }
    flag = 0;
    flag += cflag[0] * 128;
    flag += cflag[1] * 64;
    flag += cflag[2] * 32;
    flag += cflag[3] * 16;
    flag += cflag[4] * 8;
    flag += cflag[5] * 4;
    flag += cflag[6] * 2;
    flag += cflag[7] * 1;
    return flag;
  }

  
    
private:

  //---YinYangの重複部分の重み計算---------------
  struct ControlVolume
  {
    kvs::Real32 f11;
    kvs::Real32 f12;
    kvs::Real32 f21;
    kvs::Real32 f22;
    size_t pattern_flag;
  };
  
  struct RangeYY
  {
    kvs::Real32 min;
    kvs::Real32 max;
    kvs::Real32 d;
  };
  
  void rtp_self_to_othr( kvs::Real32 tht_self, kvs::Real32 phi_self, kvs::Real32& tht_othr, kvs::Real32& phi_othr ) const
  {
    kvs::Real32 x_self, y_self, z_self;
    kvs::Real32 x_othr, y_othr, z_othr, r_othr;

    x_self = sin( tht_self ) * cos( phi_self );
    y_self = sin( tht_self ) * sin( phi_self );
    z_self = cos( tht_self );
    x_othr = -x_self;
    y_othr =  z_self;
    z_othr =  y_self;
    r_othr = sqrt( x_othr*x_othr + y_othr*z_othr + z_othr*z_othr );
    tht_othr = acos( z_othr / r_othr );
    phi_othr = atan2( y_othr, x_othr );
  }
  
  size_t overlap_flag( kvs::Real32 f ) const
  {
    size_t flag = 0;

    if( f >= 0.0f )
      flag = 1;
    else
      flag = 0;

    return flag;
  }
  
  kvs::Real32 pyramid( kvs::Real32 tht, kvs::Real32 phi, kvs::Real32 tht_middle, kvs::Real32 phi_middle, kvs::Real32 tht_halfspan, kvs::Real32 phi_halfspan ) const
  {
    kvs::Real32 py;
    kvs::Real32 fabs_x, fabs_y;
    fabs_x = fabs( ( phi - phi_middle ) / phi_halfspan );
    fabs_y = fabs( ( tht - tht_middle ) / tht_halfspan );
    
    py = ( 1 - fabs_x ) * step_function( fabs_x - fabs_y ) + (1 - fabs_y ) * step_function( fabs_y - fabs_x );
    return py;
  }
  
  kvs::Real32 how_deep_in_other_system( const kvs::Real32 tht_self, const kvs::Real32 phi_self, kvs::Real32 tht_middle, kvs::Real32 phi_middle, kvs::Real32 tht_halfspan, kvs::Real32 phi_halfspan ) const
  {
    kvs::Real32 tht_othr, phi_othr;
    
    rtp_self_to_othr( tht_self, phi_self, tht_othr, phi_othr );
    return pyramid( tht_othr, phi_othr, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
  }

  kvs::Real32 overlap_weight( ControlVolume cv ) const
  {
    if( cv.pattern_flag == 0 ) //0000
      {
	//std::cout << "----" << std::endl;
	return 1.0f;
      }

    if( cv.pattern_flag == 15 ) //1111
       {
	 //std::cout << "++++" << std::endl;
	 return 0.0f;
       }

    kvs::Real32 ol_weight = 0.0f;
    kvs::Real32 a = 0.0f, b = 0.0f, s = 0.0f;
     if( cv.pattern_flag == 14 ) //1110
       {
	 a = -cv.f22 / ( cv.f12 - cv.f22 );
	 b = -cv.f22 / ( cv.f21 - cv.f22 );
	 s = a * b / 2.0f;
       }
    else if( cv.pattern_flag == 13 ) //1101
       {
	 a = -cv.f21 / ( cv.f11 - cv.f21 );
	 b = -cv.f21 / ( cv.f22 - cv.f21 );
	 s = a * b / 2.0f;
       }
    else if( cv.pattern_flag == 11 ) //1011
       {
	 a = -cv.f12 / ( cv.f11 - cv.f12 );
	 b = -cv.f12 / ( cv.f22 - cv.f12 );
	 s = a * b / 2.0f;
       }
    else if( cv.pattern_flag == 7 ) //0111
       {
	 a = -cv.f11 / ( cv.f12 - cv.f11 );
	 b = -cv.f11 / ( cv.f21 - cv.f11 );
	 s = a * b / 2.0f;
       }
    else if( cv.pattern_flag == 12 ) //1100
       {
	 a = -cv.f21 / ( cv.f11 - cv.f21 );
	 b = -cv.f22 / ( cv.f12 - cv.f22 );
	 s = ( a + b ) / 2.0f;
       }
    else if( cv.pattern_flag == 3 ) //0011
       {
	 a = -cv.f11 / ( cv.f21 - cv.f11 );
	 b = -cv.f12 / ( cv.f22 - cv.f12 );
	 s = ( a + b ) / 2.0f;
       }
    else if( cv.pattern_flag == 10 ) //1010
       {
	 a = -cv.f12 / ( cv.f11 - cv.f12 );
	 b = -cv.f22 / ( cv.f21 - cv.f22 );
	 s = ( a + b ) / 2.0f;
       }	
    else if( cv.pattern_flag == 5 ) //0101
       {
	 a = -cv.f11 / ( cv.f12 - cv.f11 );
	 b = -cv.f21 / ( cv.f22 - cv.f11 );
	 s = ( a + b ) / 2.0f;
       }
    else if( cv.pattern_flag == 8 ) //1000
       {
	 a = cv.f11 / ( cv.f11 - cv.f12 );
	 b = cv.f11 / ( cv.f11 - cv.f21 );
	 s = 1.0f - a * b / 2.0f;
       }	
    else if( cv.pattern_flag == 4 ) //0100
       {
	 a = cv.f12 / ( cv.f12 - cv.f11 );
	 b = cv.f12 / ( cv.f12 - cv.f22 );
	 s = 1.0f - a * b / 2.0f;
       }	
    else if( cv.pattern_flag == 2 ) //0010
       {
	 a = cv.f21 / ( cv.f21 - cv.f11 );
	 b = cv.f21 / ( cv.f21 - cv.f22 );
	 s = 1.0f - a * b / 2.0f;
       }	
    else if( cv.pattern_flag == 1 ) //0001
       {
	 a = cv.f22 / ( cv.f22 - cv.f12 );
	 b = cv.f22 / ( cv.f22 - cv.f21 );
	 s = 1.0f - a * b / 2.0f;
       }	

     ol_weight = s;// + ( 1.0f - s ) / 2.0f;
    
    return ol_weight;
  }

  kvs::Real32 step_function( kvs::Real32 x ) const
  {
    kvs::Real32 st;
    
    if( x > 0.0f )
      st = 1.0f;
    else if( x < 0.0f )
      st = 0.0f;
    else
      st = 0.5f;

    return st;
  }

  kvs::Real32 calc_yinyang_overlap_weight( const local::YinYangVolumeObject* object ) const
  {
    kvs::Real32 r[5], tht[5], phi[5], x[5], y[5], z[5], tht_n, tht_s, phi_w, phi_e;
    std::string c11, c12, c21, c22;

    kvs::Real32 tht_ctr_max;
    kvs::Real32 tht_ctr_min;
    kvs::Real32 tht_middle;
    kvs::Real32 tht_halfspan;
    kvs::Real32 phi_ctr_max;
    kvs::Real32 phi_ctr_min;
    kvs::Real32 phi_middle;
    kvs::Real32 phi_halfspan;

    const size_t dim_theta = object->dimTheta();
    const size_t dim_phi = object->dimPhi();

    ControlVolume cv;
    RangeYY range_theta = { object->rangeTheta().min, object->rangeTheta().max, object->rangeTheta().d };
    RangeYY range_phi = { object->rangePhi().min, object->rangePhi().max, object->rangePhi().d };

    tht_ctr_min = range_theta.min + range_theta.d * ( kvs::Real32 )( 0 - 1 );
    tht_ctr_max = range_theta.min + range_theta.d * ( kvs::Real32 )( ( dim_theta -1 ) - 1 );
    tht_middle = ( tht_ctr_max + tht_ctr_min ) / 2.0f;
    tht_halfspan = ( tht_ctr_max - tht_ctr_min ) / 2.0f;

    phi_ctr_min = range_phi.min + range_phi.d * ( kvs::Real32 )( 0 - 2 );
    phi_ctr_max = range_phi.min + range_phi.d * ( kvs::Real32 )( ( dim_phi -1 ) - 2 );
    phi_middle = ( phi_ctr_max + phi_ctr_min ) / 2.0f;
    phi_halfspan = ( phi_ctr_max - phi_ctr_min ) / 2.0f;

    for( int i = 0; i < 5; i++ )
      {
	x[i] = m_grid->coord(i).x();
	y[i] = m_grid->coord(i).y();
	z[i] = m_grid->coord(i).z();
	r[i] = sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i] );
	tht[i] = acos( z[i] / r[i] );
	phi[i] = atan2( y[i], x[i] );
      }
    
    tht_n = tht[0];
    tht_s = tht[3];
    phi_w = phi[4];
    phi_e = phi[0];

    cv.f11 = pyramid( tht_n, phi_w, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f12 = pyramid( tht_n, phi_e, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f21 = pyramid( tht_s, phi_w, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f22 = pyramid( tht_s, phi_e, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.pattern_flag = 0;
    cv.pattern_flag += overlap_flag( cv.f11 ) * 8; //1000
    cv.pattern_flag += overlap_flag( cv.f12 ) * 4; //0100
    cv.pattern_flag += overlap_flag( cv.f21 ) * 2; //0010
    cv.pattern_flag += overlap_flag( cv.f22 );     //0001

    return overlap_weight( cv );
  }

  size_t judge_yinyang_overlap( const local::YinYangVolumeObject* object, kvs::Vec3 coords ) const
  {
    kvs::Real32 r = sqrt( coords.x()*coords.x() + coords.y()*coords.y() + coords.z()*coords.z() );
    kvs::Real32 tht = acos( coords.z() / r );
    kvs::Real32 phi = atan2( coords.y(), coords.x() );
    kvs::Real32 tht_max;
    kvs::Real32 tht_min;
    kvs::Real32 phi_max;
    kvs::Real32 phi_min;

    const size_t dim_theta = object->dimTheta();
    const size_t dim_phi = object->dimPhi();

    RangeYY range_theta = { object->rangeTheta().min, object->rangeTheta().max, object->rangeTheta().d };
    RangeYY range_phi = { object->rangePhi().min, object->rangePhi().max, object->rangePhi().d };

    tht_min = range_theta.min + range_theta.d * ( 0 - 1 );
    tht_max = range_theta.min + range_theta.d * ( ( dim_theta - 1 ) - 1 );
    phi_min = range_phi.min + range_phi.d * ( 0 - 2 );
    phi_max = range_phi.min + range_phi.d * ( ( dim_phi - 1 ) - 2 );

    if( tht_min < tht && tht < tht_max && phi_min < phi && phi < phi_max )
      {
	return 1;
      }
    else
      {
	return 0;
      }
  }

  kvs::Real32 number_of_flags( size_t flag ) const
  {
    size_t counter = 0;
    if( ( flag & 1 ) == 1 ) counter++;
    if( ( flag & 2 ) == 2 ) counter++;
    if( ( flag & 4 ) == 4 ) counter++;
    if( ( flag & 8 ) == 8 ) counter++;
    if( ( flag & 16 ) == 16 ) counter++;
    if( ( flag & 32 ) == 32 ) counter++;
    if( ( flag & 64 ) == 64 ) counter++;
    if( ( flag & 128 )== 128 ) counter++;
    return counter;
  }

  kvs::Real32 averaged_scalar() const
  {
    return (
            m_grid->value(0) +
            m_grid->value(1) +
            m_grid->value(2) +
            m_grid->value(3) +
            m_grid->value(4) +
            m_grid->value(5) +
            m_grid->value(6) +
            m_grid->value(7) ) / 8.0f;
  }

  size_t number_of_particles( const kvs::Real32 density, const kvs::Real32 volume ) const
  {
    const kvs::Real32 R = ::RandomNumber();
    const kvs::Real32 N = density * volume;
    size_t n = static_cast<size_t>( N );
    if ( N - n > R ) { ++n; }
    return n;
  }

    };

}


namespace local
{
YinYangGridSampling::YinYangGridSampling(
    const kvs::VolumeObjectBase* volume,
    const size_t subpixel_level,
    const kvs::Real32 sampling_step,
    const kvs::TransferFunction& transfer_function,
    const kvs::Real32 object_depth ):
    kvs::MapperBase( transfer_function ),
    kvs::PointObject(),
    m_camera( 0 )
{
    this->setSubpixelLevel( subpixel_level );
    this->setSamplingStep( sampling_step );
    this->setObjectDepth( object_depth );
    this->exec( volume );
}

YinYangGridSampling::YinYangGridSampling(
    const kvs::Camera* camera,
    const kvs::VolumeObjectBase* volume,
    const size_t subpixel_level,
    const kvs::Real32 sampling_step,
    const kvs::TransferFunction& transfer_function,
    const kvs::Real32 object_depth )
{
    this->attachCamera( camera );
    this->setSubpixelLevel( subpixel_level );
    this->setSamplingStep( sampling_step );
    this->setObjectDepth( object_depth );
    this->exec( volume );
}

YinYangGridSampling::SuperClass* YinYangGridSampling::exec( const kvs::ObjectBase* object )
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
        kvsMessageError("Input object is not volume data.");
        return NULL;
    }

    bool delete_camera = false;
    if ( !m_camera )
    {
        m_camera = new kvs::Camera();
        delete_camera = true;
    }

    const local::YinYangVolumeObject* yin_yang_object = local::YinYangVolumeObject::DownCast( volume );
    if ( yin_yang_object->gridType() == yin_yang_object->gridYin() )
      {
        this->mapping_metro_yin( local::YinYangVolumeObject::DownCast( volume ) );
      }
    else
      {
	this->mapping_metro_yin( local::YinYangVolumeObject::DownCast( volume ) );
      }

    if ( delete_camera )
    {
        delete m_camera;
        m_camera = 0;
    }

    return this;
}

  void YinYangGridSampling::mapping_metro_yin( const local::YinYangVolumeObject* volume ) //Yin-grid Metropolis sampling.
  {
    KVS_ASSERT( volume != NULL );
    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    local::YinYangGrid grid( volume );
    local::DensityMap density_map;
    density_map.setSubpixelLevel( m_subpixel_level );
    density_map.setSamplingStep( m_sampling_step );
    density_map.attachCamera( m_camera );
    density_map.attachObject( volume );
    density_map.create( BaseClass::transferFunction().opacityMap() );

    ::Sampler sampler( &grid, &density_map );
    const size_t dim_r = volume->dimR(); // radius
    const size_t dim_theta = volume->dimTheta(); // latitude
    const size_t dim_phi = volume->dimPhi(); // longitude
    const kvs::ColorMap color_map( BaseClass::transferFunction().colorMap() );
    kvs::Timer timer;
    kvs::Timer timer2;
    float sum=0.0f;
	
    size_t size = ( dim_phi - 1 ) * ( dim_theta - 1 ) * ( dim_r - 1 );
    float* nparticles = new float[size];
    
    timer.start();    
    size_t index = 0;
    for ( size_t k = 0; k < dim_phi - 1; k++ )
      {
    	for ( size_t j = 0; j < dim_theta - 1; j++ )
    	  {
    	    for ( size_t i = 0; i < dim_r - 1; i++, index++ )
    	      {
    		sampler.bind( kvs::Vec3ui( i, j, k ) );
    		nparticles[index] = sampler.numberOfParticles( volume );
    	      }
    	  }
      }
    timer.stop();
    std::cout << "Calculate number of particles." << timer.sec() << " sec" << std::endl;

    index = 0;
    timer2.start();
    for ( size_t k = 0; k < dim_phi - 1; k++ )
      {
    	for ( size_t j = 0; j < dim_theta - 1; j++ )
    	  {
    	    for ( size_t i = 0; i < dim_r - 1; i++, index++ )
    	      {
		sampler.bind( kvs::Vec3ui( i, j, k ) );
    		const size_t max_loops = nparticles[index] * 10;
    		if ( nparticles[index] == 0 ) continue;
    		size_t nduplications = 0;
    		size_t counter = 0;
    		kvs::Real32 density =  sampler.sample( max_loops );
    		while ( counter < nparticles[index] )
    		  {
    		    const kvs::Real32 density_trial = sampler.trySample();
    		    const kvs::Real32 ratio = density_trial / density;
    		    if ( ratio >= 1.0f )
    		      {
    			sampler.acceptTrial( color_map );
    			density = density_trial;
    			counter++;
    		      }
    		    else
    		      {
    			if ( ratio >= ::RandomNumber() )
    			  {
    			    sampler.acceptTrial( color_map );
    			    density = density_trial;
    			    counter++;
    			  }
    			else
    			  {
                            #ifdef DUPLICATION
    			    sampler.accept( color_map );
    			    counter++;
    			    #else
    			    if ( ++nduplications > max_loops ) { break; }
    			    #endif
    			  }
    		      }
    		  } // end of while-loop
    	      } // end of i-loop
    	  } // end of j-loop
      } // end of k-loop
    timer2.stop();
    delete[] nparticles;
    std::cout << std::endl << "Particle generation time for loop: " << timer2.sec() - sum << " [sec]" << std::endl;
    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
 
  }

  void YinYangGridSampling::mapping_metro_yang( const local::YinYangVolumeObject* volume ) //Yang-grid Metropolis sampling.
  {
    KVS_ASSERT( volume != NULL );
    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    local::YinYangGrid grid( volume );
    local::DensityMap density_map;
    density_map.setSubpixelLevel( m_subpixel_level );
    density_map.setSamplingStep( m_sampling_step );
    density_map.attachCamera( m_camera );
    density_map.attachObject( volume );
    density_map.create( BaseClass::transferFunction().opacityMap() );

    ::Sampler sampler( &grid, &density_map );
    const size_t dim_r = volume->dimR(); // radius
    const size_t dim_theta = volume->dimTheta(); // latitude
    const size_t dim_phi = volume->dimPhi(); // longitude
    const kvs::ColorMap color_map( BaseClass::transferFunction().colorMap() );
    kvs::Timer timer;
    kvs::Timer timer2;
	
    size_t size = ( dim_phi - 1 ) * ( dim_theta - 1 ) * ( dim_r - 1 );
    float* nparticles = new float[size];
    size_t* overlap_flag = new size_t[size];
    
    timer.start();    
    size_t index = 0;
    for ( size_t k = 0; k < dim_phi - 1; k++ )
      {
    	for ( size_t j = 0; j < dim_theta - 1; j++ )
    	  {
    	    for ( size_t i = 0; i < dim_r - 1; i++, index++ )
    	      {
    		sampler.bind( kvs::Vec3ui( i, j, k ) );
		overlap_flag[index] = sampler.checkOverlapFlag( volume );
    		nparticles[index] = sampler.numberOfParticles( volume, overlap_flag[index] );
    	      }
    	  }
      }
    timer.stop();
    std::cout << "Calculate number of particles." << timer.sec() << " sec" << std::endl;
    timer2.start();
    index = 0;
    for ( size_t k = 0; k < dim_phi - 1; k++ )
      {
	for ( size_t j = 0; j < dim_theta - 1; j++ )
	  {
	    for ( size_t i = 0; i < dim_r - 1; i++, index++ )
	      {
		sampler.bind( kvs::Vec3ui( i, j, k ) );
		
		const size_t max_loops = nparticles[index] * 10;
		if ( nparticles[index] == 0 ) continue;
		
		size_t nduplications = 0;
		size_t counter = 0;
		size_t sample_stoper = 0;
		kvs::Real32 density;
		kvs::Real32 density_trial;
		kvs::Real32 ratio;
		
		if( overlap_flag[index] == 0 || overlap_flag[index] == 255 ) // overlap flag == 0.
		  {
		    density =  sampler.sample( max_loops );
		    while ( counter < nparticles[index] )
		      {
			density_trial = sampler.trySample();
			ratio = density_trial / density;
			if ( ratio >= 1.0f )
			  {
			    sampler.acceptTrial( color_map );
			    density = density_trial;
			    counter++;
			  }
			else
			  {
			    if ( ratio >= ::RandomNumber() )
			      {
				sampler.acceptTrial( color_map );
				density = density_trial;
				counter++;
			      }
			    else
			      {
                                #ifdef DUPLICATION
				sampler.accept( color_map );
				counter++;
                                #else
				if ( ++nduplications > max_loops ) { break; }
                                #endif
			      }
			  }
		      } // end of while-loop
		    
		  } // end if. overlap flag == 0
		
		else
		  {
		    //std::cout << "checkflag = " << overlap_flag << std::endl;
		    sample_stoper = 0;
		    density =  sampler.sampleOverlap( max_loops, volume, sample_stoper );
		    if( sample_stoper == 1 ) continue;
		    while ( counter < nparticles[index] )
		      {
			sample_stoper = 0;
			density_trial = sampler.trySampleOverlap( volume, sample_stoper );
			if( sample_stoper == 1 ) continue;
			ratio = density_trial / density;
			if ( ratio >= 1.0f )
			  {
			    sampler.acceptTrial( color_map );
			    density = density_trial;
			    counter++;
			  }
			else
			  {
			    if ( ratio >= ::RandomNumber() )
			      {
				sampler.acceptTrial( color_map );
				density = density_trial;
				counter++;
			      }
			    else
			      {
                                #ifdef DUPLICATION
				sampler.accept( color_map );
				counter++;
                                #else
				if ( ++nduplications > max_loops ) { break; }
                                #endif
			      }
			  }
		      } // end of while-loop
		  } //end else. 1 <= overlap flag <= 254
		
	      } // end of i-loop
	  } // end of j-loop
      } // end of k-loop
    timer2.stop();
    delete[] nparticles;
    delete[] overlap_flag;
    std::cout << std::endl << "Particle generation time for loop: " << timer2.sec() << " [sec]" << std::endl;    
    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
  }
  
  void YinYangGridSampling::mapping_uniform( const local::YinYangVolumeObject* volume )
  {
    KVS_ASSERT( volume != NULL );
    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    local::YinYangGrid grid( volume );
    local::DensityMap density_map;
    density_map.setSubpixelLevel( m_subpixel_level );
    density_map.setSamplingStep( m_sampling_step );
    density_map.attachCamera( m_camera );
    density_map.attachObject( volume );
    density_map.create( BaseClass::transferFunction().opacityMap() );

    ::Sampler sampler( &grid, &density_map );
    const size_t dim_r = volume->dimR(); // radius
    const size_t dim_theta = volume->dimTheta(); // latitude
    const size_t dim_phi = volume->dimPhi(); // longitude
    const kvs::ColorMap color_map( BaseClass::transferFunction().colorMap() );
    size_t overlap_flag = 0;
    size_t nparticles = 0;
    size_t sample_stoper = 0;
    
    for ( size_t k = 0; k < dim_phi - 1; k++ )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++ )
        {
	  for ( size_t i = 0; i < dim_r - 1; i++ )
            {
                sampler.bind( kvs::Vec3ui( i, j, k ) );

		overlap_flag = sampler.checkOverlapFlag( volume );

		nparticles = sampler.numberOfParticles( volume );
                if ( nparticles == 0 ) continue;

                for ( size_t i = 0; i < nparticles; ++i )
		  {
		    if( overlap_flag == 0 || overlap_flag == 255 || volume->gridType() == volume->gridYin() )
		      {
			sampler.sample();
			sampler.accept( color_map );
		      }
		    else //1 <=  flag <= 254
		      {
			sample_stoper = 0;
			sampler.sampleOverlap( volume, sample_stoper );
			if( sample_stoper == 1 ) continue;
			sampler.accept( color_map );
		      }
                }
            }
        }
    }

    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
}


} // end of namespace local
