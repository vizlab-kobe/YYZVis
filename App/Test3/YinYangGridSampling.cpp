#include "YinYangGridSampling.h"
#include "YinYangVolumeObject.h"
#include "YinYangGrid.h"
#include "DensityMap.h"
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
    kvs::Real32 m_overlap_weight; //overlap weight of point( j, k ) : ( j = theta, k = phi )

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
	kvs::Real32 overlapweight = averaged_overlapweight( object );
        const kvs::Real32 density = m_density_map->at( scalar, overlapweight );
        const kvs::Real32 volume = m_grid->volume();
        return this->number_of_particles( density, volume );
    }

  kvs::Vec3 getGridCoords( size_t index)
    {
      return m_grid->coord( index );
    }

  kvs::Real32 sample() // const local::YinYangVolumeObject* object )
    {
        m_current.coord = m_grid->randomSampling();
        m_current.normal = -m_grid->gradientVector();
        m_current.scalar = m_grid->scalar();
	//kvs::Real32 overlapweight = averaged_overlapweight( object );
        return m_density_map->at( m_current.scalar ); //, overlapweight );
    }

    kvs::Real32 sample( const kvs::Vec3& coord, const kvs::Vec3& normal, const kvs::Real32 scalar )
    {
        m_current.coord = coord;
        m_current.normal = normal;
        m_current.scalar = scalar;
        return m_density_map->at( m_current.scalar );
    }

  kvs::Real32 sample( /*const local::YinYangVolumeObject* object,*/ const size_t max_loops )
    {
      kvs::Real32 density = this->sample();// object );
      if ( kvs::Math::IsZero( density ) )
        {
	  for ( size_t i = 0; i < max_loops; i++ )
            {
	      density = this->sample();// object );
	      if ( !kvs::Math::IsZero( density ) ) { break; }
            }
        }

        return density;
    }

  kvs::Real32 trySample() // const local::YinYangVolumeObject* object )
    {
        m_trial.coord = m_grid->randomSampling();
        m_trial.normal = -m_grid->gradientVector();
        m_trial.scalar = m_grid->scalar();
	//kvs::Real32 overlapweight = averaged_overlapweight( object );
        return m_density_map->at( m_current.scalar ); //, overlapweight );
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
    kvs::Real32 x_othr, y_othr, z_othr;

    x_self = sin( tht_self ) * cos( phi_self );
    y_self = sin( tht_self ) * sin( phi_self );
    z_self = cos( tht_self );
    x_othr = -x_self;
    y_othr =  z_self;
    z_othr =  y_self;
    tht_othr = acos( z_othr );
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
    
    py = ( 1 - fabs_x ) * step( fabs_x - fabs_y ) + (1 - fabs_y ) * step( fabs_y - fabs_x );
    return py;
  }
  
  
  
  kvs::Real32 how_deep_in_other_system( const kvs::Real32 tht_self, const kvs::Real32 phi_self, kvs::Real32 tht_middle, kvs::Real32 phi_middle, kvs::Real32 tht_halfspan, kvs::Real32 phi_halfspan ) const
  {
    kvs::Real32 how_deep;
    kvs::Real32 tht_othr, phi_othr;
    
    rtp_self_to_othr( tht_self, phi_self, tht_othr, phi_othr );
    how_deep = pyramid( tht_othr, phi_othr, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    
    return how_deep;
  }

  kvs::Real32 overlap_weight( ControlVolume cv ) const
  {
    kvs::Real32 ol_weight;
    kvs::Real32 a, b, s;
    if( cv.pattern_flag == 14 ) //1110
       {
	 a = -cv.f22;
	 b = -cv.f22;
	 s = a * b / 2;
       }
    else if( cv.pattern_flag == 13 ) //1101
       {
	 a = -cv.f21;
	 b = -cv.f21;
	 s = a * b / 2;
       }
    else if( cv.pattern_flag == 11 ) //1011
       {
	 a = -cv.f12;
	 b = -cv.f12;
	 s = a * b / 2;
       }
    else if( cv.pattern_flag == 7 ) //0111
       {
	 a = -cv.f11;
	 b = -cv.f11;
	 s = a * b / 2;
       }
    else if( cv.pattern_flag == 12 ) //1100
       {
	 a = -cv.f21;
	 b = -cv.f22;
	 s = ( a + b ) / 2;
       }
    else if( cv.pattern_flag == 3 ) //0011
       {
	 a = -cv.f11;
	 b = -cv.f12;
	 s = ( a + b ) / 2;
       }
    else if( cv.pattern_flag == 10 ) //1010
       {
	 a = -cv.f12;
	 b = -cv.f22;
	 s = ( a + b ) / 2;
       }	
    else if( cv.pattern_flag == 5 ) //0101
       {
	 a = -cv.f11;
	 b = -cv.f21;
	 s = ( a + b ) / 2;
       }
    else if( cv.pattern_flag == 8 ) //1000
       {
	 a = -cv.f11;
	 b = -cv.f11;
	 s = 1 - a * b / 2;
       }	
    else if( cv.pattern_flag == 4 ) //0100
       {
	 a = -cv.f12;
	 b = -cv.f12;
	 s = 1 - a * b / 2;
       }	
    else if( cv.pattern_flag == 2 ) //0010
       {
	 a = -cv.f21;
	 b = -cv.f21;
	 s = 1 - a * b / 2;
       }	
    else if( cv.pattern_flag == 1 ) //0001
       {
	 a = -cv.f22;
	 b = -cv.f22;
	 s = 1 - a * b / 2;
       }	
    else if( cv.pattern_flag == 15 ) //1111
       {
	s = 0.0f;
       }	
    else if( cv.pattern_flag == 0 ) //0000
       {
	s = 1.0f;
       }
     else
       {
	return false;
       }

    ol_weight = s + ( 1 - s ) / 2;
    
    return ol_weight;
  }

  kvs::Real32 step( kvs::Real32 x ) const
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

  void debug( std::string x, kvs::Real32 y ) const
  {
    std::cout << x << " = " << y << std::endl;
  }

  void debug( std::string x, size_t y ) const
  {
    std::cout << x << " = " << y << std::endl;
  }	
  
  kvs::Real32 calc_yinyang_overlap_weight( const local::YinYangVolumeObject* object, int j, int k ) const
  {
    kvs::Real32 tht, phi, tht_n, tht_s, phi_w, phi_e;
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

    tht_ctr_min = range_theta.min + range_theta.d * ( 0 - 1 ) + range_theta.d * 0.5f;
    tht_ctr_max = range_theta.min + range_theta.d * ( dim_theta - 1 ) + range_theta.d * 0.5f;
    tht_middle = ( tht_ctr_max + tht_ctr_min ) / 2;
    tht_halfspan = ( tht_ctr_max - tht_ctr_min ) / 2;

    phi_ctr_min = range_phi.min + range_phi.d * ( 0 - 2 ) + range_phi.d * 0.5f;
    phi_ctr_max = range_phi.min + range_phi.d * ( dim_phi - 2 ) + range_phi.d * 0.5f;
    phi_middle = ( phi_ctr_max + phi_ctr_min ) / 2;
    phi_halfspan = ( phi_ctr_max - phi_ctr_min ) / 2;

    if( j == 0 && 1 <= k && k <= (int)dim_phi - 2 )
      j = 1;
    else if( j == (int)dim_theta - 1 && 1 <= k && k <= (int)dim_phi - 2 )
      j = dim_theta - 2;
    else if( 0 <= j && j <= (int)dim_theta - 1 && k == 0 )
      k = 1;
    else if( 0 <= j && j <= (int)dim_theta - 1 && k == (int)dim_phi - 1 )
      k = dim_phi - 2;
    
    tht = range_theta.min + range_theta.d * ( j - 1 );
    phi = range_phi.min + range_phi.d * ( k - 2 );
    tht_n = tht - range_theta.d / 2;
    tht_s = tht + range_theta.d / 2;
    phi_w = phi - range_phi.d / 2;
    phi_e = phi + range_phi.d / 2;
    cv.f11 = how_deep_in_other_system( tht_n, phi_w, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f12 = how_deep_in_other_system( tht_n, phi_e, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f21 = how_deep_in_other_system( tht_s, phi_w, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.f22 = how_deep_in_other_system( tht_s, phi_e, tht_middle, phi_middle, tht_halfspan, phi_halfspan );
    cv.pattern_flag |= overlap_flag( cv.f11 ) * 8; //1000
    cv.pattern_flag |= overlap_flag( cv.f12 ) * 4; //0100
    cv.pattern_flag |= overlap_flag( cv.f21 ) * 2; //0010
    cv.pattern_flag |= overlap_flag( cv.f22 );     //0001
    //std::cout << "( j, k ) = ( " << j << ", " << k << " ) cv.pattern_flag = " << cv.pattern_flag << std::endl;

    return overlap_weight( cv );
  }

  kvs::Real32 averaged_overlapweight( const local::YinYangVolumeObject* object ) const
  {
    const kvs::Vec3ui base_index = m_grid->baseIndex();
    int j = (int)base_index[1];
    int k = (int)base_index[2];
    kvs::Real32 ow = ( calc_yinyang_overlap_weight( object, j, k ) +
		       calc_yinyang_overlap_weight( object, j, k + 1 ) +
		       calc_yinyang_overlap_weight( object, j + 1, k ) +
		       calc_yinyang_overlap_weight( object, j + 1, k + 1 ) ) / 4.0f;
    return ow;

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

    this->mapping( local::YinYangVolumeObject::DownCast( volume ) );

    if ( delete_camera )
    {
        delete m_camera;
        m_camera = 0;
    }

    return this;
}

void YinYangGridSampling::mapping( const local::YinYangVolumeObject* volume )
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

    for ( size_t k = 0; k < dim_phi - 1; k++ )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++ )
        {
	  for ( size_t i = 0; i < dim_r - 1; i++ )
            {
                sampler.bind( kvs::Vec3ui( i, j, k ) );
		const size_t nparticles = sampler.numberOfParticles( volume );
		const size_t max_loops = nparticles * 10;
		if ( nparticles == 0 ) continue;

                size_t nduplications = 0;
                size_t counter = 0;
                kvs::Real32 density = sampler.sample( /*volume,*/ max_loops );
                while ( counter < nparticles )
                {
		  const kvs::Real32 density_trial = sampler.trySample(); // volume );
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

    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
}

} // end of namespace local