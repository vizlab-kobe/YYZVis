#include "ZhongGridSampling.h"
#include "ZhongVolumeObject.h"
#include "ZhongGrid.h"
#include "DensityMap.h"

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
    local::ZhongGrid* m_grid;
    local::DensityMap* m_density_map;
    Particles m_particles; ///< particles
    Particle m_current; ///< current sampled point
    Particle m_trial; ///< trial point

public:
    Sampler(
        local::ZhongGrid* grid,
        local::DensityMap* density_map ):
        m_grid( grid ),
        m_density_map( density_map ) {}

    const Particles& particles() const { return m_particles; }

    void bind( const kvs::Vec3ui& base_index )
    {
        m_grid->bind( base_index );
    }

    size_t numberOfParticles(const local::ZhongVolumeObject* object )
    {
        const kvs::Real32 scalar = this->averaged_scalar();
       	kvs::Real32 overlapweight = averaged_overlapweight( object );
        const kvs::Real32 density = m_density_map->at( scalar, overlapweight );
        const kvs::Real32 volume = m_grid->volume();
        return this->number_of_particles( density, volume );
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
        return m_density_map->at( m_trial.scalar );
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
  struct Range
  {
    kvs::Real32 min;
    kvs::Real32 max;
    kvs::Real32 d;
  };
  
  kvs::Real32 calc_zhong_overlap_weight( const local::ZhongVolumeObject* object, int i, int j, int k ) const
  {
    kvs::Real32 x, y, z;
    kvs::Real32 rad;

    const size_t dim = object->dim();
    Range range_r = { object->rangeR().min, object->rangeR().max, object->rangeR().d };

    const kvs::Real32 inner_r = range_r.min + range_r.d * 2;
    const kvs::Real32 inner_max = inner_r;
    const kvs::Real32 inner_min = -inner_r;
    const kvs::Real32 inner_d = ( inner_max - inner_min ) / ( dim - 1 );
    
    x = inner_min + inner_d * i;
    y = inner_min + inner_d * j;
    z = inner_min + inner_d * k;
    rad = sqrt( x*x + y*y + z*z );

    if( rad < inner_r )
      {
	return 1.0f;
      }
    else
      {
	return 0.0f;
      }
  }

  kvs::Real32 averaged_overlapweight( const local::ZhongVolumeObject* object ) const
  {
    const kvs::Vec3ui base_index = m_grid->baseIndex();
    int i = (int)base_index[0];
    int j = (int)base_index[1];
    int k = (int)base_index[2];
    kvs::Real32 ow = ( calc_zhong_overlap_weight( object, i, j, k ) +
		       calc_zhong_overlap_weight( object, i + 1, j, k ) +
		       calc_zhong_overlap_weight( object, i + 1, j + 1, k ) +
		       calc_zhong_overlap_weight( object, i, j + 1, k) +
		       calc_zhong_overlap_weight( object, i, j, k + 1 ) +
		       calc_zhong_overlap_weight( object, i + 1, j, k + 1 ) +
		       calc_zhong_overlap_weight( object, i + 1, j + 1, k + 1 ) +
		       calc_zhong_overlap_weight( object, i, j + 1, k + 1 ) ) / 8.0f;
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

ZhongGridSampling::ZhongGridSampling(
    const kvs::VolumeObjectBase* volume,
    const size_t subpixel_level,
    const float sampling_step,
    const kvs::TransferFunction& transfer_function,
    const float object_depth ):
    kvs::MapperBase( transfer_function ),
    kvs::PointObject(),
    m_camera( 0 )
{
    this->setSubpixelLevel( subpixel_level );
    this->setSamplingStep( sampling_step );
    this->setObjectDepth( object_depth );
    this->exec( volume );
}

ZhongGridSampling::ZhongGridSampling(
    const kvs::Camera* camera,
    const kvs::VolumeObjectBase* volume,
    const size_t subpixel_level,
    const float sampling_step,
    const kvs::TransferFunction& transfer_function,
    const float object_depth )
{
    this->attachCamera( camera );
    this->setSubpixelLevel( subpixel_level );
    this->setSamplingStep( sampling_step );
    this->setObjectDepth( object_depth );
    this->exec( volume );
}

ZhongGridSampling::SuperClass* ZhongGridSampling::exec( const kvs::ObjectBase* object )
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

    this->mapping( local::ZhongVolumeObject::DownCast( volume ) );

    if ( delete_camera )
    {
        delete m_camera;
        m_camera = 0;
    }

    return this;
}

void ZhongGridSampling::mapping( const local::ZhongVolumeObject* volume )
{
    KVS_ASSERT( volume != NULL );

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    local::ZhongGrid grid( volume );
    local::DensityMap density_map;
    density_map.setSubpixelLevel( m_subpixel_level );
    density_map.setSamplingStep( m_sampling_step );
    density_map.attachCamera( m_camera );
    density_map.attachObject( volume );
    density_map.create( BaseClass::transferFunction().opacityMap() );

    ::Sampler sampler( &grid, &density_map );
    const size_t dim = volume->dim(); // resolution
    const kvs::ColorMap color_map( BaseClass::transferFunction().colorMap() );

    for ( size_t k = 0; k < dim - 1; k++ )
    {
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                sampler.bind( kvs::Vec3ui( i, j, k ) );

                const size_t nparticles = sampler.numberOfParticles( volume );
                const size_t max_loops = nparticles * 10;
                if ( nparticles == 0 ) continue;

                size_t nduplications = 0;
                size_t counter = 0;
                kvs::Real32 density = sampler.sample( max_loops );
                while ( counter < nparticles )
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

    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
}

} // end of namespace local
