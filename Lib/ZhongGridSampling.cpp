#include "ZhongGridSampling.h"
#include "ZhongVolumeObject.h"
#include "ZhongGrid.h"
#include "DensityMap.h"
#include <kvs/Timer>

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
    YinYangVis::ZhongGrid* m_grid;
    YinYangVis::DensityMap* m_density_map;
    Particles m_particles; ///< particles
    Particle m_current; ///< current sampled point
    Particle m_trial; ///< trial point

public:
    Sampler(
        YinYangVis::ZhongGrid* grid,
        YinYangVis::DensityMap* density_map ):
        m_grid( grid ),
        m_density_map( density_map ) {}

    const Particles& particles() const { return m_particles; }

    void bind( const kvs::Vec3ui& base_index )
    {
        m_grid->bind( base_index );
    }

    size_t numberOfParticles(const YinYangVis::ZhongVolumeObject* object )
    {
        const kvs::Real32 scalar = this->averaged_scalar();
        kvs::Real32 overlapweight = averaged_overlapweight( object );
        const kvs::Real32 density = m_density_map->at( scalar , overlapweight );
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

    kvs::Real32 sampleOverlap( const YinYangVis::ZhongVolumeObject* object )
    {
        bool flag = false;
        size_t counter = 0;
        while ( flag == false )
        {
            m_current.coord = m_grid->randomSampling();
            counter++;
            if ( judge_zhong_overlap( object, m_current.coord ) == 0 || counter > 5 ) { flag = true; }
        }
        m_current.normal = -m_grid->gradientVector();
        m_current.scalar = m_grid->scalar();
        return m_density_map->at( m_current.scalar );
    }

    kvs::Real32 sampleOverlap( const size_t max_loops, const YinYangVis::ZhongVolumeObject* object )
    {
        kvs::Real32 density = this->sampleOverlap( object );
        if ( kvs::Math::IsZero( density ) )
        {
            for ( size_t i = 0; i < max_loops; i++ )
            {
                density = this->sampleOverlap( object );
                if ( !kvs::Math::IsZero( density ) ) { break; }
            }
        }

        return density;
    }

    kvs::Real32 trySampleOverlap( const YinYangVis::ZhongVolumeObject* object )
    {
        bool flag = false;
        size_t counter = 0;
        while( flag == false )
        {
            m_trial.coord = m_grid->randomSampling();
            counter++;
            if ( judge_zhong_overlap( object, m_trial.coord ) == 0 || counter > 5 ) { flag = true; }
        }
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

    void checkOverlapFlag( const YinYangVis::ZhongVolumeObject* object, size_t& flag ) const
    {
        kvs::Real32 l;
        size_t cflag[8];
        for ( int i = 0; i < 8; i++ )
        {
            l = sqrt( m_grid->coord(i).x() * m_grid->coord(i).x() + m_grid->coord(i).y() * m_grid->coord(i).y() + m_grid->coord(i).z() * m_grid->coord(i).z() );
            if ( l >= object->rangeR().min )
                cflag[i] = 1;
            else
                cflag[i] = 0;
        }
        flag = 0;
        flag |= cflag[0] * 128;
        flag |= cflag[1] * 64;
        flag |= cflag[2] * 32;
        flag |= cflag[3] * 16;
        flag |= cflag[4] * 8;
        flag |= cflag[5] * 4;
        flag |= cflag[6] * 2;
        flag |= cflag[7] * 1;
    }

private:
    struct Range
    {
        kvs::Real32 min;
        kvs::Real32 max;
        kvs::Real32 d;
    };

    kvs::Real32 calc_zhong_overlap_weight( const YinYangVis::ZhongVolumeObject* object, int i, int j, int k ) const
    {
        kvs::Real32 x, y, z, l;

        const size_t dim = object->dim();
        Range range_r = { object->rangeR().min, object->rangeR().max, object->rangeR().d };

        const kvs::Real32 inner_r = range_r.min + range_r.d * 2;
        const kvs::Real32 inner_max = inner_r;
        const kvs::Real32 inner_min = -inner_r;
        const kvs::Real32 inner_d = ( inner_max - inner_min ) / ( dim - 1 );

        x = inner_min + inner_d * i;
        y = inner_min + inner_d * j;
        z = inner_min + inner_d * k;
        l = sqrt( x*x + y*y + z*z );

        if ( l < inner_r )
        {
            return 1.0f;
        }
        else
        {
            return 0.0f;
        }
    }

    kvs::Real32 averaged_overlapweight( const YinYangVis::ZhongVolumeObject* object ) const
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

    size_t judge_zhong_overlap( const YinYangVis::ZhongVolumeObject* object, kvs::Vec3 coords ) const
    {
        kvs::Real32 l;
        kvs::Real32 outercore_min_r;
        l = sqrt( coords.x() * coords.x() + coords.y() * coords.y() + coords.z() * coords.z() );
        outercore_min_r = object->rangeR().min;
        if ( l >= outercore_min_r )
        {
            return 1;
        }
        else
        {
            return 0;
        }
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


namespace YinYangVis
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

    this->mapping( YinYangVis::ZhongVolumeObject::DownCast( volume ) );

    if ( delete_camera )
    {
        delete m_camera;
        m_camera = 0;
    }

    return this;
}

void ZhongGridSampling::mapping( const YinYangVis::ZhongVolumeObject* volume )
{
    KVS_ASSERT( volume != NULL );

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    YinYangVis::ZhongGrid grid( volume );
    YinYangVis::DensityMap density_map;
    density_map.setSubpixelLevel( m_subpixel_level );
    density_map.setSamplingStep( m_sampling_step );
    density_map.attachCamera( m_camera );
    density_map.attachObject( volume );
    density_map.create( BaseClass::transferFunction().opacityMap() );

    ::Sampler sampler( &grid, &density_map );
    const size_t dim = volume->dim(); // resolution
    const kvs::ColorMap color_map( BaseClass::transferFunction().colorMap() );
    size_t overlap_flag = 0; 
    kvs::Timer timer;
    kvs::Timer timer2;
    float sum=0.0f;
    timer.start();
    for ( size_t k = 0; k < dim - 1; k++ )
    {
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                sampler.bind( kvs::Vec3ui( i, j, k ) );
                timer2.start();
                const size_t nparticles = sampler.numberOfParticles( volume );
                timer2.stop();
                sum += timer2.sec();
                const size_t max_loops = nparticles * 10;
                if ( nparticles == 0 ) continue;

                size_t nduplications = 0;
                size_t counter = 0;
                kvs::Real32 density;
                if( overlap_flag == 0 || overlap_flag == 255 )
                {
                    density =  sampler.sample( max_loops );
                }
                else //1 <=  flag <= 254
                {
                    //std::cout << "checkflag = " << overlap_flag << std::endl;
                    density = sampler.sampleOverlap( max_loops, volume );
                }
                while ( counter < nparticles )
                {
                    kvs::Real32 density_trial;
                    if ( overlap_flag == 0 || overlap_flag == 255 )
                    {
                        density_trial = sampler.trySample();
                    }
                    else //1 <=  flag <= 254
                    {
                        density_trial = sampler.trySampleOverlap( volume );
                    }

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
    timer.stop();
//    std::cout << "calc nparticle time" << sum << std::endl;
//    std::cout << "generate particle time " << timer.sec() - sum << std::endl;
    SuperClass::setCoords( sampler.particles().coords() );
    SuperClass::setColors( sampler.particles().colors() );
    SuperClass::setNormals( sampler.particles().normals() );
    SuperClass::setSize( 1.0f );
}

} // end of namespace YinYangVis
