#include "UniformGridMerger.h"
#include "YinYangVolumeObjectBase.h"
#include "ZhongGrid.h"
#include "YinYangGrid.h"


namespace
{

inline kvs::Vec3 Cart2Polar( const kvs::Vec3 cart )
{
    const float r = cart.length();
    const float t = std::acos( cart.z() / r );
    const float p = std::atan2( cart.y(), cart.x() );
    return kvs::Vec3( r, t, p );
}

inline bool IsInsideOf( const YYZVis::ZhongVolumeObject* zng_volume, const kvs::Vec3 rtp )
{
    return rtp.x() <= zng_volume->rangeR().min;
}

inline bool IsInsideOf( const YYZVis::YinVolumeObject* yin_volume, const kvs::Vec3 rtp )
{
    const float t = rtp.y();
    const float p = rtp.z();
    return
        yin_volume->rangeTheta().min <= t && t <= yin_volume->rangeTheta().max &&
        yin_volume->rangePhi().min <= p && p <= yin_volume->rangePhi().max;
}

inline kvs::Vec3 IndexOf( const YYZVis::YinYangVolumeObjectBase* volume, const kvs::Vec3 rtp )
{
    const kvs::Vec3 min_range( volume->rangeR().min, volume->rangeTheta().min, volume->rangePhi().min );
    const kvs::Vec3 max_range( volume->rangeR().max, volume->rangeTheta().max, volume->rangePhi().max );
    const kvs::Vec3 dim( float( volume->dimR() ), float( volume->dimTheta() ), float( volume->dimPhi() ) );
    return ( dim - kvs::Vec3::Constant( 1.0f ) ) * ( rtp - min_range ) / ( max_range - min_range );
}

inline kvs::Vec3 IndexOf( const YYZVis::ZhongVolumeObject* volume, const kvs::Vec3 xyz )
{
    const float r_min = volume->rangeR().min;
    const float scale = ( volume->dim() - 1.0f ) / ( r_min * 2.0f );
    return ( xyz + kvs::Vec3::Constant( r_min ) ) * scale;
}

} // end of namespace


namespace YYZVis
{

UniformGridMerger::SuperClass* UniformGridMerger::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input object is NULL." << std::endl;
        return NULL;
    }

    if ( const YinVolumeObject* yin_volume = YinVolumeObject::DownCast( object ) )
    {
        if ( m_yin_volume != yin_volume ) { m_yin_volume = yin_volume; }
    }

    if ( const YngVolumeObject* yng_volume = YngVolumeObject::DownCast( object ) )
    {
        if ( m_yng_volume != yng_volume ) { m_yng_volume = yng_volume; }
    }

    if ( const ZngVolumeObject* zng_volume = ZngVolumeObject::DownCast( object ) )
    {
        if ( m_zng_volume != zng_volume ) { m_zng_volume = zng_volume; }
    }

    const kvs::Vec3ui resolution( m_dim, m_dim, m_dim );
    const size_t veclen = m_yin_volume->veclen(); // Supported only veclen of 1

    const kvs::Vec3 min_coord = kvs::Vec3::Constant( -m_yin_volume->rangeR().max );
    const kvs::Vec3 max_coord = kvs::Vec3::Constant(  m_yin_volume->rangeR().max );
    const kvs::Vec3 d = ( max_coord - min_coord ) / float( m_dim - 1 );

    YYZVis::YinYangGrid yin_grid( m_yin_volume );
    YYZVis::YinYangGrid yng_grid( m_yng_volume );
    YYZVis::ZhongGrid zng_grid( m_zng_volume );

    kvs::ValueArray<kvs::Real32> values( m_dim * m_dim * m_dim * veclen );
    values.fill(0);
    for ( size_t k = 0, index = 0; k < m_dim; ++k )
    {
        const float z = min_coord.z() + d.z() * k;
        for ( size_t j = 0; j < m_dim; ++j )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t i = 0; i < m_dim; ++i, ++index )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 xyz( x, y, z );
                const kvs::Vec3 rtp = ::Cart2Polar( xyz );

                // Out of region
                if ( rtp.x() > m_yin_volume->rangeR().max ) { continue; }

                // Inside zhong volume region.
                if ( ::IsInsideOf( m_zng_volume, rtp ) )
                {
                    const kvs::Vec3 zng_index = ::IndexOf( m_zng_volume, xyz );
                    const kvs::Vec3ui zng_base_index( zng_index );
                    const kvs::Vec3 zng_local_index = zng_index - kvs::Vec3( zng_base_index );
                    zng_grid.bind( zng_base_index );
                    zng_grid.setLocalPoint( zng_local_index );
                    values[index] = zng_grid.scalar();
                }
                // Inside yin volume region.
                else if ( ::IsInsideOf( m_yin_volume, rtp ) )
                {
                    const kvs::Vec3 yin_index = ::IndexOf( m_yin_volume, rtp );
                    const kvs::Vec3ui yin_base_index( yin_index );
                    const kvs::Vec3 yin_local_index = yin_index - kvs::Vec3( yin_base_index );
                    yin_grid.bind( yin_base_index );
                    yin_grid.setLocalPoint( yin_local_index );
                    values[index] = yin_grid.scalar();
                }
                // Inside yang volume region.
                else
                {
                    const kvs::Vec3 yng_index = ::IndexOf( m_yng_volume, ::Cart2Polar( kvs::Vec3( -x, z, y ) ) );
                    const kvs::Vec3ui yng_base_index( yng_index );
                    const kvs::Vec3 yng_local_index = yng_index - kvs::Vec3( yng_base_index );
                    yng_grid.bind( yng_base_index );
                    yng_grid.setLocalPoint( yng_local_index );
                    values[index] = yng_grid.scalar();
                }
            }
        }
    }

    SuperClass::setGridTypeToUniform();
    SuperClass::setVeclen( veclen );
    SuperClass::setResolution( resolution );
    SuperClass::setValues( values );

    return this;
}

} // end of namespace YYZVis
