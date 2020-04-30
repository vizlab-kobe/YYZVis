#include "Convert.h"
#include <YYZVis/Lib/YinYangVolumeObjectBase.h>
#include <YYZVis/Lib/ZhongGrid.h>
#include <YYZVis/Lib/YinYangGrid.h>


namespace
{

inline kvs::Vec3 Cart2Polar( const kvs::Vec3 cart )
{
    const float r = cart.length();
    const float t = std::acos( cart.z() / r );
    const float p = std::atan2( cart.y(), cart.x() );
    return kvs::Vec3( r, t, p );
}

inline bool IsInsideOf( const YYZVis::ZhongVolumeObject& zng_volume, const kvs::Vec3 rtp )
{
    return rtp.x() <= zng_volume.rangeR().min;
}

inline bool IsInsideOf( const YYZVis::YinVolumeObject& yin_volume, const kvs::Vec3 rtp )
{
    const float t = rtp.y();
    const float p = rtp.z();
    return
        yin_volume.rangeTheta().min <= t && t <= yin_volume.rangeTheta().max &&
        yin_volume.rangePhi().min <= p && p <= yin_volume.rangePhi().max;
}

inline kvs::Vec3 IndexOf( const YYZVis::YinYangVolumeObjectBase& volume, const kvs::Vec3 rtp )
{
    const kvs::Vec3 min_range( volume.rangeR().min, volume.rangeTheta().min, volume.rangePhi().min );
    const kvs::Vec3 max_range( volume.rangeR().max, volume.rangeTheta().max, volume.rangePhi().max );
    const kvs::Vec3 dim( float( volume.dimR() ), float( volume.dimTheta() ), float( volume.dimPhi() ) );
    return ( dim - kvs::Vec3::Constant( 1.0f ) ) * ( rtp - min_range ) / ( max_range - min_range );
}

inline kvs::Vec3 IndexOf( const YYZVis::ZhongVolumeObject& volume, const kvs::Vec3 xyz )
{
    const float r_min = volume.rangeR().min;
    const float scale = ( volume.dim() - 1.0f ) / ( r_min * 2.0f );
    return ( xyz + kvs::Vec3::Constant( r_min ) ) * scale;
}

}

namespace local
{

kvs::StructuredVolumeObject Convert(
    YYZVis::YinVolumeObject& yin_volume,
    YYZVis::YangVolumeObject& yng_volume,
    YYZVis::ZhongVolumeObject& zng_volume )
{
    const size_t dim = 200;
    const kvs::Vec3ui resolution( dim, dim, dim );
    const size_t veclen = yin_volume.veclen();

    const kvs::Vec3 min_coord = kvs::Vec3::Constant( -yin_volume.rangeR().max );
    const kvs::Vec3 max_coord = kvs::Vec3::Constant(  yin_volume.rangeR().max );
    const kvs::Vec3 d = ( max_coord - min_coord ) / float( dim - 1 );

    YYZVis::YinYangGrid yin_grid( &yin_volume );
    YYZVis::YinYangGrid yng_grid( &yng_volume );
    YYZVis::ZhongGrid zng_grid( &zng_volume );

    kvs::ValueArray<kvs::Real32> values( dim * dim * dim * veclen );
    values.fill(0);
    for ( size_t k = 0, index = 0; k < dim; ++k )
    {
        const float z = min_coord.z() + d.z() * k;
        for ( size_t j = 0; j < dim; ++j )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t i = 0; i < dim; ++i, ++index )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 xyz( x, y, z );
                const kvs::Vec3 rtp = ::Cart2Polar( xyz );

                // Out of region
                if ( rtp.x() > yin_volume.rangeR().max ) { continue; }

                // Inside zhong volume region.
                if ( ::IsInsideOf( zng_volume, rtp ) )
                {
                    const kvs::Vec3 zng_index = ::IndexOf( zng_volume, xyz );
                    const kvs::Vec3ui zng_base_index( zng_index );
                    const kvs::Vec3 zng_local_index = zng_index - kvs::Vec3( zng_base_index );
                    zng_grid.bind( zng_base_index );
                    zng_grid.setLocalPoint( zng_local_index );
                    values[index] = zng_grid.scalar();
                }
                // Inside yin volume region.
                else if ( ::IsInsideOf( yin_volume, rtp ) )
                {
                    const kvs::Vec3 yin_index = ::IndexOf( yin_volume, rtp );
                    const kvs::Vec3ui yin_base_index( yin_index );
                    const kvs::Vec3 yin_local_index = yin_index - kvs::Vec3( yin_base_index );
                    yin_grid.bind( yin_base_index );
                    yin_grid.setLocalPoint( yin_local_index );
                    values[index] = yin_grid.scalar();
                }
                // Inside yang volume region.
                else
                {
                    const kvs::Vec3 yng_index = ::IndexOf( yng_volume, ::Cart2Polar( kvs::Vec3( -x, z, y ) ) );
                    const kvs::Vec3ui yng_base_index( yng_index );
                    const kvs::Vec3 yng_local_index = yng_index - kvs::Vec3( yng_base_index );
                    yng_grid.bind( yng_base_index );
                    yng_grid.setLocalPoint( yng_local_index );
                    values[index] = yng_grid.scalar();
                }
            }
        }
    }

    kvs::StructuredVolumeObject object;
    object.setGridTypeToUniform();
    object.setVeclen( veclen );
    object.setResolution( resolution );
    object.setValues( values );
    return object;
}

} // end of namespace local
