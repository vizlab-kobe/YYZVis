#include "Program.h"
#include "Input.h"
#include "Convert.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
#include <kvs/Indent>


namespace
{

YinYangVis::YinYangVolumeObject ImportYinVolume( const local::Input& input )
{
    const size_t dim_rad = input.dim_rad;
    const size_t dim_lat = input.dim_lat;
    const size_t dim_lon = input.dim_lon;
    const std::string& filename = input.filename_yin;

    YinYangVis::YinYangVolumeObject object;
    object.setGridTypeToYin();
    object.setDimR( dim_rad );
    object.setDimTheta( dim_lat );
    object.setDimPhi( dim_lon );
    object.setVeclen( 1 );
    object.calculateCoords();
    object.readValues( filename );
    object.updateMinMaxCoords();
    object.updateMinMaxValues();
    return object;
}

YinYangVis::YinYangVolumeObject ImportYangVolume( const local::Input& input )
{
    const size_t dim_rad = input.dim_rad;
    const size_t dim_lat = input.dim_lat;
    const size_t dim_lon = input.dim_lon;
    const std::string& filename = input.filename_yang;

    YinYangVis::YinYangVolumeObject object;
    object.setGridTypeToYang();
    object.setDimR( dim_rad );
    object.setDimTheta( dim_lat );
    object.setDimPhi( dim_lon );
    object.setVeclen( 1 );
    object.calculateCoords();
    object.readValues( filename );
    object.updateMinMaxCoords();
    object.updateMinMaxValues();
    return object;
}

YinYangVis::ZhongVolumeObject ImportZhongVolume( const local::Input& input )
{
    const size_t dim_rad = input.dim_rad;
    const size_t dim_zhong = input.dim_zhong;
    const std::string& filename = input.filename_zhong;

    YinYangVis::ZhongVolumeObject object;
    object.setDimR( dim_rad );
    object.setDim( dim_zhong );
    object.setVeclen( 1 );
    object.calculateCoords();
    object.readValues( filename );
    object.updateMinMaxCoords();
    object.updateMinMaxValues();
    return object;
}

void UpdateMinMaxValues(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume )
{
    const kvs::Real32 min_value0 = yin_volume.minValue();
    const kvs::Real32 min_value1 = yang_volume.minValue();
    const kvs::Real32 min_value2 = zhong_volume.minValue();
    const kvs::Real32 max_value0 = yin_volume.maxValue();
    const kvs::Real32 max_value1 = yang_volume.maxValue();
    const kvs::Real32 max_value2 = zhong_volume.maxValue();
    const kvs::Real32 min_value = kvs::Math::Min( min_value0, min_value1, min_value2 );
    const kvs::Real32 max_value = kvs::Math::Max( max_value0, max_value1, max_value2 );
    yin_volume.setMinMaxValues( min_value, max_value );
    yang_volume.setMinMaxValues( min_value, max_value );
    zhong_volume.setMinMaxValues( min_value, max_value );
}

void UpdateMinMaxCoords(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume )
{
    const kvs::Vec3& min_coord0 = yin_volume.minObjectCoord();
    const kvs::Vec3& min_coord1 = yang_volume.minObjectCoord();
    const kvs::Vec3& min_coord2 = zhong_volume.minObjectCoord();
    const kvs::Vec3& max_coord0 = yin_volume.maxObjectCoord();
    const kvs::Vec3& max_coord1 = yang_volume.maxObjectCoord();
    const kvs::Vec3& max_coord2 = zhong_volume.maxObjectCoord();
    const kvs::Real32 min_x = kvs::Math::Min( min_coord0.x(), min_coord1.x(), min_coord2.x() );
    const kvs::Real32 min_y = kvs::Math::Min( min_coord0.y(), min_coord1.y(), min_coord2.y() );
    const kvs::Real32 min_z = kvs::Math::Min( min_coord0.z(), min_coord1.z(), min_coord2.z() );
    const kvs::Real32 max_x = kvs::Math::Min( max_coord0.x(), max_coord1.x(), max_coord2.x() );
    const kvs::Real32 max_y = kvs::Math::Min( max_coord0.y(), max_coord1.y(), max_coord2.y() );
    const kvs::Real32 max_z = kvs::Math::Min( max_coord0.z(), max_coord1.z(), max_coord2.z() );
    const kvs::Vec3 min_coord( min_x, min_y, min_z );
    const kvs::Vec3 max_coord( max_x, max_y, max_z );
    yin_volume.setMinMaxObjectCoords( min_coord, max_coord );
    yin_volume.setMinMaxExternalCoords( min_coord, max_coord );
    yang_volume.setMinMaxObjectCoords( min_coord, max_coord );
    yang_volume.setMinMaxExternalCoords( min_coord, max_coord );
    zhong_volume.setMinMaxObjectCoords( min_coord, max_coord );
    zhong_volume.setMinMaxExternalCoords( min_coord, max_coord );
}

} // end of namespace


namespace local
{

int Program::exec( int argc, char** argv )
{
    local::Input input( argc, argv );
    if ( !input.parse() ) { return 1; }

    std::cout << "IMPORT VOLUMES ..." << std::endl;
    YinYangVis::YinYangVolumeObject yin_volume = ::ImportYinVolume( input );
    YinYangVis::YinYangVolumeObject yang_volume = ::ImportYangVolume( input );
    YinYangVis::ZhongVolumeObject zhong_volume = ::ImportZhongVolume( input );
    ::UpdateMinMaxValues( yin_volume, yang_volume, zhong_volume );
    ::UpdateMinMaxCoords( yin_volume, yang_volume, zhong_volume );

    const kvs::Indent indent( 4 );
    yin_volume.print( std::cout << "YIN VOLUME DATA" << std::endl, indent );
    yang_volume.print( std::cout << "YANG VOLUME DATA" << std::endl, indent );
    zhong_volume.print( std::cout << "ZHONG VOLUME DATA" << std::endl, indent );

    std::cout << "CONVERT VOLUMES ..." << std::endl;
    const bool ascii = false;
    kvs::StructuredVolumeObject cart_volume = local::Convert( yin_volume, yang_volume, zhong_volume );
    cart_volume.print( std::cout << "STRUCTURED VOLUME DATA" << std::endl, indent );
    // cart_volume.write( input.filename_output, ascii );

    return 0;
}

} // end of namespace local
