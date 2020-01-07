#include <string>
#include <kvs/Math>
#include <kvs/CommandLine>
#include <kvs/TransferFunction>
#include <kvs/DivergingColorMap>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/Bounds>
#include <kvs/PolygonRenderer>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/Axis3D>
#include <kvs/ColorMapBar>
#include <kvs/Label>
#include <kvs/Slider>
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>
#include "Isosurface.h"
#include "SlicePlane.h"


int main( int argc, char** argv )
{
    // Parse commandline options.
    kvs::CommandLine commandline( argc, argv );
    commandline.addOption( "dim_rad", "Dimension in the radial direction.", 1, true );
    commandline.addOption( "dim_lat", "Dimension in the latitude direction.", 1, true );
    commandline.addOption( "dim_lon", "Dimension in the longitude direction.", 1, true );
    commandline.addOption( "dim_zhong", "Dimension of the cubic zhong grid.", 1, true );
    commandline.addOption( "yin", "Filename of yin volume data.", 1, true );
    commandline.addOption( "yang", "Filename of yang volume data.", 1, true );
    commandline.addOption( "zhong", "Filename of zhong volume data.", 1, true );
    commandline.addHelpOption();
    if ( !commandline.parse() ) { return 1; }

    // Read commandline options.
    const std::string filename_yin = commandline.optionValue<std::string>("yin");
    const std::string filename_yang = commandline.optionValue<std::string>("yang");
    const std::string filename_zhong = commandline.optionValue<std::string>("zhong");
    const size_t dim_rad = commandline.optionValue<size_t>("dim_rad");
    const size_t dim_lat = commandline.optionValue<size_t>("dim_lat");
    const size_t dim_lon = commandline.optionValue<size_t>("dim_lon");
    const size_t dim_zhong = commandline.optionValue<size_t>("dim_zhong");
    const kvs::TransferFunction tfunc( kvs::DivergingColorMap::CoolWarm( 256 ) );

    // Read Yin grid data.
    YinYangVis::YinYangVolumeObject yin_volume;
    yin_volume.setGridTypeToYin();
    yin_volume.setDimR( dim_rad );
    yin_volume.setDimTheta( dim_lat );
    yin_volume.setDimPhi( dim_lon );
    yin_volume.setVeclen( 1 );
    yin_volume.calculateCoords();
    yin_volume.readValues( filename_yin );
    yin_volume.updateMinMaxCoords();
    yin_volume.updateMinMaxValues();

    // Read Yang grid data.
    YinYangVis::YinYangVolumeObject yang_volume;
    yang_volume.setGridTypeToYang();
    yang_volume.setDimR( dim_rad );
    yang_volume.setDimTheta( dim_lat );
    yang_volume.setDimPhi( dim_lon );
    yang_volume.setVeclen( 1 );
    yang_volume.calculateCoords();
    yang_volume.readValues( filename_yang );
    yang_volume.updateMinMaxCoords();
    yang_volume.updateMinMaxValues();

    // Read Zhong grid data.
    YinYangVis::ZhongVolumeObject zhong_volume;
    zhong_volume.setDimR( dim_rad );
    zhong_volume.setDim( dim_zhong );
    zhong_volume.setVeclen( 1 );
    zhong_volume.calculateCoords();
    zhong_volume.readValues( filename_zhong );
    zhong_volume.updateMinMaxCoords();
    zhong_volume.updateMinMaxValues();

    // Update min/max values.
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

    // Update min/max coords.
    const kvs::Vec3 min_coord0 = yin_volume.minObjectCoord();
    const kvs::Vec3 min_coord1 = yang_volume.minObjectCoord();
    const kvs::Vec3 min_coord2 = zhong_volume.minObjectCoord();
    const kvs::Vec3 max_coord0 = yin_volume.maxObjectCoord();
    const kvs::Vec3 max_coord1 = yang_volume.maxObjectCoord();
    const kvs::Vec3 max_coord2 = zhong_volume.maxObjectCoord();
    const kvs::Real32 min_x = kvs::Math::Min( min_coord0.x(), min_coord1.x(), min_coord2.x() );
    const kvs::Real32 min_y = kvs::Math::Min( min_coord0.y(), min_coord1.y(), min_coord2.y() );
    const kvs::Real32 min_z = kvs::Math::Min( min_coord0.z(), min_coord1.z(), min_coord2.z() );
    const kvs::Real32 max_x = kvs::Math::Max( max_coord0.x(), max_coord1.x(), max_coord2.x() );
    const kvs::Real32 max_y = kvs::Math::Max( max_coord0.y(), max_coord1.y(), max_coord2.y() );
    const kvs::Real32 max_z = kvs::Math::Max( max_coord0.z(), max_coord1.z(), max_coord2.z() );
    //const kvs::Vec3 min_coord( min_x, min_y, min_z );
    //const kvs::Vec3 max_coord( max_x, max_y, max_z );
    const kvs::Vec3 min_coord( -1.0f, -1.0f, -1.0f );
    const kvs::Vec3 max_coord( 1.0f, 1.0f, 1.0f );
    yin_volume.setMinMaxObjectCoords( min_coord, max_coord );
    yin_volume.setMinMaxExternalCoords( min_coord, max_coord );
    yang_volume.setMinMaxObjectCoords( min_coord, max_coord );
    yang_volume.setMinMaxExternalCoords( min_coord, max_coord );
    zhong_volume.setMinMaxObjectCoords( min_coord, max_coord );
    zhong_volume.setMinMaxExternalCoords( min_coord, max_coord );

    // Extract isosurfaces
    const kvs::PolygonObject::NormalType n = kvs::PolygonObject::PolygonNormal;
    const bool d = true;
    const float isovalue1 = kvs::Math::Mix( min_value, max_value, 0.6f );
    kvs::PolygonObject* yin_iso1 = new YinYangVis::Isosurface( &yin_volume, isovalue1, n, d, tfunc );
    kvs::PolygonObject* yang_iso1 = new YinYangVis::Isosurface( &yang_volume, isovalue1, n, d, tfunc );
    kvs::PolygonObject* zhong_iso1 = new YinYangVis::Isosurface( &zhong_volume, isovalue1, n, d, tfunc );

    // Set up viewer application.
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setSize( 800, 600 );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    typedef kvs::glsl::PolygonRenderer Renderer;
    screen.registerObject( yin_iso1, new Renderer );
    screen.registerObject( yang_iso1, new Renderer );
    screen.registerObject( zhong_iso1, new Renderer );
    //screen.registerObject( yin_iso1, new kvs::Bounds );

    kvs::Axis3D* axis = new kvs::Axis3D();
    axis->setGridDrawMode( kvs::Axis3D::FarGrid );
    axis->setGridlinePattern( kvs::Axis3D::DashedLine );
    axis->setNumberOfGridlines( kvs::Vec3u( 4, 4, 4 ) );
    axis->setEnabledAntiAliasing( true );
    axis->setBackgroundColor( kvs::RGBColor( 220, 220, 220 ) );
    screen.registerObject( yin_iso1, axis );

    // Extract slice plane
    const kvs::Vec3 plane_point( 0, 0, 0 );
    if ( true )
    {
        const kvs::Vec3 plane_normal( 1, 0, 0 );
        kvs::PolygonObject* yin_slice = new YinYangVis::SlicePlane( &yin_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* yang_slice = new YinYangVis::SlicePlane( &yang_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* zhong_slice = new YinYangVis::SlicePlane( &zhong_volume, plane_point, plane_normal, tfunc );
        screen.registerObject( yin_slice, new Renderer );
        screen.registerObject( yang_slice, new Renderer );

        Renderer* renderer = new Renderer;
        renderer->setPolygonOffset( -0.001 );
        screen.registerObject( zhong_slice, renderer );
    }
    if ( true )
    {
        const kvs::Vec3 plane_normal( 0, 1, 0 );
        kvs::PolygonObject* yin_slice = new YinYangVis::SlicePlane( &yin_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* yang_slice = new YinYangVis::SlicePlane( &yang_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* zhong_slice = new YinYangVis::SlicePlane( &zhong_volume, plane_point, plane_normal, tfunc );
        screen.registerObject( yin_slice, new Renderer );
        screen.registerObject( yang_slice, new Renderer );

        Renderer* renderer = new Renderer;
        renderer->setPolygonOffset( -0.001 );
        screen.registerObject( zhong_slice, renderer );
    }
    if ( true )
    {
        const kvs::Vec3 plane_normal( 0, 0, 1 );
        kvs::PolygonObject* yin_slice = new YinYangVis::SlicePlane( &yin_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* yang_slice = new YinYangVis::SlicePlane( &yang_volume, plane_point, plane_normal, tfunc );
        kvs::PolygonObject* zhong_slice = new YinYangVis::SlicePlane( &zhong_volume, plane_point, plane_normal, tfunc );
        screen.registerObject( yin_slice, new Renderer );
        screen.registerObject( yang_slice, new Renderer );

        Renderer* renderer = new Renderer;
        renderer->setPolygonOffset( -0.001 );
        screen.registerObject( zhong_slice, renderer );
    }

    screen.show();

    kvs::Label label1( &screen );
    label1.setFont( kvs::Font( kvs::Font::Sans, kvs::Font::Bold, 25 ) );
    label1.setMargin( 10 );
    label1.setText( "MHD relaxation in a unit sphere" );
    label1.show();

    kvs::Label label2( &screen );
    label2.setFont( kvs::Font( kvs::Font::Sans, 20 ) );
    label2.setPosition( 0, 30 );
    label2.setMargin( 10 );
    label2.setText( "Jun28b.000.wyin.vx.n000550000.t00067" );
    label2.addText( "Jun28b.000.wyng.vx.n000550000.t00067" );
    label2.addText( "Jun28b.000.icore_3d.vx.n000550000.t00067" );
    label2.show();

    kvs::Label cmapbar_label( &screen );
    cmapbar_label.setFont( kvs::Font( kvs::Font::Sans, kvs::Font::Bold, 22 ) );
    cmapbar_label.setPosition( 610, 520 );
    cmapbar_label.setText( "Velocity (x)" );
    cmapbar_label.show();

    kvs::ColorMapBar cmap_bar( &screen );
    cmap_bar.setPosition( 600, 520 );
    cmap_bar.setColorMap( tfunc.colorMap() );
    cmap_bar.setCaption( " " );
    cmap_bar.setRange( -0.06, 0.06 );
    cmap_bar.setEnabledAntiAliasing( true );
    cmap_bar.show();

    kvs::Label slider_label( &screen );
    slider_label.setFont( kvs::Font( kvs::Font::Sans, kvs::Font::Bold, 22 ) );
    slider_label.setPosition( 610, 10 );
    slider_label.setText( "Isosurface" );
    slider_label.show();

    kvs::Slider slider( &screen );
    slider.setCaption( "" );
    slider.setPosition( 600, 0 );
    slider.setRange( -0.06, 0.06 );
    slider.setValue( 0.0 );
    slider.setWidth( 200 );
    slider.show();

    kvs::Light::SetModelTwoSide( true );

    return app.run();
}
