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
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <YYZVis/Lib/YinYangGridSampling.h>
#include <YYZVis/Lib/ZhongGridSampling.h>
#include <kvs/PointObject>
#include <kvs/LineObject>
#include <kvs/glut/TransferFunctionEditor>
#include <kvs/OrientationAxis>
#include <kvs/StochasticRenderingCompositor>
#include <kvs/StochasticLineRenderer>
#include <kvs/ParticleBasedRenderer>


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

    kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );
    kvs::OpacityMap omap( 256 );
    omap.addPoint(   0, 1.0 );
    omap.addPoint( 150, 0.0 );
    omap.addPoint( 160, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();
    kvs::TransferFunction tfunc( cmap, omap );

    // Set up viewer application.
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setSize( 512, 512 );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    // Read Yin grid data.
    YYZVis::YinYangVolumeObject yin_volume;
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
    YYZVis::YinYangVolumeObject yang_volume;
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
    YYZVis::ZhongVolumeObject zhong_volume;
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

    // Generate particles.
    const size_t repeats = 100;
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;
    kvs::PointObject* yin_pnt = new YYZVis::YinYangGridSampling( &yin_volume, level, step, tfunc );
    kvs::PointObject* yang_pnt = new YYZVis::YinYangGridSampling( &yang_volume, level, step, tfunc );
    kvs::PointObject* zhong_pnt = new YYZVis::ZhongGridSampling( &zhong_volume, level, step, tfunc );

//    yin_pnt->setNormals( kvs::ValueArray<float>() );
//    yang_pnt->setNormals( kvs::ValueArray<float>() );
//    zhong_pnt->setNormals( kvs::ValueArray<float>() );

    typedef kvs::glsl::ParticleBasedRenderer Renderer;
    Renderer* yin_rend = new Renderer();
    Renderer* yang_rend = new Renderer();
    Renderer* zhong_rend = new Renderer();

    yin_rend->disableShading();
    yang_rend->disableShading();
    zhong_rend->disableShading();

//    screen.registerObject( yin_pnt );
//    screen.registerObject( yang_pnt );
//    screen.registerObject( zhong_pnt );
    screen.registerObject( yin_pnt, yin_rend );
    screen.registerObject( yang_pnt, yang_rend );
    screen.registerObject( zhong_pnt, zhong_rend );

    // Bounds
    kvs::Bounds bounds;
    kvs::LineObject* bounds_line = bounds.outputLineObject( yin_pnt );

//    screen.registerObject( bounds_line );
    screen.registerObject( bounds_line, new kvs::StochasticLineRenderer );

    // Set up stochastic rendering compositor.
    kvs::StochasticRenderingCompositor compositor( screen.scene() );
    compositor.setRepetitionLevel( repeats );
    compositor.enableLODControl();
    screen.setEvent( &compositor );

    screen.show();


    const kvs::Xform S = kvs::Xform::Scaling( kvs::Vec3::All( 0.8 ) );
//    const kvs::Xform R = kvs::Xform::Rotation( kvs::Mat3::RotationX( 15 ) * kvs::Mat3::RotationY( -25 ) );
    const kvs::Xform R = kvs::Xform::Rotation( kvs::Mat3::RotationX( 30 ) * kvs::Mat3::RotationY( -40 ) );
    const kvs::Xform X = S * R;
    yin_pnt->multiplyXform( X );
    yang_pnt->multiplyXform( X );
    zhong_pnt->multiplyXform( X );
    bounds_line->multiplyXform( X );

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
    cmapbar_label.setPosition( 310, 440 );
    cmapbar_label.setText( "Velocity (x)" );
    cmapbar_label.show();

    kvs::ColorMapBar cmap_bar( &screen );
    cmap_bar.setPosition( 300, 440 );
    cmap_bar.setColorMap( tfunc.colorMap() );
    cmap_bar.setCaption( " " );
    cmap_bar.setRange( -0.06, 0.06 );
    cmap_bar.setEnabledAntiAliasing( true );
    cmap_bar.show();

//    kvs::Light::SetModelTwoSide( true );

    kvs::glut::TransferFunctionEditor editor( &screen );
    editor.setVolumeObject( &yin_volume );
    editor.setTransferFunction( tfunc );
    editor.show();

//    kvs::OrientationAxis orientation_axis( &screen, screen.scene() );
    kvs::OrientationAxis orientation_axis( &screen, bounds_line );
//    orientation_axis.setBoxType( kvs::OrientationAxis::SolidBox );
    orientation_axis.setPosition( 0, 512 - 80 );
    orientation_axis.setAxisLineWidth( 2.0 );
    orientation_axis.setAxisLength( 4.0 );
    orientation_axis.setEnabledAntiAliasing( true );
    orientation_axis.show();

//    kvs::Light::SetModelTwoSide( true );

    return app.run();
}
