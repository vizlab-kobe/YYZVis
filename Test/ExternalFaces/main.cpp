#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/ColorMap>
#include <kvs/KeyPressEventListener>
#include <YYZVis/Lib/YinYangVolumeObject.h>
#include <YYZVis/Lib/ZhongVolumeObject.h>
#include <YYZVis/Lib/ExternalFaces.h>


class KeyPressEvent : public kvs::KeyPressEventListener
{
    void update( kvs::KeyEvent* event )
    {
        switch ( event->key() )
        {
        case kvs::Key::One:
        {
            auto* object = scene()->object( "Yin" );
            if ( object->isShown() ) { object->hide(); }
            else { object->show(); }
            break;
        }
        case kvs::Key::Two:
        {
            auto* object = scene()->object( "Yang" );
            if ( object->isShown() ) { object->hide(); }
            else { object->show(); }
            break;
        }
        case kvs::Key::Three:
        {
            auto* object = scene()->object( "Zhong" );
            if ( object->isShown() ) { object->hide(); }
            else { object->show(); }
            break;
        }
        default: break;
        }
    }
};

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setTitle( "YYZVis::ExternalFaces" );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    // Get input values.
    const size_t dim_rad = atoi( argv[1] );
    const size_t dim_lat = atoi( argv[2] );
    const size_t dim_lon = atoi( argv[3] );
    const size_t dim_zhong = atoi( argv[4] );
    const std::string yin_file( argv[5] );
    const std::string yang_file( argv[6] );
    const std::string zhong_file( argv[7] );

    // Import Yin data.
    YYZVis::YinYangVolumeObject yin_volume;
    yin_volume.setGridTypeToYin();
    yin_volume.setDimR( dim_rad );
    yin_volume.setDimTheta( dim_lat );
    yin_volume.setDimPhi( dim_lon );
    yin_volume.setVeclen( 1 );
    yin_volume.calculateCoords();
    yin_volume.readValues( yin_file );
    yin_volume.updateMinMaxCoords();
    yin_volume.updateMinMaxValues();

    // Import Yang data.
    YYZVis::YinYangVolumeObject yang_volume;
    yang_volume.setGridTypeToYang();
    yang_volume.setDimR( dim_rad );
    yang_volume.setDimTheta( dim_lat );
    yang_volume.setDimPhi( dim_lon );
    yang_volume.setVeclen( 1 );
    yang_volume.calculateCoords();
    yang_volume.readValues( yang_file );
    yang_volume.updateMinMaxCoords();
    yang_volume.updateMinMaxValues();

    // Import Zhong data.
    YYZVis::ZhongVolumeObject zhong_volume;
    zhong_volume.setDimR( dim_rad );
    zhong_volume.setDim( dim_zhong );
    zhong_volume.setVeclen( 1 );
    zhong_volume.calculateCoords();
    zhong_volume.readValues( zhong_file );
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

    // Dump.
    const kvs::Indent indent( 4 );
    yin_volume.print( std::cout << "YIN VOLUME DATA" << std::endl, indent );
    yang_volume.print( std::cout << "YANG VOLUME DATA" << std::endl, indent );
    zhong_volume.print( std::cout << "ZHONG VOLUME DATA" << std::endl, indent );

    // Extract faces.
    const kvs::ColorMap cmap = kvs::ColorMap::BrewerSpectral();
    kvs::PolygonObject* yin_object = new YYZVis::ExternalFaces( &yin_volume, cmap );
    kvs::PolygonObject* yang_object = new YYZVis::ExternalFaces( &yang_volume, cmap );
    kvs::PolygonObject* zhong_object = new YYZVis::ExternalFaces( &zhong_volume, cmap );

    yin_object->setName( "Yin" );
    yang_object->setName( "Yang" );
    zhong_object->setName( "Zhong" );

    screen.registerObject( yin_object );
    screen.registerObject( yang_object );
    screen.registerObject( zhong_object );

    // Key press event.
    KeyPressEvent key_event;
    screen.addEvent( &key_event );

    screen.show();

    return app.run();
}
