#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/ColorMap>
#include <kvs/KeyPressEventListener>
#include <YYZVis/Lib/YinVolumeImporter.h>
#include <YYZVis/Lib/YangVolumeImporter.h>
#include <YYZVis/Lib/ZhongVolumeImporter.h>
#include <YYZVis/Lib/SlicePlane.h>


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
    screen.setTitle( "YYZVis::SlicePlane" );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    // Import YYZ data.
    const std::string input_file( argv[1] );
    auto yin_volume = YYZVis::YinVolumeImporter( input_file );
    auto yng_volume = YYZVis::YangVolumeImporter( input_file );
    auto zng_volume = YYZVis::ZhongVolumeImporter( input_file );

    // Update min/max values.
    const kvs::Real32 min_value0 = yin_volume.minValue();
    const kvs::Real32 min_value1 = yng_volume.minValue();
    const kvs::Real32 min_value2 = zng_volume.minValue();
    const kvs::Real32 max_value0 = yin_volume.maxValue();
    const kvs::Real32 max_value1 = yng_volume.maxValue();
    const kvs::Real32 max_value2 = zng_volume.maxValue();
    const kvs::Real32 min_value = kvs::Math::Min( min_value0, min_value1, min_value2 );
    const kvs::Real32 max_value = kvs::Math::Max( max_value0, max_value1, max_value2 );
    yin_volume.setMinMaxValues( min_value, max_value );
    yng_volume.setMinMaxValues( min_value, max_value );
    zng_volume.setMinMaxValues( min_value, max_value );

    // Dump.
    const kvs::Indent indent( 4 );
    yin_volume.print( std::cout << "YIN VOLUME DATA" << std::endl, indent );
    yng_volume.print( std::cout << "YANG VOLUME DATA" << std::endl, indent );
    zng_volume.print( std::cout << "ZHONG VOLUME DATA" << std::endl, indent );

    // Extract slice planes.
    const kvs::Vec3 point( 0.0f, 0.0f, 0.0f );
    const kvs::Vec3 normal( 0.0f, 0.0f, 1.0f );
    const kvs::ColorMap cmap = kvs::ColorMap::BrewerSpectral();
    kvs::PolygonObject* yin_object = new YYZVis::SlicePlane( &yin_volume, point, normal, cmap );
    kvs::PolygonObject* yng_object = new YYZVis::SlicePlane( &yng_volume, point, normal, cmap );
    kvs::PolygonObject* zng_object = new YYZVis::SlicePlane( &zng_volume, point, normal, cmap );

    yin_object->setName( "Yin" );
    yng_object->setName( "Yang" );
    zng_object->setName( "Zhong" );

    screen.registerObject( yin_object );
    screen.registerObject( yng_object );
    screen.registerObject( zng_object );

    // Key press event.
    KeyPressEvent key_event;
    screen.addEvent( &key_event );

    screen.show();

    return app.run();
}
