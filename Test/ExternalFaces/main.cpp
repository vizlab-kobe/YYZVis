#include <kvs/Application>
#include <kvs/Screen>
#include <kvs/ColorMap>
#include <kvs/KeyPressEventListener>
#include <YYZVis/Lib/YinVolumeImporter.h>
#include <YYZVis/Lib/YangVolumeImporter.h>
#include <YYZVis/Lib/ZhongVolumeImporter.h>
#include <YYZVis/Lib/UpdateMinMaxValues.h>
#include <YYZVis/Lib/ExternalFaces.h>


int main( int argc, char** argv )
{
    kvs::Application app( argc, argv );
    kvs::Screen screen( &app );
    screen.setTitle( "YYZVis::ExternalFaces" );
    screen.setBackgroundColor( kvs::RGBColor::White() );
    screen.create();

    // Import YYZ data.
    const std::string input_file( argv[1] );
    auto* yin_volume = new YYZVis::YinVolumeImporter( input_file );
    auto* yng_volume = new YYZVis::YangVolumeImporter( input_file );
    auto* zng_volume = new YYZVis::ZhongVolumeImporter( input_file );
    YYZVis::UpdateMinMaxValues( yin_volume, yng_volume, zng_volume );

    // Dump.
    const kvs::Indent indent( 4 );
    yin_volume->print( std::cout << "YIN VOLUME DATA" << std::endl, indent );
    yng_volume->print( std::cout << "YANG VOLUME DATA" << std::endl, indent );
    zng_volume->print( std::cout << "ZHONG VOLUME DATA" << std::endl, indent );

    // Extract faces.
    const kvs::ColorMap cmap = kvs::ColorMap::BrewerSpectral();
    auto* yin_object = new YYZVis::ExternalFaces( yin_volume, cmap );
    auto* yng_object = new YYZVis::ExternalFaces( yng_volume, cmap );
    auto* zng_object = new YYZVis::ExternalFaces( zng_volume, cmap );
    delete yin_volume;
    delete yng_volume;
    delete zng_volume;

    screen.registerObject( yin_object );
    screen.registerObject( yng_object );
    screen.registerObject( zng_object );

    // Key press event.
    kvs::KeyPressEventListener key_event( [&] ( kvs::KeyEvent* e )
    {
        auto onoff = [] ( kvs::ObjectBase* o )
        {
            if ( o->isVisible() ) { o->hide(); }
            else { o->show(); }
        };

        switch ( e->key() )
        {
        case kvs::Key::One: { onoff( yin_object ); break; }
        case kvs::Key::Two: { onoff( yng_object ); break; }
        case kvs::Key::Three: { onoff( zng_object ); break; }
        default: break;
        }

        screen.redraw();
    } );
    screen.addEvent( &key_event );

    return app.run();
}
