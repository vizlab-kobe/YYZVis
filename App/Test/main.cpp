#include <kvs/Message>
#include <kvs/PolygonObject>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/Endian>
#include <kvs/PointObject>
#include <kvs/CellByCellMetropolisSampling>
#include <kvs/ParticleBasedRenderer>
#include <kvs/TransferFunction>
#include <kvs/StochasticPolygonRenderer>
#include <kvs/StochasticTetrahedraRenderer>
#include <kvs/StochasticRenderingCompositor>
#include <kvs/ExternalFaces>
#include <kvs/Bounds>
#include <kvs/PolygonRenderer>
#include <kvs/Scene>
#include <kvs/ObjectManager>
#include <kvs/RendererManager>
#include <kvs/TargetChangeEvent>
#include <kvs/Timer>

#include <iostream>
#include <fstream>

#include "YinYangVolumeObject.h"
#include "YinYangGridSampling.h"
#include "ZhongVolumeObject.h"
#include "CellByCellMetropolisSampling.h"


void ExternalFaces( kvs::glut::Screen& screen, local::YinYangVolumeObject* volume )
{
    kvs::UnstructuredVolumeObject* temp = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume );
    temp->print( std::cout << std::endl );

    kvs::PolygonObject* object = new kvs::ExternalFaces( temp );
    object->print( std::cout << std::endl );
    delete temp;

    screen.registerObject( object, new kvs::PolygonRenderer() );
    screen.registerObject( object, new kvs::Bounds() );
}

void ParticleBasedRendering( kvs::glut::Screen& screen, local::YinYangVolumeObject* volume, size_t repeats = 1 )
{
    kvs::UnstructuredVolumeObject* temp = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume );
    temp->print( std::cout << std::endl );

    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 90, 0.0 );
    omap.addPoint( 180, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    const kvs::TransferFunction tfunc( omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new local::CellByCellMetropolisSampling( temp, level, step, tfunc );
    timer.stop();
    object->print( std::cout << std::endl );
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;
    delete temp;

    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

    kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
    compositor->setRepetitionLevel( repeats );
    compositor->enableLODControl();
    screen.setEvent( compositor );
}

void ParticleBasedRenderingYinYang( kvs::glut::Screen& screen, local::YinYangVolumeObject* volume, size_t repeats = 1 )
{
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 90, 0.0 );
    omap.addPoint( 180, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    const kvs::TransferFunction tfunc( omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new local::YinYangGridSampling( volume, level, step, tfunc );
    timer.stop();
    object->print( std::cout << std::endl );
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;

    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

    kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
    compositor->setRepetitionLevel( repeats );
    compositor->enableLODControl();
    screen.setEvent( compositor );
}

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );

    const size_t rad_n = 201;
    const size_t lat_n = 204;
    const size_t lon_n = 608;

    const std::string filename( argv[1] );
    local::YinYangVolumeObject* volume_yin = new local::YinYangVolumeObject();
    volume_yin->setDimR( rad_n );
    volume_yin->setDimTheta( lat_n );
    volume_yin->setDimPhi( lon_n );
    volume_yin->setVeclen( 1 );
    volume_yin->setGridTypeToYin();
    volume_yin->calculateCoords();
    volume_yin->readValues( filename );
    volume_yin->updateMinMaxCoords();
    volume_yin->updateMinMaxValues();
    volume_yin->print( std::cout << std::endl );

    //ExternalFaces( screen, volume_yin );

    size_t repeats = 16;
//    ParticleBasedRendering( screen, volume_yin, repeats );
    ParticleBasedRenderingYinYang( screen, volume_yin, repeats );
    delete volume_yin;

    kvs::TargetChangeEvent event;
    screen.addEvent( &event );

    screen.show();

    kvs::Light::SetModelTwoSide( true );

    return app.run();
}
