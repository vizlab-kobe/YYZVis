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
#include <kvs/DivergingColorMap>

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "YinYangVolumeObject.h"
#include "YinYangGridSampling.h"
#include "ZhongVolumeObject.h"
#include "ZhongGridSampling.h"
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
        
    const kvs::TransferFunction tfunc( 256 );
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

  std::cout << "repeats = " << repeats << std::endl;
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 140, 0.0 );
    omap.addPoint( 160, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    const kvs::TransferFunction tfunc( cmap, omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new local::YinYangGridSampling( volume, level, step, tfunc );
    object->setName( volume->name() );
    timer.stop();
    object->print( std::cout << std::endl );
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;

    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

}

void ParticleBasedRenderingZhong( kvs::glut::Screen& screen, local::ZhongVolumeObject* volume, size_t repeats = 1 )
{

  std::cout << "repeats = " << repeats << std::endl;
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 140, 0.0 );
    omap.addPoint( 160, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    const kvs::TransferFunction tfunc( cmap, omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new local::ZhongGridSampling( volume, level, step, tfunc );
    object->setName( volume->name() );
    timer.stop();
    object->print( std::cout << std::endl );
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;

    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

}


void SetMinMax( local::YinYangVolumeObject* volume_yin, local::YinYangVolumeObject* volume_yang, local::ZhongVolumeObject* volume_zhong )
{
    const float min_yyz_value  = kvs::Math::Min( volume_yin->minValue(), volume_yang->minValue(), volume_zhong->minValue() );
    const float max_yyz_value  = kvs::Math::Max( volume_yin->maxValue(), volume_yang->maxValue(), volume_zhong->maxValue() );

    volume_yin->setMinMaxValues( min_yyz_value, max_yyz_value );
    volume_yang->setMinMaxValues( min_yyz_value, max_yyz_value );
    volume_zhong->setMinMaxValues( min_yyz_value, max_yyz_value );

    const kvs::Vec3 min_yin_coord = volume_yin->minObjectCoord();
    const kvs::Vec3 min_yang_coord = volume_yang->minObjectCoord();
    const kvs::Vec3 min_zhong_coord = volume_zhong->minObjectCoord();
    const kvs::Vec3 max_yin_coord = volume_yin->maxObjectCoord();
    const kvs::Vec3 max_yang_coord = volume_yang->maxObjectCoord();
    const kvs::Vec3 max_zhong_coord = volume_zhong->maxObjectCoord();

    const float min_x_coord = kvs::Math::Min( min_yin_coord.x(), min_yang_coord.x(), min_zhong_coord.x() );
    const float min_y_coord = kvs::Math::Min( min_yin_coord.y(), min_yang_coord.y(), min_zhong_coord.y() );
    const float min_z_coord = kvs::Math::Min( min_yin_coord.z(), min_yang_coord.z(), min_zhong_coord.z() );
    const kvs::Vec3 min_coord = kvs::Vec3( min_x_coord, min_y_coord, min_z_coord );

    const float max_x_coord = kvs::Math::Max( max_yin_coord.x(), max_yang_coord.x(), max_zhong_coord.x() );
    const float max_y_coord = kvs::Math::Max( max_yin_coord.y(), max_yang_coord.y(), max_zhong_coord.y() );
    const float max_z_coord = kvs::Math::Max( max_yin_coord.z(), max_yang_coord.z(), max_zhong_coord.z() );
    const kvs::Vec3 max_coord = kvs::Vec3( max_x_coord, max_y_coord, max_z_coord );

    volume_yin->setMinMaxObjectCoords( min_coord, max_coord );
    volume_yang->setMinMaxObjectCoords( min_coord, max_coord );
    volume_zhong->setMinMaxObjectCoords( min_coord, max_coord );

    volume_yin->setMinMaxExternalCoords( min_coord, max_coord );
    volume_yang->setMinMaxExternalCoords( min_coord, max_coord );
    volume_zhong->setMinMaxExternalCoords( min_coord, max_coord );
}

void SetMinMax( local::YinYangVolumeObject* volume_yin, local::YinYangVolumeObject* volume_yang )
{
    const float min_yyz_value  = kvs::Math::Min( volume_yin->minValue(), volume_yang->minValue() );
    const float max_yyz_value  = kvs::Math::Max( volume_yin->maxValue(), volume_yang->maxValue() );

    volume_yin->setMinMaxValues( min_yyz_value, max_yyz_value );
    volume_yang->setMinMaxValues( min_yyz_value, max_yyz_value );

    const kvs::Vec3 min_yin_coord = volume_yin->minObjectCoord();
    const kvs::Vec3 min_yang_coord = volume_yang->minObjectCoord();
    const kvs::Vec3 max_yin_coord = volume_yin->maxObjectCoord();
    const kvs::Vec3 max_yang_coord = volume_yang->maxObjectCoord();

    const float min_x_coord = kvs::Math::Min( min_yin_coord.x(), min_yang_coord.x() );
    const float min_y_coord = kvs::Math::Min( min_yin_coord.y(), min_yang_coord.y() );
    const float min_z_coord = kvs::Math::Min( min_yin_coord.z(), min_yang_coord.z() );
    const kvs::Vec3 min_coord = kvs::Vec3( min_x_coord, min_y_coord, min_z_coord );

    const float max_x_coord = kvs::Math::Max( max_yin_coord.x(), max_yang_coord.x() );
    const float max_y_coord = kvs::Math::Max( max_yin_coord.y(), max_yang_coord.y() );
    const float max_z_coord = kvs::Math::Max( max_yin_coord.z(), max_yang_coord.z() );
    const kvs::Vec3 max_coord = kvs::Vec3( max_x_coord, max_y_coord, max_z_coord );
    
    volume_yin->setMinMaxObjectCoords( min_coord, max_coord );
    volume_yang->setMinMaxObjectCoords( min_coord, max_coord );
    
    volume_yin->setMinMaxExternalCoords( min_coord, max_coord );
    volume_yang->setMinMaxExternalCoords( min_coord, max_coord );
    }


void SetVolumeYin( local::YinYangVolumeObject* volume, size_t rad_n, size_t lat_n, size_t lon_n, std::string filename )
{
    volume->setDimR( rad_n );
    volume->setDimTheta( lat_n );
    volume->setDimPhi( lon_n );
    volume->setVeclen( 1 );
    volume->setGridTypeToYin();
    volume->calculateCoords();
    volume->readValues( filename );
    volume->updateMinMaxCoords();
    volume->updateMinMaxValues();
    volume->print( std::cout << std::endl );
}

void SetVolumeYang( local::YinYangVolumeObject* volume, size_t rad_n, size_t lat_n, size_t lon_n, std::string filename )
{
    volume->setDimR( rad_n );
    volume->setDimTheta( lat_n );
    volume->setDimPhi( lon_n );
    volume->setVeclen( 1 );
    volume->setGridTypeToYang();
    volume->calculateCoords();
    volume->readValues( filename );
    volume->updateMinMaxCoords();
    volume->updateMinMaxValues();
    volume->print( std::cout << std::endl );
}

void SetVolumeZhong( local::ZhongVolumeObject* volume, size_t zhong_n, size_t rad_n, std::string filename )
{
  volume->setDimR( rad_n );
    volume->setDim( zhong_n );
    volume->setVeclen( 1 );
    volume->calculateCoords();
    volume->readValues( filename );
    volume->updateMinMaxCoords();
    volume->updateMinMaxValues();
    volume->print( std::cout << std::endl );
}

#include <kvs/KeyPressEventListener>
#include <kvs/Xform>
#include <kvs/XformControl>

class KeyPress : public kvs::KeyPressEventListener
{
private:
  kvs::Xform x;
  
  void update( kvs::KeyEvent* event )
  {
    switch ( event->key() )
      {
      case kvs::Key::One:
	{
	  std::cout << "aaa" << std::endl;
	  if ( !scene()->hasObject("Yin") ) break;
	  if ( scene()->object("Yin")->isShown() ) scene()->object("Yin")->hide();
	  else scene()->object("Yin")->show();
	  break;
	}
      case kvs::Key::Two:
	{
	  std::cout << "bbb" << std::endl;
	  if ( !scene()->hasObject("Yang") ) break;
	  if ( scene()->object("Yang")->isShown() ) scene()->object("Yang")->hide();
	  else scene()->object("Yang")->show();
	  break;
	}
      case kvs::Key::Three:
	{
	  std::cout << "ccc" << std::endl;
	  if ( !scene()->hasObject("Zhong") ) break;
	  if ( scene()->object("Zhong")->isShown() ) scene()->object("Zhong")->hide();
	  else scene()->object("Zhong")->show();
	  break;
	}
      case kvs::Key::Four:
	{
	  x = kvs::Xform::Rotation( kvs::Mat3::RotationX(45) );
	  scene()->object("Yin")->multiplyXform( x );
	  scene()->object("Yang")->multiplyXform( x );
	  scene()->object("Zhogn")->multiplyXform( x );
	}
      case kvs::Key::Five:
	{
	  x = kvs::Xform::Rotation( kvs::Mat3::RotationY(45) );
	  scene()->object("Yin")->multiplyXform( x );
	  scene()->object("Yang")->multiplyXform( x );
	  scene()->object("Zhogn")->multiplyXform( x );
	}
      case kvs::Key::Six:
	{
	  x = kvs::Xform::Rotation( kvs::Mat3::RotationZ(45) );
	  scene()->object("Yin")->multiplyXform( x );
	  scene()->object("Yang")->multiplyXform( x );
	  scene()->object("Zhogn")->multiplyXform( x );
	}
      default: break;
      }
  }
};

int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    const size_t rad_n = 201;
    const size_t lat_n = 204;
    const size_t lon_n = 608;
    const size_t zhong_n = 222;

    //const std::string filename_yang( argv[2] );
    const std::string filename_yang( "../../../bx_vx/oct09b.011.wyng.vx.n000250000.t00302" );
    local::YinYangVolumeObject* volume_yang = new local::YinYangVolumeObject();
    volume_yang->setName("Yang");
    SetVolumeYang( volume_yang, rad_n, lat_n, lon_n, filename_yang );

    //const std::string filename_yin( argv[1] );
    const std::string filename_yin( "../../../bx_vx/oct09b.011.wyin.vx.n000250000.t00302" );
    local::YinYangVolumeObject* volume_yin = new local::YinYangVolumeObject();
    volume_yin->setName("Yin");
    SetVolumeYin( volume_yin, rad_n, lat_n, lon_n, filename_yin );

    //const std::string filename_zhong( argv[3] );
    const std::string filename_zhong( "../../../bx_vx/oct09b.011.icore_3d.vx.n000250000.t00302" );
    local::ZhongVolumeObject* volume_zhong = new local::ZhongVolumeObject();
    volume_zhong->setName("Zhong");
    SetVolumeZhong( volume_zhong,zhong_n, rad_n, filename_zhong );
    
    SetMinMax( volume_yin, volume_yang, volume_zhong );

    size_t repeats = atoi(argv[1]);

    //ParticleBasedRendering( screen, volume_yang, repeats );
    ParticleBasedRenderingYinYang( screen, volume_yang, repeats );
    delete volume_yang;
    ParticleBasedRenderingYinYang( screen, volume_yin, repeats );
    delete volume_yin;
    ParticleBasedRenderingZhong( screen, volume_zhong, repeats );
    delete volume_zhong;
    
    kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
    compositor->setRepetitionLevel( repeats );
    compositor->enableLODControl();
    screen.setEvent( compositor );

    KeyPress key_press;
    screen.addEvent( &key_press );
    
    kvs::TargetChangeEvent event;
    screen.addEvent( &event );
    screen.show();
    kvs::Light::SetModelTwoSide( true );

    return app.run();
}
