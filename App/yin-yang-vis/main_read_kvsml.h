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
#include <kvs/Xform>
#include <kvs/XformControl>
#include <kvs/ScreenCaptureEvent>
#include <kvs/PaintEventListener>

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "YinYangVolumeObject.h"
#include "YinYangGridSampling.h"
#include "ZhongVolumeObject.h"
#include "ZhongGridSampling.h"

class FrameRate : public kvs::PaintEventListener
{
  kvs::StochasticRenderingCompositor* m_compositor;

public:
  FrameRate( kvs::StochasticRenderingCompositor* compositor ):
    m_compositor( compositor ) {}

  void update()
  {
    const int n = 50;

    static int counter = 1;
    static float accum = 0.0f;
    accum += m_compositor->timer().fps();
    if ( counter++ == n )
      {
	const float fps = accum / n;
	std::cout << fps << " [fps]" << std::endl;
	counter = 1;
	accum = 0.0f;
      }
  }
};

#include <kvs/KeyPressEventListener>

class KeyPress : public kvs::KeyPressEventListener
{
public:
  kvs::StochasticRenderingCompositor* m_compositor;
  KeyPress( kvs::StochasticRenderingCompositor* compositor )
  {
    m_compositor = compositor;
  }
  void update( kvs::KeyEvent* event )
  {
    switch ( event->key() )
      {
      case kvs::Key::One:
	{
	  if ( !scene()->hasObject("Yin") ) break;
	  if ( scene()->object("Yin")->isShown() )
	    {
	      std::cout << "Yin grid hide" << std::endl;
	      scene()->object("Yin")->hide();
	    }
	  else
	    {
	      std::cout << "Yin grid show" << std::endl;
	      scene()->object("Yin")->show();
	    }
	  break;
	}
      case kvs::Key::Two:
	{
	  if ( !scene()->hasObject("Yang") ) break;
	  if ( scene()->object("Yang")->isShown() )
	    {
	      std::cout << "Yang grid hide" << std::endl;
	      scene()->object("Yang")->hide();
	    }
	  else
	    {
	      std::cout << "Yang grid show" << std::endl;
	      scene()->object("Yang")->show();
	    }
	  break;
	}
      case kvs::Key::Three:
	{
	  if ( !scene()->hasObject("Zhong") ) break;
	  if ( scene()->object("Zhong")->isShown() )
	    {
	      std::cout << "Zhong grid hide" << std::endl;
	      scene()->object("Zhong")->hide();
	    }
	  else
	    {
	      std::cout << "Zhong grid show" << std::endl;
	      scene()->object("Zhong")->show();
	    }
	  break;
	}
      // case kvs::Key::l:
      // 	{
      // 	  if( m_compositor->isEnabledLODControl() )
      // 	    {
      // 	      m_compositor->disableLODControl();
      // 	    }
      // 	  else
      // 	    {
      // 	      m_compositor->enableLODControl();
      // 	    }
      // 	  screen()->redraw();
      // 	  break;
      // 	}
      default: break;
      }
  }
};

#include <kvs/glut/LegendBar>

class LegendBar : public kvs::glut::LegendBar
{
public:
  LegendBar( kvs::ScreenBase* screen ):
    kvs::glut::LegendBar( screen )
  {
    this->setWidth( 250 );
    this->setHeight( 70 );
  }

  void screenResized()
  {
    this->setX( this->screen()->width() - this->width() );
    this->setY( this->screen()->height() - this->height() );
  }
};

#include <kvs/glut/OrientationAxis>

class OrientationAxis : public kvs::glut::OrientationAxis
{
public:
  OrientationAxis( kvs::glut::Screen* screen ):
    kvs::glut::OrientationAxis( screen )
  {
    this->setMargin( 10 );
    this->setSize( 90 );
    this->setBoxType( kvs::glut::OrientationAxis::SolidBox );
    this->enableAntiAliasing();
  }

  void screenResized()
  {
    this->setY( this->screen()->height() - this->height() );
  }
};

int main_read_kvsml( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    kvs::Real32 min_value = -3.77136e-05;
    kvs::Real32 max_value = 2.29182e-05;
    
    size_t repeats = atoi( argv[1] );
    std::string dir_name( argv[2] );

    kvs::PointObject* object_yin = new kvs::PointObject;
    object_yin->read( dir_name  + "/yin.kvsml" );
    object_yin->setName("Yin");
    kvs::PointObject* object_yang = new kvs::PointObject;
    object_yang->read( dir_name + "/yang.kvsml" );
    object_yang->setName("Yang");
    kvs::PointObject* object_zhong = new kvs::PointObject;
    object_zhong->read( dir_name + "/zhong.kvsml" );
    object_zhong->setName("Zhong");
  
    kvs::glsl::ParticleBasedRenderer* renderer_yin = new kvs::glsl::ParticleBasedRenderer();
    renderer_yin->disableShading();
    kvs::glsl::ParticleBasedRenderer* renderer_yang = new kvs::glsl::ParticleBasedRenderer();
    renderer_yang->disableShading();
    kvs::glsl::ParticleBasedRenderer* renderer_zhong = new kvs::glsl::ParticleBasedRenderer();
    renderer_zhong->disableShading();

    screen.registerObject( object_yin, renderer_yin );
    screen.registerObject( object_yang, renderer_yang );
    screen.registerObject( object_zhong, renderer_zhong );
    
    kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
    compositor->setRepetitionLevel( repeats );
    if( atoi( argv[3] ) == 0 )
      {
	compositor->disableLODControl();
      }
    else
      {
	compositor->enableLODControl();
      }
    screen.setEvent( compositor );

    KeyPress key_press( compositor );
    screen.addEvent( &key_press );

    kvs::ScreenCaptureEvent capture;
    capture.setKey( kvs::Key::s );
    //capture.setFilename( "filename" );
    capture.setBasename( "ParticleImage" );
    screen.addEvent( &capture );
    
    kvs::TargetChangeEvent event;
    screen.addEvent( &event );
    screen.show();
    kvs::Light::SetModelTwoSide( true );

    FrameRate frame_rate( compositor );
    screen.addEvent( &frame_rate );
    
    kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 90, 0.0 );
    omap.addPoint( 180, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    kvs::TransferFunction tfunc( cmap, omap );
    LegendBar legend( &screen );
    legend.setCaption( "vx value" ); // 変数名を指定してください。
    legend.setColorMap( tfunc.colorMap() );
    legend.setRange( min_value, max_value ); // データ値の範囲を指定してください。
    legend.show();

    OrientationAxis orientation( &screen );
    orientation.show();

    return app.run();
}

