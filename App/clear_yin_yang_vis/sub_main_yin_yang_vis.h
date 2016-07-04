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
#include <kvs/LineObject>
#include <kvs/LineRenderer>
#include <kvs/StochasticLineRenderer>

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <Lib/YinYangVolumeObject.h>
#include <Lib/YinYangGridSampling.h>
#include <Lib/ZhongVolumeObject.h>
#include <Lib/ZhongGridSampling.h>

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

void ParticleBasedRenderingYinYang( kvs::glut::Screen& screen, YinYangVis::YinYangVolumeObject* volume, kvs::ColorMap cmap, kvs::OpacityMap omap, size_t repeats = 1 )
{

  std::cout << "repeats = " << repeats << std::endl;
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    const kvs::TransferFunction tfunc( cmap, omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new YinYangVis::YinYangGridSampling( volume, level, step, tfunc );
    timer.stop();
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;
    object->setName( volume->name() );
    kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
    object->multiplyXform( x );
    object->print( std::cout << std::endl );
      
    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

}

void ParticleBasedRenderingZhong( kvs::glut::Screen& screen, YinYangVis::ZhongVolumeObject* volume, kvs::ColorMap cmap, kvs::OpacityMap omap, size_t repeats = 1 )
{
  std::cout << "repeats = " << repeats << std::endl;
    const size_t subpixels = 1; // fixed to '1'
    const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
    const float step = 0.1f;

    const kvs::TransferFunction tfunc( cmap, omap );
    kvs::Timer timer( kvs::Timer::Start );
    kvs::PointObject* object = new YinYangVis::ZhongGridSampling( volume, level, step, tfunc );
    object->setName( volume->name() );
    kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
    object->multiplyXform( x );

    timer.stop();
    object->print( std::cout << std::endl );
    std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;

    kvs::glsl::ParticleBasedRenderer* renderer = new kvs::glsl::ParticleBasedRenderer();
    renderer->disableShading();

    screen.registerObject( object, renderer );

}


void SetMinMax( YinYangVis::YinYangVolumeObject* volume_yin, YinYangVis::YinYangVolumeObject* volume_yang, YinYangVis::ZhongVolumeObject* volume_zhong )
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

void SetMinMax( YinYangVis::YinYangVolumeObject* volume_yin, YinYangVis::YinYangVolumeObject* volume_yang )
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


void SetVolumeYin( YinYangVis::YinYangVolumeObject* volume, size_t rad_n, size_t lat_n, size_t lon_n, std::string filename )
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

void SetVolumeYang( YinYangVis::YinYangVolumeObject* volume, size_t rad_n, size_t lat_n, size_t lon_n, std::string filename )
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

void SetVolumeZhong( YinYangVis::ZhongVolumeObject* volume, size_t zhong_n, size_t rad_n, std::string filename )
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

kvs::LineObject* CreateLineObject( YinYangVis::YinYangVolumeObject* volume, size_t rad_n, size_t lat_n, size_t lon_n )
{
  const YinYangVis::YinYangVolumeObject::Range range_r = volume->rangeR();
  const YinYangVis::YinYangVolumeObject::Range range_theta = volume->rangeTheta();
  const YinYangVis::YinYangVolumeObject::Range range_phi = volume->rangePhi();

  const size_t nnodes = ( lat_n - 2 ) * 4 + ( lon_n - 4 - 2 ) * 4;
  kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
  kvs::Real32* pcoords = coords.data();

  for ( int k = 0; k < lon_n; k += lon_n - 5 )
  {
    const float phi = range_phi.min + range_phi.d * ( k - 2 );
    const float sin_phi = std::sin( phi );
    const float cos_phi = std::cos( phi );
    for ( int j = 0; j < rad_n; j += rad_n - 1 )
    {
      const float r = range_r.min + range_r.d * ( j - 1 );
      for ( int i = 0; i < lat_n - 2; i++ )
      {
        const float theta = range_theta.min + range_theta.d * i;
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );

        const float x = r * sin_theta * cos_phi;
        const float y = r * sin_theta * sin_phi;
        const float z = r * cos_theta;

        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x : -x;
        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y : z;
        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z : y;
      }
    }
  }

  for ( int k = 0; k < lat_n; k += lat_n - 3)
  {
    const float theta = range_theta.min + range_theta.d * ( k - 2 );
    const float sin_theta = std::sin( theta );
    const float cos_theta = std::cos( theta );
    for ( int j = 0; j < rad_n; j += rad_n - 1)
    {
      const float r = range_r.min + range_r.d * ( j - 1 );
      for ( int i = 1; i < lon_n - 5; i++ )
      {
        const float phi = range_phi.min + range_phi.d * i;
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );

        const float x = r * sin_theta * cos_phi;
        const float y = r * sin_theta * sin_phi;
        const float z = r * cos_theta;

        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x : -x;
        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y : z;
        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z : y;
      }
    }
  }

  const size_t nconnections = ( lat_n - 3 ) * 4 + ( lon_n - 5 ) * 4 + 4;
  kvs::ValueArray<kvs::UInt32> connections( nconnections * 2 );
  kvs::UInt32* pconnections = connections.data();

  size_t index;
  for ( int j = 0; j < lat_n * 3; j += lat_n - 2, index = j )
  {
    for ( int i = 0; i < lat_n - 3; i++ )
    {
      *(pconnections++) = i + j;
      *(pconnections++) = (i + 1) + j;
    }
  }

  size_t index2 = index;
  *(pconnections++) = 0;
  *(pconnections++) = index;
  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
  {
    *(pconnections++) = i;
    *(pconnections++) = i + 1;
  }
  *(pconnections++) = index;
  *(pconnections++) = 404;
  index++;
  *(pconnections++) = 202;
  *(pconnections++) = index;
  index2 = index;
  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
  {
    *(pconnections++) = i;
    *(pconnections++) = i + 1;
  }
  *(pconnections++) = index;
  *(pconnections++) = 606;
  index++;

  *(pconnections++) = 0 + 201;
  *(pconnections++) = index;
  index2 = index;
  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
  {
    *(pconnections++) = i;
    *(pconnections++) = i + 1;
  }
  *(pconnections++) = index;
  *(pconnections++) = 404 + 201;
  index++;
  *(pconnections++) = 202 + 201;
  *(pconnections++) = index;
  index2 = index;
  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
  {
    *(pconnections++) = i;
    *(pconnections++) = i + 1;
  }
  *(pconnections++) = index;
  *(pconnections++) = 605 + 201;

  *(pconnections++) = 0;
  *(pconnections++) = 202;
  *(pconnections++) = 201;
  *(pconnections++) = 403;
  *(pconnections++) = 404;
  *(pconnections++) = 606;
  *(pconnections++) = 605;
  *(pconnections++) = 807;

  kvs::ValueArray<kvs::UInt8> colors( nnodes * 3 );
  kvs::UInt8* pcolors = colors.data();
  for ( int i = 0; i < nnodes; i++ )
  {
    *(pcolors++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? 0 : 255;
    *(pcolors++) = 0;
    *(pcolors++) = 0;
  }

  kvs::LineObject* object = new kvs::LineObject();
  object->setCoords( coords );
  object->setColors( colors );
  object->setConnections( connections );
  object->setLineType( kvs::LineObject::Segment );
  object->setColorType( kvs::LineObject::VertexColor );

  kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
  object->multiplyXform( x );
  
  return object;
}

//kvs::LineObject* CreateLineObject( YinYangVis::ZhongVolumeObject* volume )
//{
  //const size_t dim_r = volume->dim();
  //const float r_min = volume->rangeR().min;
  //const float r_d = volume->rangeR().d;

//  const size_t nnodes = ( lat_n - 2 ) * 4 + ( lon_n - 4 - 2 ) * 4;
//  kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
//  kvs::Real32* pcoords = coords.data();
//
//  for ( int k = 0; k < lon_n; k += lon_n - 5 )
//  {
//    const float phi = range_phi.min + range_phi.d * (k-2);
//    const float sin_phi = std::sin( phi );
//    const float cos_phi = std::cos( phi );
//    for ( int j = 0; j < rad_n; j += rad_n - 1 )
//    {
//      const float r = range_r.min + range_r.d * (j-1);
//      for ( int i = 0; i < lat_n - 2; i++ )
//      {
//        const float theta = range_theta.min + range_theta.d * i;
//        const float sin_theta = std::sin( theta );
//        const float cos_theta = std::cos( theta );
//
//        const float x = r * sin_theta * cos_phi;
//        const float y = r * sin_theta * sin_phi;
//        const float z = r * cos_theta;
//
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x : -x;
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y : z;
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z : y;
//      }
//    }
//  }
//
//  for ( int k = 0; k < lat_n; k += lat_n - 3)
//  {
//    const float theta = range_theta.min + range_theta.d * k;
//    const float sin_theta = std::sin( theta );
//    const float cos_theta = std::cos( theta );
//    for ( int j = 0; j < rad_n; j += rad_n - 1)
//    {
//      const float r = range_r.min + range_r.d * j;
//      for ( int i = 1; i < lon_n - 5; i++ )
//      {
//        const float phi = range_phi.min + range_phi.d * i;
//        const float sin_phi = std::sin( phi );
//        const float cos_phi = std::cos( phi );
//
//        const float x = r * sin_theta * cos_phi;
//        const float y = r * sin_theta * sin_phi;
//        const float z = r * cos_theta;
//
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? x : -x;
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? y : z;
//        *(pcoords++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? z : y;
//      }
//    }
//  }
//
//  const size_t nconnections = ( lat_n - 3 ) * 4 + ( lon_n - 5 ) * 4 + 4;
//  kvs::ValueArray<kvs::UInt32> connections( nconnections * 2 );
//  kvs::UInt32* pconnections = connections.data();
//
//  size_t index;
//  for ( int j = 0; j < lat_n * 3; j += lat_n - 2, index = j )
//  {
//    for ( int i = 0; i < lat_n - 3; i++ )
//    {
//      *(pconnections++) = i + j;
//      *(pconnections++) = (i + 1) + j;
//    }
//  }
//
//  size_t index2 = index;
//  *(pconnections++) = 0;
//  *(pconnections++) = index;
//  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
//  {
//    *(pconnections++) = i;
//    *(pconnections++) = i + 1;
//  }
//  *(pconnections++) = index;
//  *(pconnections++) = 404;
//  index++;
//  *(pconnections++) = 202;
//  *(pconnections++) = index;
//  index2 = index;
//  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
//  {
//    *(pconnections++) = i;
//    *(pconnections++) = i + 1;
//  }
//  *(pconnections++) = index;
//  *(pconnections++) = 606;
//  index++;
//
//  *(pconnections++) = 0 + 201;
//  *(pconnections++) = index;
//  index2 = index;
//  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
//  {
//    *(pconnections++) = i;
//    *(pconnections++) = i + 1;
//  }
//  *(pconnections++) = index;
//  *(pconnections++) = 404 + 201;
//  index++;
//  *(pconnections++) = 202 + 201;
//  *(pconnections++) = index;
//  index2 = index;
//  for ( int i = index; i < index2 + lon_n - 7; i++, index = i )
//  {
//    *(pconnections++) = i;
//    *(pconnections++) = i + 1;
//  }
//  *(pconnections++) = index;
//  *(pconnections++) = 605 + 201;
//
//  *(pconnections++) = 0;
//  *(pconnections++) = 202;
//  *(pconnections++) = 201;
//  *(pconnections++) = 403;
//  *(pconnections++) = 404;
//  *(pconnections++) = 606;
//  *(pconnections++) = 605;
//  *(pconnections++) = 807;
//
//  kvs::ValueArray<kvs::UInt8> colors( nnodes * 3 );
//  kvs::UInt8* pcolors = colors.data();
//  for ( int i = 0; i < nnodes; i++ )
//  {
//    *(pcolors++) = 0;
//    *(pcolors++) = ( volume->gridType() == YinYangVis::YinYangVolumeObject::Yin ) ? 0 : 128;
//    *(pcolors++) = 0;
//  }
//
//  kvs::LineObject* object = new kvs::LineObject();
//  object->setCoords( coords );
//  object->setColors( colors );
//  object->setConnections( connections );
//  object->setLineType( kvs::LineObject::Segment );
//  object->setColorType( kvs::LineObject::VertexColor );
//
//  kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
//  object->multiplyXform( x );
//  
//  return object;
//}

#include <kvs/KeyPressEventListener>

class KeyPress : public kvs::KeyPressEventListener
{
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
      case kvs::Key::Four:
	{
	  if ( !scene()->hasObject("YinLine") ) break;
	  if ( scene()->object("YinLine")->isShown() )
	    {
	      std::cout << "Yin line hide" << std::endl;
	      scene()->object("YinLine")->hide();
	    }
	  else
	    {
	      std::cout << "Yin line show" << std::endl;
	      scene()->object("YinLine")->show();
	    }
	  break;
	}
      case kvs::Key::Five:
	{
	  if ( !scene()->hasObject("YangLine") ) break;
	  if ( scene()->object("YangLine")->isShown() )
	    {
	      std::cout << "Yang line hide" << std::endl;
	      scene()->object("YangLine")->hide();
	    }
	  else
	    {
	      std::cout << "Yang line show" << std::endl;
	      scene()->object("YangLine")->show();
	    }
	  break;
	}
      default: break;
      }
  }
};

int main_yin_yang_vis( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );
    kvs::glut::Screen screen( &app );
    screen.setBackgroundColor( kvs::RGBColor::White() );

    const size_t rad_n = atoi( argv[1] );
    const size_t lat_n = atoi( argv[2] );
    const size_t lon_n = atoi( argv[3] );
    const size_t zhong_n = atoi( argv[4] );
    
    /* const size_t rad_n = 201; */
    /* const size_t lat_n = 204; */
    /* const size_t lon_n = 608; */
    /* const size_t zhong_n = 222; */

    size_t repeats = atoi(argv[5]);
    
    const std::string filename_yin( argv[6] );
    YinYangVis::YinYangVolumeObject* volume_yin = new YinYangVis::YinYangVolumeObject();
    volume_yin->setName("Yin");
    SetVolumeYin( volume_yin, rad_n, lat_n, lon_n, filename_yin );

    const std::string filename_yang( argv[7] );
    YinYangVis::YinYangVolumeObject* volume_yang = new YinYangVis::YinYangVolumeObject();
    volume_yang->setName("Yang");
    SetVolumeYang( volume_yang, rad_n, lat_n, lon_n, filename_yang );
    
    const std::string filename_zhong( argv[8] );
    YinYangVis::ZhongVolumeObject* volume_zhong = new YinYangVis::ZhongVolumeObject();
    volume_zhong->setName("Zhong");
    SetVolumeZhong( volume_zhong,zhong_n, rad_n, filename_zhong );

    SetMinMax( volume_yin, volume_yang, volume_zhong );

    kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );

    kvs::OpacityMap omap( 256 );
    omap.addPoint( 0, 1.0 );
    omap.addPoint( 90, 0.0 );
    omap.addPoint( 180, 0.0 );
    omap.addPoint( 255, 1.0 );
    omap.create();

    kvs::LineObject* yinLineObject = CreateLineObject( volume_yin, rad_n, lat_n, lon_n );
    yinLineObject->setName("YinLine");
    kvs::StochasticLineRenderer* yinLineRenderer = new kvs::StochasticLineRenderer();
    screen.registerObject( yinLineObject, yinLineRenderer );

    kvs::LineObject* yangLineObject = CreateLineObject( volume_yang, rad_n, lat_n, lon_n );
    yangLineObject->setName("YangLine");
    kvs::StochasticLineRenderer* yangLineRenderer = new kvs::StochasticLineRenderer();
    screen.registerObject( yangLineObject, yangLineRenderer );

    ParticleBasedRenderingYinYang( screen, volume_yin, cmap, omap, repeats );
    delete volume_yin;
    ParticleBasedRenderingYinYang( screen, volume_yang, cmap, omap, repeats );
    delete volume_yang;
    ParticleBasedRenderingZhong( screen, volume_zhong, cmap, omap, repeats );
    delete volume_zhong;

    
    kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
    compositor->setRepetitionLevel( repeats );
    compositor->disableLODControl();
    screen.setEvent( compositor );

    KeyPress key_press;
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

    return app.run();
}
