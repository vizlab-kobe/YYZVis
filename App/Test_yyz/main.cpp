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
#include <kvs/Vector3>
#include <kvs/StochasticPolygonRenderer>
#include <kvs/StochasticTetrahedraRenderer>
#include <kvs/StochasticRenderingCompositor>
#include <kvs/ExternalFaces>
#include <kvs/ExtractEdges>
#include <kvs/Bounds>
#include <kvs/PolygonRenderer>
#include <kvs/Scene>
#include <kvs/ObjectManager>
#include <kvs/RendererManager>

#include <iostream>
#include <fstream>

#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"


int main( int argc, char** argv )
{
  kvs::glut::Application app( argc, argv );

  const size_t c_nr = 201;
  const size_t c_nt_half = 201;

  const std::string filename_yin( "../../../bx_vx/oct09b.011.wyin.vx.n000250000.t00302" );
  const std::string filename_yang( "../../../bx_vx/oct09b.011.wyng.vx.n000250000.t00302" );
  const std::string filename_zhong( "../../../bx_vx/oct09b.011.icore_3d.vx.n000250000.t00302" );
  //const std::string filename_yin( argv[1] );
  //const std::string filename_yang( argv[2] );
  //const std::string filename_zhong( argv[3] );

  //---Yin Volume-----
  local::YinYangVolumeObject* volume_yin = new local::YinYangVolumeObject();
  volume_yin->setDimR( c_nr );  //dim_r = c_nr, range_min = 0.35f, range_max = 1.0f
  volume_yin->setDimTheta( c_nt_half );  //dim_theta = c_nt_half + 3, overwrap = 2
  volume_yin->setDimPhi( c_nt_half );  //dim_phi = c_nt_half * 3 + 5, overwrap = 4
  volume_yin->setVeclen( 1 );
  volume_yin->setGridTypeToYin();
  volume_yin->readValues( filename_yin );
  volume_yin->updateMinMaxValues();

  //---Yang Volume-----
  local::YinYangVolumeObject* volume_yang = new local::YinYangVolumeObject();
  volume_yang->setDimR( c_nr );  //dim_r = c_nr, range_min = 0.35f, range_max = 1.0f
  volume_yang->setDimTheta( c_nt_half );  //dim_theta = c_nt_half + 3, overwrap = 2
  volume_yang->setDimPhi( c_nt_half );  //dim_phi   = c_nt_half * 3 + 5, overwrap = 4
  volume_yang->setVeclen( 1 );
  volume_yang->setGridTypeToYang();
  volume_yang->readValues( filename_yang );
  volume_yang->updateMinMaxValues();

  //---Zhong Volume-----
  local::ZhongVolumeObject* volume_zhong = new local::ZhongVolumeObject();
  volume_zhong->setDimR( c_nr );  //dim_r = c_nr, range_min = 0.35f, range_max = 1.0f
  volume_zhong->setDim( volume_zhong );
  volume_zhong->setVeclen( 1 );
  volume_zhong->readValues( filename_zhong );
  volume_zhong->updateMinMaxValues();
 
  //---To Unstructured Volume Object-----
  kvs::UnstructuredVolumeObject* unsvolume_yin = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume_yin );
  //kvs::StructuredVolumeObject* strvolume_yin = local::YinYangVolumeObject::ToStructuredVolumeObject( volume_yin );
  unsvolume_yin->print( std::cout << std::endl );
  delete volume_yin;

  kvs::UnstructuredVolumeObject* unsvolume_yang = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume_yang );
  //kvs::StructuredVolumeObject* strvolume_yang = local::YinYangVolumeObject::ToStructuredVolumeObject( volume_yang );
  unsvolume_yang->print( std::cout << std::endl );
  delete volume_yang;

  kvs::UnstructuredVolumeObject* unsvolume_zhong = local::ZhongVolumeObject::ToUnstructuredVolumeObject( volume_zhong );
  //kvs::StructuredVolumeObject* strvolume_zhong = local::YinYangVolumeObject::ToStructuredVolumeObject( volume_zhong );
  unsvolume_zhong->print( std::cout << std::endl );
  delete volume_zhong;

  //---ParticleBasedRenderer Parameter------
  const size_t repetitions = 80;
  const size_t subpixels = 1;
  const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repetitions ) ) );
  const float step = 0.1f;

  kvs::OpacityMap omap( 256 );
  omap.addPoint( 0, 1.0 );
  omap.addPoint( 90, 0.0 );
  omap.addPoint( 180, 0.0 );
  omap.addPoint( 255, 1.0 );
  omap.create();

  const kvs::TransferFunction tfunc( omap );

  //---Set Min Max-----
  const float min_yyz_value  = kvs::Math::Min( unsvolume_yin->minValue(), unsvolume_yang->minValue() , unsvolume_zhong->minValue() );
  const float max_yyz_value  = kvs::Math::Max( unsvolume_yin->maxValue(), unsvolume_yang->maxValue() , unsvolume_zhong->maxValue() );

  unsvolume_yin->setMinMaxValues( min_yyz_value, max_yyz_value );
  unsvolume_yang->setMinMaxValues( min_yyz_value, max_yyz_value );
  unsvolume_zhong->setMinMaxValues( min_yyz_value, max_yyz_value );
  
  const kvs::Vec3 min_yin_coord = unsvolume_yin->minObjectCoord();
  const kvs::Vec3 min_yang_coord = unsvolume_yang->minObjectCoord();
  const kvs::Vec3 min_zhong_coord = unsvolume_zhong->minObjectCoord();
  const kvs::Vec3 max_yin_coord = unsvolume_yin->maxObjectCoord();
  const kvs::Vec3 max_yang_coord = unsvolume_yang->maxObjectCoord();
  const kvs::Vec3 max_zhong_coord = unsvolume_zhong->maxObjectCoord();

  const float min_x_coord = kvs::Math::Min( min_yin_coord.x(), min_yang_coord.x() , min_zhong_coord.x() );
  const float min_y_coord = kvs::Math::Min( min_yin_coord.y(), min_yang_coord.y() , min_zhong_coord.y() );
  const float min_z_coord = kvs::Math::Min( min_yin_coord.z(), min_yang_coord.z() , min_zhong_coord.z() );
  const kvs::Vec3 min_coord = kvs::Vec3( min_x_coord, min_y_coord, min_z_coord );

  const float max_x_coord = kvs::Math::Max( max_yin_coord.x(), max_yang_coord.x() , max_zhong_coord.x() );
  const float max_y_coord = kvs::Math::Max( max_yin_coord.y(), max_yang_coord.y() , max_zhong_coord.y() );
  const float max_z_coord = kvs::Math::Max( max_yin_coord.z(), max_yang_coord.z() , max_zhong_coord.z() );
  const kvs::Vec3 max_coord = kvs::Vec3( max_x_coord, max_y_coord, max_z_coord );

  unsvolume_yin->setMinMaxObjectCoords( min_coord, max_coord );
  unsvolume_yang->setMinMaxObjectCoords( min_coord, max_coord );
  unsvolume_zhong->setMinMaxObjectCoords( min_coord, max_coord );
  
  unsvolume_yin->setMinMaxExternalCoords( min_coord, max_coord );
  unsvolume_yang->setMinMaxExternalCoords( min_coord, max_coord );
  unsvolume_zhong->setMinMaxExternalCoords( min_coord, max_coord );

  //--To Perticle Based Volume Object-----
  kvs::PointObject* object_yin = new kvs::CellByCellMetropolisSampling( unsvolume_yin, level, step, tfunc );
  //kvs::PolygonObject* object_yin = new kvs::ExternalFaces( unsvolume_yin );
  object_yin->print( std::cout << std::endl );
  delete unsvolume_yin;

  kvs::PointObject* object_yang = new kvs::CellByCellMetropolisSampling( unsvolume_yang, level, step, tfunc );
  //kvs::PolygonObject* object_yang = new kvs::ExternalFaces( unsvolume_yang );
  object_yang->print( std::cout << std::endl );
  delete unsvolume_yang;
  
  kvs::PointObject* object_zhong = new kvs::CellByCellMetropolisSampling( unsvolume_zhong, level, step, tfunc );
  //kvs::PolygonObject* object_zhong = new kvs::ExternalFaces( unsvolume_zhong );
  object_zhong->print( std::cout << std::endl );
  delete unsvolume_zhong;

  //---Renderer-----
  kvs::glsl::ParticleBasedRenderer* renderer_yin = new kvs::glsl::ParticleBasedRenderer();
  renderer_yin->enableShading();

  kvs::glsl::ParticleBasedRenderer* renderer_yang = new kvs::glsl::ParticleBasedRenderer();
  renderer_yang->enableShading();

  kvs::glsl::ParticleBasedRenderer* renderer_zhong = new kvs::glsl::ParticleBasedRenderer();
  renderer_zhong->enableShading();

  kvs::glut::Screen screen( &app );
  screen.registerObject( object_yin, renderer_yin );
  screen.registerObject( object_yang, renderer_yang );
  screen.registerObject( object_zhong, renderer_zhong );
  screen.show();

  //---Compositor-----
    kvs::StochasticRenderingCompositor compositor( screen.scene() );
  compositor.setRepetitionLevel( repetitions);
  compositor.enableLODControl();
  screen.setEvent( &compositor );

  return app.run();
}
