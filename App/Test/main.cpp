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
#include <kvs/Scene>
#include <kvs/ObjectManager>
#include <kvs/RendererManager>
#include <iostream>
#include <fstream>
#include "YinYangVolumeObject.h"


int main( int argc, char** argv )
{
  kvs::glut::Application app( argc, argv );

  const size_t rad_n = 201;
  const size_t lat_n = 204;
  const size_t lon_n = 608;
  const std::string filename_yin_value( "../bx_vx/oct09b.011.wyin.vx.n000250000.t00302" );
  local::YinYangVolumeObject* volume_yin = new local::YinYangVolumeObject();
//  volume_yin->setResolution( kvs::Vec3ui( rad_n, lat_n, lon_n ) );
  volume_yin->setDimR( rad_n );
  volume_yin->setDimTheta( lat_n );
  volume_yin->setDimPhi( lon_n );
  volume_yin->setVeclen( 1 );
  volume_yin->setGridTypeToYin();
  volume_yin->calculateCoords();
  volume_yin->readValues( filename_yin_value );

  kvs::UnstructuredVolumeObject* unsvolume_yin = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume_yin );
//  kvs::StructuredVolumeObject* strvolume_yin = local::YinYangVolumeObject::ToStructuredVolumeObject( volume_yin );

  kvs::PolygonObject* object_yin = new kvs::ExternalFaces( unsvolume_yin );
//  kvs::PolygonObject* object_yin = new kvs::ExternalFaces( strvolume_yin );

  kvs::glut::Screen screen( &app );
  screen.registerObject( object_yin );
  screen.show();

  return app.run();
}
