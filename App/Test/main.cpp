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

  const size_t rad_n = 201;
  const size_t lat_n = 204;
  const size_t lon_n = 608;
//  const std::string filename( "../bx_vx/oct09b.011.wyin.vx.n000250000.t00302" );
  const std::string filename( argv[1] );
  local::YinYangVolumeObject* volume1 = new local::YinYangVolumeObject();
  volume1->setDimR( rad_n );
  volume1->setDimTheta( lat_n );
  volume1->setDimPhi( lon_n );
  volume1->setVeclen( 1 );
  volume1->setGridTypeToYin();
  volume1->readValues( filename );

  kvs::UnstructuredVolumeObject* volume2 = local::YinYangVolumeObject::ToUnstructuredVolumeObject( volume1 );
  //kvs::StructuredVolumeObject* volume2 = local::YinYangVolumeObject::ToStructuredVolumeObject( volume1 );
  volume2->print( std::cout << std::endl );
  delete volume1;

  kvs::PolygonObject* object = new kvs::ExternalFaces( volume2 );
  object->print( std::cout << std::endl );
  delete volume2;

  kvs::glut::Screen screen( &app );
  screen.registerObject( object, new kvs::PolygonRenderer() );
  screen.registerObject( object, new kvs::Bounds() );
  screen.show();

  return app.run();
}
