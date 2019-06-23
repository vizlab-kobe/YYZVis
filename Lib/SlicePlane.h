#include <kvs/PolygonObject>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/MapperBase>
#include <kvs/Module>

#include <kvs/UnstructuredVolumeObject>
#include <kvs/ValueArray>
#include <kvs/Vector3>
#include <kvs/Vector4>

namespace kvs
{
  class SlicePlane : public kvs::MapperBase, public kvs::PolygonObject
  {
    kvsModule( kvs::SlicePlane, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

  private:
    kvs::Vec4 coef;
    
  public:
    SlicePlane();
    SlicePlane(
	        const kvs::UnstructuredVolumeObject* volume,
		  const kvs::Vector3f point,
		const kvs::Vector3f normal);
    ~SlicePlane();

     SuperClass* exec( const kvs::ObjectBase* object );

     
    void setPlane( const kvs::Vec3 point, const kvs::Vec3 normal );

    void extract_plane( const kvs::ObjectBase* volume );
  
  kvs::Vector3f interpolate_vertex(  const kvs::Vector3f& vertex0,
    const kvs::Vector3f& vertex1 );

    int  calculate_table_index( const kvs::UnstructuredVolumeObject* volume ,size_t local_index[8] );

    float plane_function( kvs::Vec3 vertex );

    double interpolate_value(
        const kvs::UnstructuredVolumeObject* volume,
        const size_t index0,
        const size_t index1 );
  };
  
}

