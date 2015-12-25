#pragma once
#include <kvs/Module>
#include <kvs/VolumeObjectBase>
//#include <kvs/Vector3>
#include <kvs/StructuredVolumeObject>
#include <kvs/UnstructuredVolumeObject>


namespace local
{

class ZhongVolumeObject : public kvs::VolumeObjectBase
{
    kvsModule( local::ZhongVolumeObject, Object );
    kvsModuleBaseClass( kvs::VolumeObjectBase );

public:
    static kvs::StructuredVolumeObject* ToStructuredVolumeObject( const local::ZhongVolumeObject* object );
    static kvs::UnstructuredVolumeObject* ToUnstructuredVolumeObject( const local::ZhongVolumeObject* object );

private:
    size_t m_dim; ///< resolution
    size_t m_dim_r; ///< resolution in radius

public:
    ZhongVolumeObject();

    void setDim( const size_t dim ) { m_dim = dim; }
    void setDimR( const size_t dim_r ) { m_dim_r = dim_r; }

    size_t dim() const { return m_dim; }
    size_t dimR() const { return m_dim_r; }

    size_t numberOfNodes() const;
    size_t numberOfCells() const;
    void calculateCoords();
    bool readValues( const std::string& filename );
};

} // end of namespace local
