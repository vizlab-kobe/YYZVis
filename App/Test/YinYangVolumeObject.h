#pragma once
#include <kvs/Module>
#include <kvs/VolumeObjectBase>
//#include <kvs/Vector3>
#include <kvs/StructuredVolumeObject>
#include <kvs/UnstructuredVolumeObject>


namespace local
{

class YinYangVolumeObject : public kvs::VolumeObjectBase
{
    kvsModule( local::YinYangVolumeObject, Object );
    kvsModuleBaseClass( kvs::VolumeObjectBase );

public:
    static kvs::StructuredVolumeObject* ToStructuredVolumeObject( const local::YinYangVolumeObject* object );
    static kvs::UnstructuredVolumeObject* ToUnstructuredVolumeObject( const local::YinYangVolumeObject* object );

public:
    enum GridType
    {
        Yin,
        Yang
    };

private:
    GridType m_grid_type; ///< grid type
//    kvs::Vec3ui m_resolution; ///< resolution in (r, theta, phi).
    size_t m_dim_r; ///< resolution in radius
    size_t m_dim_theta; ///< resolution in latitude
    size_t m_dim_phi; ///< resolution in longitude

public:
    YinYangVolumeObject();

    void setGridType( const GridType grid_type ) { m_grid_type = grid_type; }
    void setGridTypeToYin() { this->setGridType( Yin ); }
    void setGridTypeToYang() { this->setGridType( Yang ); }
//    void setResolution( const kvs::Vec3ui& resolution ) { m_resolution = resolution; }
    void setDimR( const size_t dim_r ) { m_dim_r = dim_r; }
    void setDimTheta( const size_t dim_theta ) { m_dim_theta = dim_theta; }
    void setDimPhi( const size_t dim_phi ) { m_dim_phi = dim_phi; }

    GridType gridType() const { return m_grid_type; }
//    const kvs::Vec3ui& resolution() const { return m_resolution; }
    size_t dimR() const { return m_dim_r; }
    size_t dimTheta() const { return m_dim_theta; }
    size_t dimPhi() const { return m_dim_phi; }

    size_t numberOfNodes() const;
    size_t numberOfCells() const;
    void calculateCoords();
    bool readValues( const std::string& filename );
};

} // end of namespace local
