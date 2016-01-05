#pragma once
#include <kvs/Module>
#include <kvs/VolumeObjectBase>
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

    struct Range
    {
        float min;
        float max;
        float d;
    };

private:
    GridType m_grid_type; ///< grid type
    size_t m_dim_r; ///< resolution in radius
    size_t m_dim_theta; ///< resolution in latitude
    size_t m_dim_phi; ///< resolution in longitude
    Range m_range_r; ///< range of radius
    Range m_range_theta; ///< range of latitude
    Range m_range_phi; ///< range of longitude

public:
    YinYangVolumeObject();

    void setGridType( const GridType grid_type ) { m_grid_type = grid_type; }
    void setGridTypeToYin() { this->setGridType( Yin ); }
    void setGridTypeToYang() { this->setGridType( Yang ); }

    void setDimR( const size_t dim_r, const float range_min = 0.35f, const float range_max = 1.0f );
    void setDimTheta( const size_t dim_theta, const size_t overwrap = 2 );
    void setDimPhi( const size_t dim_phi, const size_t overwrap = 4 );

    GridType gridType() const { return m_grid_type; }
    size_t dimR() const { return m_dim_r; }
    size_t dimTheta() const { return m_dim_theta; }
    size_t dimPhi() const { return m_dim_phi; }
    Range rangeR() const { return m_range_r; }
    Range rangeTheta() const { return m_range_theta; }
    Range rangePhi() const { return m_range_phi; }

    size_t numberOfNodes() const;
    size_t numberOfCells() const;
    bool readValues( const std::string& filename );
};

} // end of namespace local
