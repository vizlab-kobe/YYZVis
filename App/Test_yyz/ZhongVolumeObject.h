#pragma once
#include <kvs/Module>
#include <kvs/VolumeObjectBase>
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

public:
    struct Range
    {
        float min;
        float max;
        float d;
    };

private:
    size_t m_dim; ///< resolution
    size_t m_dim_r; ///< resolution in radius
    Range m_range_r; ///< range of radius

public:
    ZhongVolumeObject();

    void setDim( const local::ZhongVolumeObject* object, const float boundary_r = 0.35f );
    void setDimR( const size_t dim_r, const float range_min = 0.35f, const float range_max = 1.0f );

    size_t dim() const { return m_dim; }
    size_t dimR() const { return m_dim_r; }
    Range rangeR() const { return m_range_r; }

    size_t numberOfNodes() const;
    size_t numberOfCells() const;
    bool readValues( const std::string& filename );
};

} // end of namespace local
