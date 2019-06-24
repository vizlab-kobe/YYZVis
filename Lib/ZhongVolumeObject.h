#pragma once
#include <ostream>
#include <kvs/Module>
#include <kvs/VolumeObjectBase>
#include <kvs/StructuredVolumeObject>
#include <kvs/UnstructuredVolumeObject>


namespace YinYangVis
{

class ZhongVolumeObject : public kvs::VolumeObjectBase
{
    kvsModule( YinYangVis::ZhongVolumeObject, Object );
    kvsModuleBaseClass( kvs::VolumeObjectBase );

public:
    static kvs::StructuredVolumeObject* ToStructuredVolumeObject( const YinYangVis::ZhongVolumeObject* object );
    static kvs::UnstructuredVolumeObject* ToUnstructuredVolumeObject( const YinYangVis::ZhongVolumeObject* object );

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
    ZhongVolumeObject( const ZhongVolumeObject& object ) { this->shallowCopy( object ); }

    void shallowCopy( const ZhongVolumeObject& object );
    void deepCopy( const ZhongVolumeObject& object );
    void print( std::ostream& os, const kvs::Indent& indent = kvs::Indent(0) ) const;

    void setDim( const size_t dim ) { m_dim = dim; }
    void setDimR( const size_t dim_r, const float range_min = 0.35f, const float range_max = 1.0f );

    size_t dim() const { return m_dim; }
    size_t dimR() const { return m_dim_r; }
    Range rangeR() const { return m_range_r; }

    size_t numberOfNodes() const;
    size_t numberOfCells() const;
    void calculateCoords();
    bool readValues( const std::string& filename );
    void updateMinMaxCoords();
};

} // end of namespace YinYangVis
