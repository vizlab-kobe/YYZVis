#pragma once
#include <kvs/ValueArray>
#include <kvs/ObjectBase>
#include <kvs/Camera>
#include <kvs/OpacityMap>


namespace local
{

class DensityMap
{
public:
    typedef kvs::ValueArray<kvs::Real32> Table;

private:
    size_t m_resolution; ///< table resolution
    kvs::Real32 m_min_value; ///< min. value
    kvs::Real32 m_max_value; ///< max. value
    Table m_table; ///< value table
    kvs::Real32 m_sampling_step; ///< length of the ray segment (dt)
    const kvs::Camera* m_camera; ///< pointer to the referenced camera
    const kvs::ObjectBase* m_object; ///< pointer to the referenced object

public:
    void setSamplingStep( const kvs::Real32 step ) { m_sampling_step = step; }
    void attachCamera( const kvs::Camera* camera ) { m_camera = camera; }
    void attachObject( const kvs::ObjectBase* object ) { m_object = object; }
    size_t resolution() const { return m_resolution; }
    kvs::Real32 minValue() const { return m_min_value; }
    kvs::Real32 maxValue() const { return m_max_value; }
    const Table& table() const { return m_table; }

    kvs::Real32 at( const kvs::Real32 value ) const;
    void create( const kvs::OpacityMap& omap );

private:
    kvs::Real32 max_density( const kvs::Real32 s0, const kvs::Real32 s1 ) const;
};

} // end of namespace local
