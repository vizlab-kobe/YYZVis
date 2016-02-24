#pragma once
#include "ZhongVolumeObject.h"


namespace local
{

class ZhongGrid
{
private:
    kvs::Vec3ui m_base_index; ///< base index of the bound grid
    kvs::Vec3 m_coords[8];
    kvs::Real32 m_values[8];
    mutable kvs::Real32 m_interpolation_functions[8];
    mutable kvs::Real32 m_differential_functions[24];
    mutable kvs::Vec3 m_local_point; ///< local point
    const local::ZhongVolumeObject* m_reference_volume;

public:
    ZhongGrid( const local::ZhongVolumeObject* volume );

    const kvs::Vec3ui& baseIndex() const { return m_base_index; }
    const kvs::Vec3& coord( const size_t index ) const { return m_coords[index]; }
    const kvs::Real32 value( const size_t index ) const { return m_values[index]; }
    const kvs::Vec3& localPoint() const { return m_local_point; }
    const kvs::Real32* interpolationFunctions() const { return m_interpolation_functions; }
    const kvs::Real32* differentialFunctions() const { return m_differential_functions; }

    void bind( const kvs::Vec3ui& base_index );
    void setLocalPoint( const kvs::Vec3& local ) const;
    void updateInterpolationFunctions( const kvs::Vec3& local ) const;
    void updateDifferentialFunctions( const kvs::Vec3& local ) const;
    const kvs::Vec3 globalPoint() const;
    const kvs::Mat3 JacobiMatrix() const;
    const kvs::Vec3 center() const;
    const kvs::Real32 volume() const;
    const kvs::Real32 scalar() const;
    const kvs::Vec3 gradientVector() const;
    const kvs::Vec3 randomSampling() const;
};

} // end of namespace local
