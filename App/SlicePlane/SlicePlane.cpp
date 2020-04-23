#include "SlicePlane.h"
#include <kvs/MarchingHexahedraTable>


namespace
{

bool HasIgnoreValue( const kvs::AnyValueArray values, const size_t* indices, const double ignore_value = 0.0 )
{
    for ( size_t i = 0; i < 8; ++i )
    {
        if ( kvs::Math::Equal( values.at<double>( indices[i] ), ignore_value ) ) { return true; }
    }
    return false;
}

} // end of namespace


namespace YYZVis
{

SlicePlane::SlicePlane():
    kvs::MapperBase(),
    kvs::PolygonObject()
{
}

SlicePlane::SlicePlane(
    const kvs::VolumeObjectBase* volume,
    const kvs::Vec4& coefficients,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->setPlane( coefficients );
    this->exec( volume );
}

SlicePlane::SlicePlane(
    const kvs::VolumeObjectBase* volume,
    const kvs::Vec3& point,
    const kvs::Vec3& normal,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject()
{
    this->setPlane( point, normal );
    this->exec( volume );
}

void SlicePlane::setPlane( const kvs::Vector4f& coefficients )
{
    m_coefficients = coefficients;
}

void SlicePlane::setPlane( const kvs::Vec3& point, const kvs::Vec3& normal )
{
    m_coefficients = kvs::Vec4( normal, -point.dot( normal ) );
}

SlicePlane::SuperClass* SlicePlane::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input object is NULL." << std::endl;
        return NULL;
    }

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input object is not volume dat." << std::endl;
        return NULL;
    }

    // Check whether the volume can be processed or not.
    if ( volume->veclen() != 1 )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "The input volume is not a sclar field data." << std::endl;
        return NULL;
    }

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );

    if ( YYZVis::ZhongVolumeObject::DownCast( volume ) )
    {
        const auto* zvolume = YYZVis::ZhongVolumeObject::DownCast( volume );
        this->mapping( zvolume );
    }
    else if ( YYZVis::YinYangVolumeObject::DownCast( volume ) )
    {
        const auto* yvolume = YYZVis::YinYangVolumeObject::DownCast( volume );
        this->mapping( yvolume );
    }

    return this;
}

void SlicePlane::mapping( const YYZVis::ZhongVolumeObject* zvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;
    // std::vector<kvs::UInt8> colors;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // Calculate coords, normals and colors.
     this->extract_plane( zvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void SlicePlane::mapping( const YYZVis::YinYangVolumeObject* yvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;
    // std::vector<kvs::UInt8> colors;
    // const kvs::ColorMap cmap = BaseClass::colorMap();

    // Calculate coords, normals and colors.
     this->extract_plane( yvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void SlicePlane::extract_plane( const YYZVis::YinYangVolumeObject* yvolume )
{
   // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;
    std::vector<kvs::UInt8> colors;

    const kvs::Real64 min_value( yvolume->minValue() );
    const kvs::Real64 max_value( yvolume->maxValue() );
    const kvs::Real64 normalize_factor( 255.0 / ( max_value - min_value ) );

    const size_t dim_r = yvolume->dimR(); // radius
    const size_t dim_theta = yvolume->dimTheta(); // latitude
    const size_t dim_phi = yvolume->dimPhi(); // longitude
    const size_t line_size = dim_r;
    const size_t slice_size = dim_r * dim_theta;

    const kvs::Real32* yvolume_coords = yvolume->coords().data();

    const kvs::ColorMap& color_map( BaseClass::transferFunction().colorMap() );

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( size_t k = 0; k < dim_phi - 1; k++, index += dim_r )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, index++ )
            {
                // Calculate the indices of the target cell.
                local_index[0] = index + line_size;
                local_index[1] = index + line_size + 1;
                local_index[2] = index + line_size + slice_size + 1;
                local_index[3] = index + line_size + slice_size;
                local_index[4] = index;
                local_index[5] = index + 1;
                local_index[6] = index + slice_size + 1;
                local_index[7] = index + slice_size;

                // Calculate the index of the reference table.
                const size_t table_index = this->calculate_hexahedra_table_index( local_index );
                if ( table_index == 0 ) continue;
                if ( table_index == 255 ) continue;

                // Calculate the triangle polygons.
                for ( size_t t = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][t] != -1; t += 3 )
                {
                    // Refer the edge IDs from the TriangleTable using the table_index.
                    const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][t];
                    const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+1];
                    const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+2];

                    // Refer indices of the coordinate array from the VertexTable using the edgeIDs.
                    const size_t c0 = local_index[ kvs::MarchingHexahedraTable::VertexID[e0][0] ];
                    const size_t c1 = local_index[ kvs::MarchingHexahedraTable::VertexID[e0][1] ];
                    const size_t c2 = local_index[ kvs::MarchingHexahedraTable::VertexID[e1][0] ];
                    const size_t c3 = local_index[ kvs::MarchingHexahedraTable::VertexID[e1][1] ];
                    const size_t c4 = local_index[ kvs::MarchingHexahedraTable::VertexID[e2][0] ];
                    const size_t c5 = local_index[ kvs::MarchingHexahedraTable::VertexID[e2][1] ];

                    // Determine vertices for each edge.
                    const kvs::Vec3 v0( yvolume_coords + 3 * c0 );
                    const kvs::Vec3 v1( yvolume_coords + 3 * c1 );

                    const kvs::Vec3 v2( yvolume_coords + 3 * c2 );
                    const kvs::Vec3 v3( yvolume_coords + 3 * c3 );

                    const kvs::Vec3 v4( yvolume_coords + 3 * c4 );
                    const kvs::Vec3 v5( yvolume_coords + 3 * c5 );

                    // Calculate coordinates of the vertices which are composed
                    // of the triangle polygon.
                    const kvs::Vec3 vertex0( this->interpolate_vertex( v0, v1 ) );
                    coords.push_back( vertex0.x() );
                    coords.push_back( vertex0.y() );
                    coords.push_back( vertex0.z() );

                    const kvs::Vec3 vertex1( this->interpolate_vertex( v2, v3 ) );
                    coords.push_back( vertex1.x() );
                    coords.push_back( vertex1.y() );
                    coords.push_back( vertex1.z() );

                    const kvs::Vec3 vertex2( this->interpolate_vertex( v4, v5 ) );
                    coords.push_back( vertex2.x() );
                    coords.push_back( vertex2.y() );
                    coords.push_back( vertex2.z() );

                    const double value0 = this->interpolate_value( yvolume, c0, c1 );
                    const double value1 = this->interpolate_value( yvolume, c2, c3 );
                    const double value2 = this->interpolate_value( yvolume, c4, c5 );

                    const kvs::UInt8 color0 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value0 - min_value ) );
                    colors.push_back( color_map[ color0 ].r() );
                    colors.push_back( color_map[ color0 ].g() );
                    colors.push_back( color_map[ color0 ].b() );

                    const kvs::UInt8 color1 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value1 - min_value ) );
                    colors.push_back( color_map[ color1 ].r() );
                    colors.push_back( color_map[ color1 ].g() );
                    colors.push_back( color_map[ color1 ].b() );

                    const kvs::UInt8 color2 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value2 - min_value ) );
                    colors.push_back( color_map[ color2 ].r() );
                    colors.push_back( color_map[ color2 ].g() );
                    colors.push_back( color_map[ color2 ].b() );

                    // Calculate a normal vector for the triangle polygon.
                    const kvs::Vec3 normal( -( vertex2 - vertex0 ).cross( vertex1 - vertex0 ) );
                    normals.push_back( normal.x() );
                    normals.push_back( normal.y() );
                    normals.push_back( normal.z() );
                } // end of loop-triangle
            }
        }
    } // end of loop-cell

    SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
    SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
    SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
    SuperClass::setOpacity( 255 );
    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );
}

void SlicePlane::extract_plane( const YYZVis::ZhongVolumeObject* zvolume )
{
    // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;
    std::vector<kvs::UInt8> colors;

    const kvs::Real64 min_value( zvolume->minValue() );
    const kvs::Real64 max_value( zvolume->maxValue() );
    const kvs::Real64 normalize_factor( 255.0 / ( max_value - min_value ) );

    const size_t dim = zvolume->dim();
    const size_t line_size = dim;
    const size_t slice_size = dim * dim;

    const kvs::Real32* zvolume_coords = zvolume->coords().data();

    const kvs::ColorMap& color_map( BaseClass::transferFunction().colorMap() );

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( size_t k = 0; k < dim - 1; k++, index += dim )
    {
        for ( size_t j = 0; j < dim - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim - 1; i++, index++ )
            {
                // Calculate the indices of the target cell.
                local_index[0] = index + line_size;
                local_index[1] = index + line_size + 1;
                local_index[2] = index + line_size + slice_size + 1;
                local_index[3] = index + line_size + slice_size;
                local_index[4] = index;
                local_index[5] = index + 1;
                local_index[6] = index + slice_size + 1;
                local_index[7] = index + slice_size;

                if ( ::HasIgnoreValue( zvolume->values(), local_index, 0.0 ) ) { continue; }

                // Calculate the index of the reference table.
                const size_t table_index = this->calculate_hexahedra_table_index( local_index );
                if ( table_index == 0 ) continue;
                if ( table_index == 255 ) continue;

                // Calculate the triangle polygons.
                for ( size_t t = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][t] != -1; t += 3 )
                {
                    // Refer the edge IDs from the TriangleTable using the table_index.
                    const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][t];
                    const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+1];
                    const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+2];

                    // Refer indices of the coordinate array from the VertexTable using the edgeIDs.
                    const size_t c0 = local_index[ kvs::MarchingHexahedraTable::VertexID[e0][0] ];
                    const size_t c1 = local_index[ kvs::MarchingHexahedraTable::VertexID[e0][1] ];
                    const size_t c2 = local_index[ kvs::MarchingHexahedraTable::VertexID[e1][0] ];
                    const size_t c3 = local_index[ kvs::MarchingHexahedraTable::VertexID[e1][1] ];
                    const size_t c4 = local_index[ kvs::MarchingHexahedraTable::VertexID[e2][0] ];
                    const size_t c5 = local_index[ kvs::MarchingHexahedraTable::VertexID[e2][1] ];

                    // Determine vertices for each edge.
                    const kvs::Vec3 v0( zvolume_coords + 3 * c0 );
                    const kvs::Vec3 v1( zvolume_coords + 3 * c1 );

                    const kvs::Vec3 v2( zvolume_coords + 3 * c2 );
                    const kvs::Vec3 v3( zvolume_coords + 3 * c3 );

                    const kvs::Vec3 v4( zvolume_coords + 3 * c4 );
                    const kvs::Vec3 v5( zvolume_coords + 3 * c5 );

                    // Calculate coordinates of the vertices which are composed
                    // of the triangle polygon.
                    const kvs::Vec3 vertex0( this->interpolate_vertex( v0, v1 ) );
                    coords.push_back( vertex0.x() );
                    coords.push_back( vertex0.y() );
                    coords.push_back( vertex0.z() );

                    const kvs::Vec3 vertex1( this->interpolate_vertex( v2, v3 ) );
                    coords.push_back( vertex1.x() );
                    coords.push_back( vertex1.y() );
                    coords.push_back( vertex1.z() );

                    const kvs::Vec3 vertex2( this->interpolate_vertex( v4, v5 ) );
                    coords.push_back( vertex2.x() );
                    coords.push_back( vertex2.y() );
                    coords.push_back( vertex2.z() );

                    const double value0 = this->interpolate_value( zvolume, c0, c1 );
                    const double value1 = this->interpolate_value( zvolume, c2, c3 );
                    const double value2 = this->interpolate_value( zvolume, c4, c5 );

                    const kvs::UInt8 color0 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value0 - min_value ) );
                    colors.push_back( color_map[ color0 ].r() );
                    colors.push_back( color_map[ color0 ].g() );
                    colors.push_back( color_map[ color0 ].b() );

                    const kvs::UInt8 color1 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value1 - min_value ) );
                    colors.push_back( color_map[ color1 ].r() );
                    colors.push_back( color_map[ color1 ].g() );
                    colors.push_back( color_map[ color1 ].b() );

                    const kvs::UInt8 color2 =
                        static_cast<kvs::UInt8>( normalize_factor * ( value2 - min_value ) );
                    colors.push_back( color_map[ color2 ].r() );
                    colors.push_back( color_map[ color2 ].g() );
                    colors.push_back( color_map[ color2 ].b() );

                    // Calculate a normal vector for the triangle polygon.
                    const kvs::Vec3 normal( -( vertex2 - vertex0 ).cross( vertex1 - vertex0 ) );
                    normals.push_back( normal.x() );
                    normals.push_back( normal.y() );
                    normals.push_back( normal.z() );
                } // end of loop-triangle
            }
        }
    } // end of loop-cell

    SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
    SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
    SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
    SuperClass::setOpacity( 255 );
    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );
}

size_t SlicePlane::calculate_hexahedra_table_index( const size_t* local_index ) const
{
    const kvs::Real32* const coords = BaseClass::volume()->coords().data();

    const kvs::Vec3 vertex0( coords + 3 * local_index[0] );
    const kvs::Vec3 vertex1( coords + 3 * local_index[1] );
    const kvs::Vec3 vertex2( coords + 3 * local_index[2] );
    const kvs::Vec3 vertex3( coords + 3 * local_index[3] );
    const kvs::Vec3 vertex4( coords + 3 * local_index[4] );
    const kvs::Vec3 vertex5( coords + 3 * local_index[5] );
    const kvs::Vec3 vertex6( coords + 3 * local_index[6] );
    const kvs::Vec3 vertex7( coords + 3 * local_index[7] );

    size_t table_index = 0;
    if ( this->substitute_plane_equation( vertex0 ) > 0.0 ) { table_index |=   1; }
    if ( this->substitute_plane_equation( vertex1 ) > 0.0 ) { table_index |=   2; }
    if ( this->substitute_plane_equation( vertex2 ) > 0.0 ) { table_index |=   4; }
    if ( this->substitute_plane_equation( vertex3 ) > 0.0 ) { table_index |=   8; }
    if ( this->substitute_plane_equation( vertex4 ) > 0.0 ) { table_index |=  16; }
    if ( this->substitute_plane_equation( vertex5 ) > 0.0 ) { table_index |=  32; }
    if ( this->substitute_plane_equation( vertex6 ) > 0.0 ) { table_index |=  64; }
    if ( this->substitute_plane_equation( vertex7 ) > 0.0 ) { table_index |= 128; }

    return table_index;
}

float SlicePlane::substitute_plane_equation( const kvs::Vec3& vertex ) const
{
    return
        m_coefficients.x() * vertex.x() +
        m_coefficients.y() * vertex.y() +
        m_coefficients.z() * vertex.z() +
        m_coefficients.w();
}

const kvs::Vec3 SlicePlane::interpolate_vertex(
    const kvs::Vec3& vertex0,
    const kvs::Vec3& vertex1 ) const
{
    const float value0 = this->substitute_plane_equation( vertex0 );
    const float value1 = this->substitute_plane_equation( vertex1 );
    const float ratio = kvs::Math::Abs( value0 / ( value1 - value0 ) );

    return ( 1.0f - ratio ) * vertex0 + ratio * vertex1;
}

double SlicePlane::interpolate_value(
    const YYZVis::YinYangVolumeObject* yvolume,
    const size_t index0,
    const size_t index1 ) const
{
    const kvs::AnyValueArray values = yvolume->values();
    const kvs::Real32* const coords = yvolume->coords().data();

    const float value0 = this->substitute_plane_equation( kvs::Vec3( coords + 3 * index0 ) );
    const float value1 = this->substitute_plane_equation( kvs::Vec3( coords + 3 * index1 ) );
    const float ratio = kvs::Math::Abs( value0 / ( value1 - value0 ) );

    return values.at<float>( index0 ) + ratio * ( values.at<float>( index1 ) - values.at<float>( index0 ) );
}

double SlicePlane::interpolate_value(
    const YYZVis::ZhongVolumeObject* zvolume,
    const size_t index0,
    const size_t index1 ) const
{
    const kvs::AnyValueArray values = zvolume->values();
    const kvs::Real32* const coords = zvolume->coords().data();

    const float value0 = this->substitute_plane_equation( kvs::Vec3( coords + 3 * index0 ) );
    const float value1 = this->substitute_plane_equation( kvs::Vec3( coords + 3 * index1 ) );
    const float ratio = kvs::Math::Abs( value0 / ( value1 - value0 ) );

    return values.at<float>( index0 ) + ratio * ( values.at<float>( index1 ) - values.at<float>( index0 ) );
}

} // end of namespace YYZVis
