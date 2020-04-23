#include "Isosurface.h"
#include <kvs/MarchingHexahedraTable>

namespace YYZVis
{

Isosurface::Isosurface():
    kvs::MapperBase(),
    kvs::PolygonObject(),
    m_isolevel( 0 ),
    m_duplication( true )
{
}

Isosurface::Isosurface(
    const kvs::VolumeObjectBase* volume,
    const double isolevel,
    const SuperClass::NormalType normal_type,
    const bool duplication,
    const kvs::TransferFunction& transfer_function ):
    kvs::MapperBase( transfer_function ),
    kvs::PolygonObject(),
    m_duplication( duplication )
{
    SuperClass::setNormalType( normal_type );
    this->setIsolevel( isolevel );
    this->exec( volume );
}

Isosurface::SuperClass* Isosurface::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is NULL.");
        return NULL;
    }

    const kvs::VolumeObjectBase* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("Input object is not volume dat.");
        return NULL;
    }

    // Check whether the volume can be processed or not.
    if ( volume->veclen() != 1 )
    {
        BaseClass::setSuccess( false );
        kvsMessageError("The input volume is not a sclar field data.");
        return NULL;
    }

    // In the case of VertexNormal-type, the duplicated vertices are forcibly deleted.
    if ( SuperClass::normalType() == kvs::PolygonObject::VertexNormal )
    {
        m_duplication = false;
    }

    BaseClass::attachVolume( volume );
    BaseClass::setRange( volume );
    BaseClass::setMinMaxCoords( volume, this );

    const kvs::Real64 min_value = BaseClass::volume()->minValue();
    const kvs::Real64 max_value = BaseClass::volume()->maxValue();
    if ( kvs::Math::Equal( min_value, max_value ) ) { return this; }

    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );

    if ( YYZVis::ZhongVolumeObject::DownCast( volume ) )
    {
        const YYZVis::ZhongVolumeObject* zvolume = YYZVis::ZhongVolumeObject::DownCast( volume );
        this->mapping( zvolume );
    }
    else if ( YYZVis::YinYangVolumeObject::DownCast( volume ) )
    {
        const YYZVis::YinYangVolumeObject* yvolume = YYZVis::YinYangVolumeObject::DownCast( volume );
        this->mapping( yvolume );
    }

    return this;
}

void Isosurface::mapping( const YYZVis::ZhongVolumeObject* zvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
     if ( m_duplication )
     {
        this->extract_zhong_surfaces_with_duplication( zvolume );
     }
    // else
    // {
    //     ...
    // }

    // const kvs::ColorMap cmap = BaseClass::colorMap();
    // const kvs::RGBColor color = COLOR_OF_ISOSURFACES;

    // SuperClass::setCoords( coords );
    // SuperClass::setColor( color );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void Isosurface::mapping( const YYZVis::YinYangVolumeObject* yvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
     if ( m_duplication )
     {
       this->extract_yinyang_surfaces_with_duplication( yvolume );
     }
    // else
    // {
    //     ...
    // }

    // const kvs::ColorMap cmap = BaseClass::colorMap();
    // const kvs::RGBColor color = COLOR_OF_ISOSURFACES;

    // SuperClass::setCoords( coords );
    // SuperClass::setColor( color );
    // SuperClass::setNormals( normals );
    // SuperClass::setOpacity( 255 );
}

void Isosurface::extract_yinyang_surfaces_with_duplication(const YYZVis::YinYangVolumeObject* yvolume )
{
    // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;

    const size_t dim_r = yvolume->dimR(); // radius
    const size_t dim_theta = yvolume->dimTheta(); // latitude
    const size_t dim_phi = yvolume->dimPhi(); // longitude

    const kvs::AnyValueArray values = yvolume->values();
    const kvs::Real32* const yvolume_coords = yvolume->coords().data();
    
    // Calculate connections.
    const size_t nnodes = yvolume->numberOfNodes();
    kvs::ValueArray<kvs::UInt32> connections( nnodes * 8 );
    kvs::UInt32* pconnections = connections.data();
    for ( size_t k = 0, index = 0; k < dim_phi - 1; k++, index += dim_r )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, index++ )
            {
                *(pconnections++) = index;
                *(pconnections++) = index + 1;
                *(pconnections++) = index + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + ( dim_r * dim_theta );
                *(pconnections++) = index + dim_r;
                *(pconnections++) = index + dim_r + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta ) + 1;
                *(pconnections++) = index + dim_r + ( dim_r * dim_theta );
            }
        }
    }
    
    const kvs::UInt32 ncells( yvolume->numberOfCells() );

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( kvs::UInt32 cell = 0; cell < ncells; ++cell, index += 8 )
    {
        // Calculate the indices of the target cell.
        local_index[0] = connections[ index + 4 ];
        local_index[1] = connections[ index + 5 ];
        local_index[2] = connections[ index + 6 ];
        local_index[3] = connections[ index + 7 ];
        local_index[4] = connections[ index + 0 ];
        local_index[5] = connections[ index + 1 ];
        local_index[6] = connections[ index + 2 ];
        local_index[7] = connections[ index + 3 ];

        // Calculate the index of the reference table.
        const size_t table_index = this->calculate_table_index( values, local_index );
        if ( table_index == 0 ) continue;
        if ( table_index == 255 ) continue;

        // Calculate the triangle polygons.
        for ( size_t i = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][i] != -1; i += 3 )
        {
            // Refer the edge IDs from the TriangleTable by using the table_index.
	  const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][i];
	  const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][i+2];
	  const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][i+1];

            // Determine vertices for each edge.
	  const int v0 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][0]];
	  const int v1 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][1]];

	  const int v2 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][0]];
	  const int v3 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][1]];

	  const int v4 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][0]];
	  const int v5 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][1]];

            // Calculate coordinates of the vertices which are composed
            // of the triangle polygon.
       	    const kvs::Vec3 vertex0( this->interpolate_vertex( values, yvolume_coords, v0, v1 ) );
            coords.push_back( vertex0.x() );
            coords.push_back( vertex0.y() );
            coords.push_back( vertex0.z() );

            const kvs::Vec3 vertex1( this->interpolate_vertex( values, yvolume_coords, v2, v3 ) );
            coords.push_back( vertex1.x() );
            coords.push_back( vertex1.y() );
            coords.push_back( vertex1.z() );

            const kvs::Vec3 vertex2( this->interpolate_vertex( values, yvolume_coords, v4, v5 ) );
            coords.push_back( vertex2.x() );
            coords.push_back( vertex2.y() );
            coords.push_back( vertex2.z() );

            // Calculate a normal vector for the triangle polygon.
            const kvs::Vec3 normal( ( vertex1 - vertex0 ).cross( vertex2 - vertex0 ) );
            normals.push_back( normal.x() );
            normals.push_back( normal.y() );
            normals.push_back( normal.z() );
        } // end of loop-triangle
    } // end of loop-cell

    // Calculate the polygon color for the isolevel.
    const kvs::RGBColor color = this->calculate_color();

    if( coords.size() > 0 ){
        SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
        SuperClass::setColor( color );
        SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
        SuperClass::setOpacity( 255 );
        SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
        SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
        SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );
    }
}

void Isosurface::extract_zhong_surfaces_with_duplication(const YYZVis::ZhongVolumeObject* zvolume )
{
    // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;

    const size_t dim_x = zvolume->dim(); // X
    const size_t dim_y = zvolume->dim(); // Y
    const size_t dim_z = zvolume->dim(); // Z

    const kvs::AnyValueArray values = zvolume->values();    
    const kvs::Real32* const zvolume_coords = zvolume->coords().data();
    
    // Calculate connections.
    const size_t nnodes = zvolume->numberOfNodes();
    kvs::ValueArray<kvs::UInt32> connections( nnodes * 8 );
    kvs::UInt32* pconnections = connections.data();
    for ( size_t k = 0, index = 0; k < dim_z - 1; k++, index += dim_x )
    {
        for ( size_t j = 0; j < dim_y - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim_x - 1; i++, index++ )
            {
                *(pconnections++) = index;
                *(pconnections++) = index + 1;
                *(pconnections++) = index + ( dim_x * dim_y ) + 1;
                *(pconnections++) = index + ( dim_x * dim_y );
                *(pconnections++) = index + dim_x;
                *(pconnections++) = index + dim_x + 1;
                *(pconnections++) = index + dim_x + ( dim_x * dim_y ) + 1;
                *(pconnections++) = index + dim_x + ( dim_x * dim_y );
            }
        }
    }
    
    const kvs::UInt32 ncells( zvolume->numberOfCells() );

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( kvs::UInt32 cell = 0; cell < ncells; ++cell, index += 8 )
    {
        // Calculate the indices of the target cell.
        local_index[0] = connections[ index + 4 ];
        local_index[1] = connections[ index + 5 ];
        local_index[2] = connections[ index + 6 ];
        local_index[3] = connections[ index + 7 ];
        local_index[4] = connections[ index + 0 ];
        local_index[5] = connections[ index + 1 ];
        local_index[6] = connections[ index + 2 ];
        local_index[7] = connections[ index + 3 ];

        // Calculate the index of the reference table.
        const size_t table_index = this->calculate_table_index( values, local_index );
        if ( table_index == 0 ) continue;
        if ( table_index == 255 ) continue;

        // Calculate the triangle polygons.
        for ( size_t i = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][i] != -1; i += 3 )
        {
            // Refer the edge IDs from the TriangleTable by using the table_index.
	  const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][i];
	  const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][i+2];
	  const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][i+1];

            // Determine vertices for each edge.
	  const int v0 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][0]];
	  const int v1 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][1]];

	  const int v2 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][0]];
	  const int v3 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][1]];

	  const int v4 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][0]];
	  const int v5 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][1]];

            // Calculate coordinates of the vertices which are composed
            // of the triangle polygon.
	  const kvs::Vec3 vertex0( this->interpolate_vertex( values, zvolume_coords, v0, v1 ) );
            coords.push_back( vertex0.x() );
            coords.push_back( vertex0.y() );
            coords.push_back( vertex0.z() );

            const kvs::Vec3 vertex1( this->interpolate_vertex( values, zvolume_coords, v2, v3 ) );
            coords.push_back( vertex1.x() );
            coords.push_back( vertex1.y() );
            coords.push_back( vertex1.z() );

            const kvs::Vec3 vertex2( this->interpolate_vertex( values, zvolume_coords, v4, v5 ) );
            coords.push_back( vertex2.x() );
            coords.push_back( vertex2.y() );
            coords.push_back( vertex2.z() );

            // Calculate a normal vector for the triangle polygon.
            const kvs::Vec3 normal( ( vertex1 - vertex0 ).cross( vertex2 - vertex0 ) );
            normals.push_back( normal.x() );
            normals.push_back( normal.y() );
            normals.push_back( normal.z() );
        } // end of loop-triangle
    } // end of loop-cell

    // Calculate the polygon color for the isolevel.
    const kvs::RGBColor color = this->calculate_color();

    if( coords.size() > 0 ){
        SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
        SuperClass::setColor( color );
        SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
        SuperClass::setOpacity( 255 );
        SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
        SuperClass::setColorType( kvs::PolygonObject::PolygonColor );
        SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );
    }
}

size_t Isosurface::calculate_table_index( const kvs::AnyValueArray values, const size_t* local_index ) const
{
    const double isolevel = m_isolevel;

    size_t table_index = 0;
    if ( values.at<double>( local_index[0] ) > isolevel ) { table_index |=   1; }
    if ( values.at<double>( local_index[1] ) > isolevel ) { table_index |=   2; }
    if ( values.at<double>( local_index[2] ) > isolevel ) { table_index |=   4; }
    if ( values.at<double>( local_index[3] ) > isolevel ) { table_index |=   8; }
    if ( values.at<double>( local_index[4] ) > isolevel ) { table_index |=  16; }
    if ( values.at<double>( local_index[5] ) > isolevel ) { table_index |=  32; }
    if ( values.at<double>( local_index[6] ) > isolevel ) { table_index |=  64; }
    if ( values.at<double>( local_index[7] ) > isolevel ) { table_index |= 128; }

    return table_index;
}


const kvs::Vec3 Isosurface::interpolate_vertex(
				                const kvs::AnyValueArray values,
						const kvs::Real32* coords,
						const int vertex0,
                                                const int vertex1 ) const
{
    const size_t coord0_index = 3 * vertex0;
    const size_t coord1_index = 3 * vertex1;

    const double v0 = values.at<double>( vertex0 );
    const double v1 = values.at<double>( vertex1 );
    const float ratio = static_cast<float>( kvs::Math::Abs( ( m_isolevel - v0 ) / ( v1 - v0 ) ) );

    const float x = coords[coord0_index]   + ratio * ( coords[coord1_index]   - coords[coord0_index] );
    const float y = coords[coord0_index+1] + ratio * ( coords[coord1_index+1] - coords[coord0_index+1] );
    const float z = coords[coord0_index+2] + ratio * ( coords[coord1_index+2] - coords[coord0_index+2] );

    return kvs::Vec3( x, y, z );
}

const kvs::RGBColor Isosurface::calculate_color()
{
    const kvs::Real64 min_value = BaseClass::volume()->minValue();
    const kvs::Real64 max_value = BaseClass::volume()->maxValue();
    const kvs::Real64 normalize_factor = 255.0 / ( max_value - min_value );
    const kvs::UInt8  index = static_cast<kvs::UInt8>( normalize_factor * ( m_isolevel - min_value ) );

    return BaseClass::transferFunction().colorMap()[ index ];
}
} // end of namespace YYZVis
