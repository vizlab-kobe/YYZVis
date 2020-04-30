#include "Isosurface.h"
#include <kvs/MarchingHexahedraTable>
#include <kvs/MarchingCubesTable>
#include <kvs/Math>


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

/*===========================================================================*/
/**
 *  @brief  Constructs a new Isosurface class.
 */
/*===========================================================================*/
Isosurface::Isosurface():
    kvs::MapperBase(),
    kvs::PolygonObject(),
    m_isolevel( 0 ),
    m_duplication( true )
{
}

/*===========================================================================*/
/**
 *  @brief  Constructs a new Isosurface class.
 *  @param  volume [in] pointer to the volume object
 *  @param  isolevel [in] level of the isosurfaces
 *  @param  normal_type [in] type of the normal vector
 *  @param  duplication [in] duplication flag
 *  @param  transfer_function [in] transfer function
 */
/*===========================================================================*/
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

/*===========================================================================*/
/**
 *  @brief  Executes the isosurface extraction process.
 *  @param  object [in] pointer to the input volume object
 *  @return pointer to the polygon object
 */
/*===========================================================================*/
Isosurface::SuperClass* Isosurface::exec( const kvs::ObjectBase* object )
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
        const auto* zvolume = YYZVis::ZhongVolumeObject::DownCast( volume );
        this->mapping( zvolume );
    }
    else if ( YYZVis::YinYangVolumeObjectBase::DownCast( volume ) )
    {
        const auto* yvolume = YYZVis::YinYangVolumeObjectBase::DownCast( volume );
        this->mapping( yvolume );
    }

    return this;
}

/*==========================================================================*/
/**
 *  @brief  Extracts isosurfaces from the specified zhong volume object.
 *  @param  volume [in] pointer to zhong volume object
 */
/*==========================================================================*/
void Isosurface::mapping( const YYZVis::ZhongVolumeObject* zvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
//     if ( m_duplication )
     {
        this->extract_surfaces_with_duplication( zvolume );
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

/*==========================================================================*/
/**
 *  @brief  Extracts isosurfaces from the specified yin/yang volume object.
 *  @param  volume [in] pointer to yin/yang volume object
 */
/*==========================================================================*/
void Isosurface::mapping( const YYZVis::YinYangVolumeObjectBase* yvolume )
{
    // std::vector<kvs::Real32> coords;
    // std::vector<kvs::Real32> normals;

    // Calculate coords and normals.
//     if ( m_duplication )
     {
       this->extract_surfaces_with_duplication( yvolume );
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

void Isosurface::extract_surfaces_with_duplication( const YYZVis::YinYangVolumeObjectBase* yvolume )
{
    // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;

    const size_t dim_r = yvolume->dimR(); // radius
    const size_t dim_theta = yvolume->dimTheta(); // latitude
    const size_t dim_phi = yvolume->dimPhi(); // longitude
    const size_t line_size = dim_r;
    const size_t slice_size = dim_r * dim_theta;

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    for ( size_t k = 0; k < dim_phi - 1; k++, index += dim_r )
    {
        for ( size_t j = 0; j < dim_theta - 1; j++, index++ )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, index++ )
            {
                // Calculate the indices of the target cell.
                local_index[0] = index;
                local_index[1] = local_index[0] + 1;
                local_index[2] = local_index[1] + line_size;
                local_index[3] = local_index[0] + line_size;
                local_index[4] = local_index[0] + slice_size;
                local_index[5] = local_index[1] + slice_size;
                local_index[6] = local_index[2] + slice_size;
                local_index[7] = local_index[3] + slice_size;

                // Calculate the index of the reference table.
                const size_t table_index = this->calculate_table_index( yvolume->values(), local_index );
                if ( table_index == 0 ) continue;
                if ( table_index == 255 ) continue;

                // Calculate the triangle polygons.
                for ( size_t t = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][t] != -1; t += 3 )
                {
                    // Refer the edge IDs from the TriangleTable by using the table_index.
                    const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][t];
                    const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+2];
                    const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+1];

                    // Determine vertices for each edge.
                    const int p0 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][0]];
                    const int p1 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][1]];

                    const int p2 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][0]];
                    const int p3 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][1]];

                    const int p4 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][0]];
                    const int p5 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][1]];

                    // Calculate coordinates of the vertices which are composed
                    // of the triangle polygon.
                    const kvs::Vec3 v0( yvolume->coords().data() + p0 * 3 );
                    const kvs::Vec3 v1( yvolume->coords().data() + p1 * 3 );
                    const kvs::Real64 s0 = yvolume->values().at<kvs::Real64>( p0 );
                    const kvs::Real64 s1 = yvolume->values().at<kvs::Real64>( p1 );
                    const kvs::Vec3 vertex0( this->interpolate_vertex( v0, v1, s0, s1 ) );
                    coords.push_back( vertex0.x() );
                    coords.push_back( vertex0.y() );
                    coords.push_back( vertex0.z() );

                    const kvs::Vec3 v2( yvolume->coords().data() + p2 * 3 );
                    const kvs::Vec3 v3( yvolume->coords().data() + p3 * 3 );
                    const kvs::Real64 s2 = yvolume->values().at<kvs::Real64>( p2 );
                    const kvs::Real64 s3 = yvolume->values().at<kvs::Real64>( p3 );
                    const kvs::Vec3 vertex1( this->interpolate_vertex( v2, v3, s2, s3 ) );
                    coords.push_back( vertex1.x() );
                    coords.push_back( vertex1.y() );
                    coords.push_back( vertex1.z() );

                    const kvs::Vec3 v4( yvolume->coords().data() + p4 * 3 );
                    const kvs::Vec3 v5( yvolume->coords().data() + p5 * 3 );
                    const kvs::Real64 s4 = yvolume->values().at<kvs::Real64>( p4 );
                    const kvs::Real64 s5 = yvolume->values().at<kvs::Real64>( p5 );
                    const kvs::Vec3 vertex2( this->interpolate_vertex( v4, v5, s4, s5 ) );
                    coords.push_back( vertex2.x() );
                    coords.push_back( vertex2.y() );
                    coords.push_back( vertex2.z() );

                    // Calculate a normal vector for the triangle polygon.
                    const kvs::Vec3 normal( ( vertex1 - vertex0 ).cross( vertex2 - vertex0 ) );
                    normals.push_back( normal.x() );
                    normals.push_back( normal.y() );
                    normals.push_back( normal.z() );
                } // end of loop-triangle
            }
        }
    }

    // Calculate the polygon color for the isolevel.
    const kvs::RGBColor color = this->calculate_color();

    if ( coords.size() > 0 )
    {
        SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
        SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
        SuperClass::setColor( color );
        SuperClass::setOpacity( 255 );
        SuperClass::setPolygonTypeToTriangle();
        SuperClass::setColorTypeToPolygon();
        SuperClass::setNormalTypeToPolygon();
    }
}

void Isosurface::extract_surfaces_with_duplication( const YYZVis::ZhongVolumeObject* zvolume )
{
    // Calculated the coordinate data array and the normal vector array.
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;

    const size_t dim = zvolume->dim();
    const size_t line_size = dim;
    const size_t slice_size = dim * dim;

    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( size_t k = 0; k < dim - 1; k++ )
    {
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                // Calculate the indices of the target cell.
                local_index[0] = index;
                local_index[1] = local_index[0] + 1;
                local_index[2] = local_index[1] + line_size;
                local_index[3] = local_index[0] + line_size;
                local_index[4] = local_index[0] + slice_size;
                local_index[5] = local_index[1] + slice_size;
                local_index[6] = local_index[2] + slice_size;
                local_index[7] = local_index[3] + slice_size;
                index++;

                if ( ::HasIgnoreValue( zvolume->values(), local_index, 0.0 ) ) { continue; }

                // Calculate the index of the reference table.
                const size_t table_index = this->calculate_table_index( zvolume->values(), local_index );
                if ( table_index == 0 ) continue;
                if ( table_index == 255 ) continue;

                // Calculate the triangle polygons.
                for ( size_t t = 0; kvs::MarchingHexahedraTable::TriangleID[ table_index ][t] != -1; t += 3 )
                {
                    // Refer the edge IDs from the TriangleTable by using the table_index.
                    const int e0 = kvs::MarchingHexahedraTable::TriangleID[table_index][t];
                    const int e1 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+2];
                    const int e2 = kvs::MarchingHexahedraTable::TriangleID[table_index][t+1];

                    // Determine vertices for each edge.
                    const int p0 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][0]];
                    const int p1 = local_index[kvs::MarchingHexahedraTable::VertexID[e0][1]];

                    const int p2 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][0]];
                    const int p3 = local_index[kvs::MarchingHexahedraTable::VertexID[e1][1]];

                    const int p4 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][0]];
                    const int p5 = local_index[kvs::MarchingHexahedraTable::VertexID[e2][1]];

                    // Calculate coordinates of the vertices which are composed
                    // of the triangle polygon.
                    const kvs::Vec3 v0( zvolume->coords().data() + p0 * 3 );
                    const kvs::Vec3 v1( zvolume->coords().data() + p1 * 3 );
                    const kvs::Real64 s0 = zvolume->values().at<kvs::Real64>( p0 );
                    const kvs::Real64 s1 = zvolume->values().at<kvs::Real64>( p1 );
                    const kvs::Vec3 vertex0( this->interpolate_vertex( v0, v1, s0, s1 ) );
                    coords.push_back( vertex0.x() );
                    coords.push_back( vertex0.y() );
                    coords.push_back( vertex0.z() );

                    const kvs::Vec3 v2( zvolume->coords().data() + p2 * 3 );
                    const kvs::Vec3 v3( zvolume->coords().data() + p3 * 3 );
                    const kvs::Real64 s2 = zvolume->values().at<kvs::Real64>( p2 );
                    const kvs::Real64 s3 = zvolume->values().at<kvs::Real64>( p3 );
                    const kvs::Vec3 vertex1( this->interpolate_vertex( v2, v3, s2, s3 ) );
                    coords.push_back( vertex1.x() );
                    coords.push_back( vertex1.y() );
                    coords.push_back( vertex1.z() );

                    const kvs::Vec3 v4( zvolume->coords().data() + p4 * 3 );
                    const kvs::Vec3 v5( zvolume->coords().data() + p5 * 3 );
                    const kvs::Real64 s4 = zvolume->values().at<kvs::Real64>( p4 );
                    const kvs::Real64 s5 = zvolume->values().at<kvs::Real64>( p5 );
                    const kvs::Vec3 vertex2( this->interpolate_vertex( v4, v5, s4, s5 ) );
                    coords.push_back( vertex2.x() );
                    coords.push_back( vertex2.y() );
                    coords.push_back( vertex2.z() );

                    // Calculate a normal vector for the triangle polygon.
                    const kvs::Vec3 normal( ( vertex1 - vertex0 ).cross( vertex2 - vertex0 ) );
                    normals.push_back( normal.x() );
                    normals.push_back( normal.y() );
                    normals.push_back( normal.z() );
                } // end of loop-triangle
            }
        }
    }

    // Calculate the polygon color for the isolevel.
    const kvs::RGBColor color = this->calculate_color();

    if ( coords.size() > 0 )
    {
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
    const kvs::Vec3& vertex0,
    const kvs::Vec3& vertex1,
    const double value0,
    const double value1 ) const
{
    const float ratio = static_cast<float>( kvs::Math::Abs( ( m_isolevel - value0 ) / ( value1 - value0 ) ) );
    return ( 1.0f - ratio ) * vertex0 + ratio * vertex1;
}

const kvs::RGBColor Isosurface::calculate_color()
{
    const kvs::Real64 min_value = BaseClass::volume()->minValue();
    const kvs::Real64 max_value = BaseClass::volume()->maxValue();
    const kvs::Real64 normalize_factor = 255.0 / ( max_value - min_value );
    const kvs::UInt8 index = static_cast<kvs::UInt8>( normalize_factor * ( m_isolevel - min_value ) );
    return BaseClass::transferFunction().colorMap()[ index ];
}

} // end of namespace YYZVis
