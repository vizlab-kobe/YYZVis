#include <kvs/UnstructuredVolumeObject>
#include <kvs/ValueArray>
#include <kvs/AnyValueArray>
#include <kvs/Vector3>
#include <kvs/Vector4>
#include <kvs/MarchingHexahedraTable>
#include "SlicePlane.h"

namespace kvs
{
  SlicePlane::SlicePlane(){
    //kvs::MapperBase(),
    kvs::PolygonObject();
  }
  
  SlicePlane::SlicePlane(
			 const kvs::UnstructuredVolumeObject* volume,
			 const kvs::Vector3f point,
			 const kvs::Vector3f normal)
  {
    //SuperClass::setNormalType( normal_type );
    this->setPlane( point, normal);
  
    this->exec( volume );
  }

  SlicePlane::~SlicePlane()
  {
  }
  

  void SlicePlane::setPlane(  const kvs::Vector3f point, const kvs::Vector3f normal )
  {
    coef = kvs::Vector4f( normal, -point.dot( normal ) );
  }
  
  SlicePlane::SuperClass* SlicePlane::exec( const kvs::ObjectBase* object )
  {

    this->extract_plane( object );

    return this;
  }
  
  void SlicePlane::extract_plane( const kvs::ObjectBase* object )
  {
    std::vector<kvs::Real32> coords;
    std::vector<kvs::Real32> normals;
    std::vector<kvs::UInt8>  colors;

    const kvs::VolumeObjectBase* volume_0 = kvs::VolumeObjectBase::DownCast( object );
    BaseClass::attachVolume( volume_0 );
    BaseClass::setRange( volume_0 );
    BaseClass::setMinMaxCoords( volume_0, this );
    const kvs::UnstructuredVolumeObject* volume =  kvs::UnstructuredVolumeObject::DownCast( volume_0 );

    if ( !volume->hasMinMaxValues() )
      {
        volume->updateMinMaxValues();
      }


    const kvs::Real64 min_value( volume->minValue() );
    const kvs::Real64 max_value( volume->maxValue() );
    const kvs::Real64 normalize_factor( 255.0 / ( max_value - min_value ) );
    
    const kvs::ColorMap& color_map( transferFunction().colorMap() );
    
    const kvs::Real32* volume_coords = volume->coords().data();
    const kvs::UInt32 ncells( volume->numberOfCells() );
    const kvs::UInt32* connections = static_cast<const kvs::UInt32*>( volume->connections().data() );


    // Extract surfaces.
    size_t index = 0;
    size_t local_index[8];
    for ( kvs::UInt32 cell = 0; cell < ncells; ++cell, index += 8 )
      {
        // Calculate the indices of the target cell.
        local_index[4] = connections[ index ];
        local_index[5] = connections[ index + 1 ];
        local_index[6] = connections[ index + 2 ];
        local_index[7] = connections[ index + 3 ];
        local_index[0] = connections[ index + 4 ];
        local_index[1] = connections[ index + 5 ];
        local_index[2] = connections[ index + 6 ];
        local_index[3] = connections[ index + 7 ];

        // Calculate the index of the reference table.
        const size_t table_index = this->calculate_table_index( volume ,local_index );

        if ( table_index == 0 ) continue;
        if ( table_index == 255 ) continue;

        // Calculate the triangle polygons.
        for ( size_t i = 0; MarchingHexahedraTable::TriangleID[ table_index ][i] != -1; i += 3 )
	  {
            // Refer the edge IDs from the TriangleTable using the table_index.
            const int e0 = MarchingHexahedraTable::TriangleID[table_index][i];
            const int e1 = MarchingHexahedraTable::TriangleID[table_index][i+1];
            const int e2 = MarchingHexahedraTable::TriangleID[table_index][i+2];

            // Refer indices of the coordinate array from the VertexTable using the edgeIDs.
            const size_t c0 = local_index[ MarchingHexahedraTable::VertexID[e0][0] ];
            const size_t c1 = local_index[ MarchingHexahedraTable::VertexID[e0][1] ];
            const size_t c2 = local_index[ MarchingHexahedraTable::VertexID[e1][0] ];
            const size_t c3 = local_index[ MarchingHexahedraTable::VertexID[e1][1] ];
            const size_t c4 = local_index[ MarchingHexahedraTable::VertexID[e2][0] ];
            const size_t c5 = local_index[ MarchingHexahedraTable::VertexID[e2][1] ];

            // Determine vertices for each edge.
            const kvs::Vector3f v0( volume_coords + 3 * c0 );
            const kvs::Vector3f v1( volume_coords + 3 * c1 );

            const kvs::Vector3f v2( volume_coords + 3 * c2 );
            const kvs::Vector3f v3( volume_coords + 3 * c3 );

            const kvs::Vector3f v4( volume_coords + 3 * c4 );
            const kvs::Vector3f v5( volume_coords + 3 * c5 );

            // Calculate coordinates of the vertices which are composed
            // of the triangle polygon.
            const kvs::Vector3f vertex0( this->interpolate_vertex( v0, v1 ) );
            coords.push_back( vertex0.x() );
            coords.push_back( vertex0.y() );
            coords.push_back( vertex0.z() );

            const kvs::Vector3f vertex1( this->interpolate_vertex( v2, v3 ) );
            coords.push_back( vertex1.x() );
            coords.push_back( vertex1.y() );
            coords.push_back( vertex1.z() );

            const kvs::Vector3f vertex2( this->interpolate_vertex( v4, v5 ) );
            coords.push_back( vertex2.x() );
            coords.push_back( vertex2.y() );
            coords.push_back( vertex2.z() );

            const double value0 = this->interpolate_value( volume, c0, c1 );
            const double value1 = this->interpolate_value( volume, c2, c3 );
            const double value2 = this->interpolate_value( volume, c4, c5 );

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
            const kvs::Vector3f normal( -( vertex2 - vertex0 ).cross( vertex1 - vertex0 ) );
            normals.push_back( normal.x() );
            normals.push_back( normal.y() );
            normals.push_back( normal.z() );

	  } // end of loop-triangle
      } // end of loop-cell

    SuperClass::setCoords( kvs::ValueArray<kvs::Real32>( coords ) );
    SuperClass::setColors( kvs::ValueArray<kvs::UInt8>( colors ) );
    SuperClass::setNormals( kvs::ValueArray<kvs::Real32>( normals ) );
    SuperClass::setOpacity( 255 );
    SuperClass::setPolygonType( kvs::PolygonObject::Triangle );
    SuperClass::setColorType( kvs::PolygonObject::VertexColor );
    SuperClass::setNormalType( kvs::PolygonObject::PolygonNormal );  
  }
  
  kvs::Vector3f SlicePlane::interpolate_vertex(  const kvs::Vector3f& vertex0,
						 const kvs::Vector3f& vertex1 )
  {
    const float value0 = this->plane_function( vertex0 );
    const float value1 = this->plane_function( vertex1 );
    const float ratio = kvs::Math::Abs( value0 / ( value1 - value0 ) );

    return ( 1.0f - ratio ) * vertex0 + ratio * vertex1;
  }
  
  int  SlicePlane::calculate_table_index( const kvs::UnstructuredVolumeObject* volume , size_t local_index[8] )
  {
    kvs::ValueArray<kvs::Real32> coords = volume->coords();

    kvs::Vec3 v0( coords[3 * local_index[0]], coords[3 * local_index[0] + 1], coords[3 * local_index[0] + 2]);
    kvs::Vec3 v1( coords[3 * local_index[1]], coords[3 * local_index[1] + 1], coords[3 * local_index[1] + 2] );
    kvs::Vec3 v2( coords[3 * local_index[2]], coords[3 * local_index[2] + 1], coords[3 * local_index[2] + 2] );
    kvs::Vec3 v3( coords[3 * local_index[3]], coords[3 * local_index[3] + 1], coords[3 * local_index[3] + 2] );
    kvs::Vec3 v4( coords[3 * local_index[4]], coords[3 * local_index[4] + 1], coords[3 * local_index[4] + 2] );
    kvs::Vec3 v5( coords[3 * local_index[5]], coords[3 * local_index[5] + 1], coords[3 * local_index[5] + 2] );
    kvs::Vec3 v6( coords[3 * local_index[6]], coords[3 * local_index[6] + 1], coords[3 * local_index[6] + 2] );
    kvs::Vec3 v7( coords[3 * local_index[7]], coords[3 * local_index[7] + 1], coords[3 * local_index[7] + 2] );
 
    float s0 = plane_function( v0 );
    float s1 = plane_function( v1 );
    float s2 = plane_function( v2 );
    float s3 = plane_function( v3 );
    float s4 = plane_function( v4 );
    float s5 = plane_function( v5 );
    float s6 = plane_function( v6 );
    float s7 = plane_function( v7 );
   
    int index = 0;
    if ( s0 > 0 ) { index |=   1; }
    if ( s1 > 0 ) { index |=   2; }
    if ( s2 > 0 ) { index |=   4; }
    if ( s3 > 0 ) { index |=   8; }
    if ( s4 > 0 ) { index |=  16; }
    if ( s5 > 0 ) { index |=  32; }
    if ( s6 > 0 ) { index |=  64; }
    if ( s7 > 0 ) { index |= 128; }

    return index;
  }
  
  float SlicePlane::plane_function( kvs::Vec3 vertex )
  {
    return coef.x() * vertex.x() + coef.y() * vertex.y() + coef.z() * vertex.z() + coef.w();
  }

  double SlicePlane::interpolate_value(
				       const kvs::UnstructuredVolumeObject* volume,
				       const size_t                         index0,
				       const size_t                         index1 ) 
  {
    const kvs::AnyValueArray values = volume->values();
    const kvs::Real32* const coords = volume->coords().data();

    const float value0 = this->plane_function( kvs::Vector3f( coords + 3 * index0 ) );
    const float value1 = this->plane_function( kvs::Vector3f( coords + 3 * index1 ) );
    const float ratio = kvs::Math::Abs( value0 / ( value1 - value0 ) );

    return values.at<float>( index0 ) + ratio * ( values.at<float>( index1 ) - values.at<float>( index0 ) );
  }
  
}

