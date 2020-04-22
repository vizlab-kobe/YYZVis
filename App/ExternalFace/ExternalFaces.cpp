#include "ExternalFaces.h"

#define SET_FACE( coord, v1, v2, v3 )           \
    *( coord++ ) = v1.x();                      \
    *( coord++ ) = v1.y();                      \
    *( coord++ ) = v1.z();                      \
    *( coord++ ) = v2.x();                      \
    *( coord++ ) = v2.y();                      \
    *( coord++ ) = v2.z();                      \
    *( coord++ ) = v3.x();                      \
    *( coord++ ) = v3.y();                      \
    *( coord++ ) = v3.z()

#define SET_NORMAL( normal, n )                 \
    *( normal++ ) = n.x();                      \
    *( normal++ ) = n.y();                      \
    *( normal++ ) = n.z()

#define SET_COLOR( color, c )                  \
    *( color++ ) = c.r();                      \
    *( color++ ) = c.g();                      \
    *( color++ ) = c.b()

namespace
{

/*===========================================================================*/
/**
 *  @brief  Calculates color indices for the given nodes.
 *  @param  values [in] node values
 *  @param  min_value [in] minimum value of the node value
 *  @param  max_value [in] maximum value of the node value
 *  @param  veclen [in] vector length of the node data
 *  @param  colormap_resolution [in] resolution of the color map
 *  @param  node_index [in] node indices
 *  @param  color_level [out] pointer to the color indices
 */
/*===========================================================================*/
template <const size_t N>
inline void GetColorIndices(
    const kvs::AnyValueArray& values,
    const kvs::Real64 min_value,
    const kvs::Real64 max_value,
    const size_t veclen,
    const size_t colormap_resolution,
    const kvs::UInt32 node_index[N],
    kvs::UInt32 (*color_index)[N] )
{
    const kvs::Real64 normalize =
        static_cast<kvs::Real64>( colormap_resolution - 1 ) / ( max_value - min_value );

    // Scalar data.
    if ( veclen == 1 )
    {
        for ( size_t i = 0; i < N; i++ )
        {
            (*color_index)[i] = kvs::UInt32( normalize * ( values.at<kvs::Real64>( node_index[i] ) - min_value ) );
        }
    }
    // Vector data.
    else
    {
        // In case of the vector component, the magnitude value is calculated.
        kvs::Real64 magnitude[N]; memset( magnitude, 0, sizeof( kvs::Real64 ) * N );
        for ( size_t i = 0; i < veclen; ++i )
        {
            for ( size_t j = 0; j < N; j++ )
            {
                magnitude[j] += kvs::Math::Square( values.at<kvs::Real64>( veclen * node_index[j] + i ) );
            }
        }

        for ( size_t i = 0; i < N; i++ )
        {
            magnitude[i] = std::sqrt( magnitude[i] );
            (*color_index)[i] = kvs::UInt32( normalize * ( magnitude[i] - min_value ) );
        }
    }
}

} // end of namespace

namespace local
{

 /*===========================================================================*/
 /**
  *  @brief  Constructs a new ExternalFaces class.
  *  @param  volume [in] pointer to yin, yang, or zhong volume object
  */
 /*===========================================================================*/
ExternalFaces::ExternalFaces( const kvs::VolumeObjectBase* volume ):
    kvs::MapperBase(),
    kvs::PolygonObject()
{
    this->exec( volume );
}

/*===========================================================================*/
/**
 *  @brief  Constructs a new ExternalFaces class.
 *  @param  volume [in] pointer to yin, yang, or zhong volume object
 *  @param  tfunc [in] transfer function
 */
/*===========================================================================*/
ExternalFaces::ExternalFaces(
    const kvs::VolumeObjectBase* volume,
    const kvs::TransferFunction& tfunc ):
    kvs::MapperBase( tfunc ),
    kvs::PolygonObject()
{
    this->exec( volume );
}

/*===========================================================================*/
/**
 *  @brief  Executes external face extraction.
 *  @param  object [in] pointer to yin, yang, or zhong volume object
 *  @return extracted external faces
 */
/*===========================================================================*/
ExternalFaces::SuperClass* ExternalFaces::exec( const kvs::ObjectBase* object )
{
    if ( !object )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input object is NULL." << std::endl;
        return NULL;
    }

    const auto* volume = kvs::VolumeObjectBase::DownCast( object );
    if ( !volume )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input object is not volume dat." << std::endl;
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

/*===========================================================================*/
/**
 *  @brief  Extracts external faces from zhong volume object.
 *  @param  zvolume [in] pointer to zhong volume object
 */
/*===========================================================================*/
void ExternalFaces::mapping( const YYZVis::ZhongVolumeObject* zvolume )
{
    // Calculate coords, colors, and normals.
    this->calculate_coords( zvolume );
    this->calculate_colors( zvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

/*===========================================================================*/
/**
 *  @brief  Extracts external faces from yin or yang volume object.
 *  @param  yvolume [in] pointer to yin or yang volume object
 */
/*===========================================================================*/
void ExternalFaces::mapping( const YYZVis::YinYangVolumeObject* yvolume )
{
    // Calculate coords, colors, and normals.
    this->calculate_yinyang_coords( yvolume );
    this->calculate_yinyang_colors( yvolume );

    // SuperClass::setCoords( coords );
    // SuperClass::setColors( colors );
    // SuperClass::setNormals( normals );
    // if ( SuperClass::numberOfOpacities() == 0 ) SuperClass::setOpacity( 255 );
}

/*===========================================================================*/
/**
 *  @brief  Calculates coordinate values of the specified zhong volume object.
 *  @param  zvolume [in] pointer to zhong volume object
 */
/*===========================================================================*/
void ExternalFaces::calculate_coords( const YYZVis::ZhongVolumeObject* zvolume )
{
    const size_t dim = zvolume->dim(); // number of grid points along each axis
    const size_t nfaces = ( dim - 1 ) * ( dim - 1 ) * 6 * 2; // number of triangle faces
    const size_t nverts = nfaces * 3; // number of vertices

    const float r_min = zvolume->rangeR().min; // min range of R direction
    const kvs::Vec3 min_coord = kvs::Vec3::Constant( -r_min ); // min coord of the volume
    const kvs::Vec3 max_coord = kvs::Vec3::Constant(  r_min ); // max coord of the volume
    const kvs::Vec3 d = max_coord * 2.0f / ( dim - 1.0f ); // distance between grid points

    kvs::ValueArray<kvs::Real32> coords( 3 * nverts );
    kvs::ValueArray<kvs::Real32> normals( 3 * nfaces );
    kvs::Real32* coord = coords.data();
    kvs::Real32* normal = normals.data();

    // XY (Z=Zmin) plane.
    {
        const float z = min_coord.z();
        const kvs::Vec3 n( 0.0f, 0.0f, -1.0f );
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x + d.x(), y, z );
                const kvs::Vec3 v2( x + d.x(), y + d.y(), z );
                const kvs::Vec3 v3( x, y + d.y(), z );
                SET_FACE( coord, v3, v2, v1 );
                SET_FACE( coord, v1, v0, v3 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    // XY (Z=Zmax) plane.
    {
        const float z = max_coord.z();
        const kvs::Vec3 n( 0.0f, 0.0f, 1.0f );
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x + d.x(), y, z );
                const kvs::Vec3 v2( x + d.x(), y + d.y(), z );
                const kvs::Vec3 v3( x, y + d.y(), z );
                SET_FACE( coord, v0, v1, v2 );
                SET_FACE( coord, v2, v3, v0 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    // YZ (X=Xmin) plane.
    {
        const float x = min_coord.x();
        const kvs::Vec3 n( -1.0f, 0.0f, 0.0f );
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t k = 0; k < dim - 1; k++ )
            {
                const float z = min_coord.z() + d.z() * k;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x, y, z + d.z() );
                const kvs::Vec3 v2( x, y + d.y(), z + d.z() );
                const kvs::Vec3 v3( x, y + d.y(), z );
                SET_FACE( coord, v0, v1, v2 );
                SET_FACE( coord, v2, v3, v0 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    // YZ (X=Xmax) plane.
    {
        const float x = max_coord.x();
        const kvs::Vec3 n( 1.0f, 0.0f, 0.0f );
        for ( size_t j = 0; j < dim - 1; j++ )
        {
            const float y = min_coord.y() + d.y() * j;
            for ( size_t k = 0; k < dim - 1; k++ )
            {
                const float z = min_coord.z() + d.z() * k;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x, y, z + d.z() );
                const kvs::Vec3 v2( x, y + d.y(), z + d.z() );
                const kvs::Vec3 v3( x, y + d.y(), z );
                SET_FACE( coord, v3, v2, v1 );
                SET_FACE( coord, v1, v0, v3 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    // XZ (Y=Ymin) plane.
    {
        const float y = min_coord.y();
        const kvs::Vec3 n( 0.0f, -1.0f, 0.0f );
        for ( size_t k = 0; k < dim - 1; k++ )
        {
            const float z = min_coord.z() + d.z() * k;
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x + d.x(), y, z );
                const kvs::Vec3 v2( x + d.x(), y, z + d.z() );
                const kvs::Vec3 v3( x, y, z + d.z() );
                SET_FACE( coord, v0, v1, v2 );
                SET_FACE( coord, v2, v3, v0 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    // XZ (Y=Ymax) plane.
    {
        const float y = max_coord.y();
        const kvs::Vec3 n( 0.0f, 1.0f, 0.0f );
        for ( size_t k = 0; k < dim - 1; k++ )
        {
            const float z = min_coord.z() + d.z() * k;
            for ( size_t i = 0; i < dim - 1; i++ )
            {
                const float x = min_coord.x() + d.x() * i;
                const kvs::Vec3 v0( x, y, z );
                const kvs::Vec3 v1( x + d.x(), y, z );
                const kvs::Vec3 v2( x + d.x(), y, z + d.z() );
                const kvs::Vec3 v3( x, y, z + d.z() );
                SET_FACE( coord, v3, v2, v1 );
                SET_FACE( coord, v1, v0, v3 );
                SET_NORMAL( normal, n );
                SET_NORMAL( normal, n );
            }
        }
    }

    SuperClass::setCoords( coords );
    SuperClass::setNormals( normals );
}

void ExternalFaces::calculate_yinyang_coords( const YYZVis::YinYangVolumeObject* yvolume )
{
    size_t dim_r = 0, dim_theta = 0, dim_phi = 0;

    float theta_from,theta_to,phi_from,phi_to,r_from,r_to;
    float r_d,theta_d,phi_d;
    size_t index = 0;

    r_from = yvolume->rangeR().min;
    r_to =  yvolume->rangeR().max;
    theta_from =  yvolume->rangeTheta().min;
    theta_to =  yvolume->rangeTheta().max;
    phi_from =  yvolume->rangePhi().min;
    phi_to =  yvolume->rangePhi().max;

    r_d = yvolume->rangeR().d;
    theta_d = yvolume->rangeTheta().d;
    phi_d = yvolume->rangePhi().d;

    dim_r = yvolume->dimR();
    dim_theta= yvolume->dimTheta();
    dim_phi = yvolume->dimPhi();

    const size_t nfaces = ((dim_r - 1) * (dim_theta - 1) * 2
                           + (dim_r - 1) * (dim_phi - 1) * 2
                           + (dim_theta - 1) * (dim_phi - 1) * 2) * 2;
    const size_t nverts = nfaces * 3;

    kvs::ValueArray<kvs::Real32> coords( 3 * nverts );
    kvs::Real32* coord = coords.data();
    kvs::ValueArray<kvs::Real32> normals( 3 * nfaces );
    kvs::Real32* normal = normals.data();
    
    // phi = 0	
    {
        const float phi = phi_from;
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );
        for ( size_t j = 0; j < dim_theta - 1; j++ )
	{
            const float theta = theta_from + theta_d * j;
            const float theta_next = theta + theta_d;
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            const float sin_theta_next = std::sin( theta_next );
            const float cos_theta_next = std::cos( theta_next );
            for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
                const float r = r_from + r_d * i;
                const float r_next = r + r_d;
	      
                // v3
                const float x3 = r * sin_theta_next * cos_phi;
                const float y3 = r * sin_theta_next * sin_phi;
                const float z3 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
                // v2
                const float x2 = r_next * sin_theta_next * cos_phi;
                const float y2 = r_next * sin_theta_next * sin_phi;
                const float z2 = r_next * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
                // v1
                const float x1 = r_next * sin_theta * cos_phi;
                const float y1 = r_next * sin_theta * sin_phi;
                const float z1 = r_next * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;

                // v1
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
                // v0
                const float x0 = r * sin_theta * cos_phi;
                const float y0 = r * sin_theta * sin_phi;
                const float z0 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
                // v3
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;

                if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
		{
                    this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
                    index++;
                    this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
                    index++;
		}
                else
		{
                    this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
                    index++;
                    this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
                    index++;
		}
	    }
	}

    }
    
    // phi = dim_phi - 1
    {
        const float phi = phi_to;  //grid_size.z() * ( dim_phi - 1 );
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );
        for ( size_t j = 0; j < dim_theta - 1; j++ )
	{
            const float theta = theta_from + theta_d * j;
            const float theta_next = theta + theta_d;
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            const float sin_theta_next = std::sin( theta_next );
            const float cos_theta_next = std::cos( theta_next );
            for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
                const float r = r_from + r_d * i;
                const float r_next = r + r_d;
	      
                // v0
                const float x0 = r * sin_theta * cos_phi;
                const float y0 = r * sin_theta * sin_phi;
                const float z0 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
                // v1
                const float x1 = r_next * sin_theta * cos_phi;
                const float y1 = r_next * sin_theta * sin_phi;
                const float z1 = r_next * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
                // v2
                const float x2 = r_next * sin_theta_next * cos_phi;
                const float y2 = r_next * sin_theta_next * sin_phi;
                const float z2 = r_next * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;

                // v2
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
                // v3
                const float x3 = r * sin_theta_next * cos_phi;
                const float y3 = r * sin_theta_next * sin_phi;
                const float z3 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
                // v0
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;

                if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
                {
                    this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
                    index++;
                    this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
                    index++;
                }
                else
                {
                    this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
                    index++;
                    this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
                    index++;
                }
	    }
	}
     
        // r = 0
        const float r = r_from;
        for  ( size_t j = 0; j < dim_theta - 1; j++ ) 
        {
            const float theta = theta_from + theta_d * j;
            const float theta_next = theta + theta_d;
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            const float sin_theta_next = std::sin( theta_next );
            const float cos_theta_next = std::cos( theta_next );
            for ( size_t k = 0; k < dim_phi - 1; k++ )
            {
                const float phi = phi_from + phi_d * k;
                const float phi_next = phi + phi_d;
		    
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );
                const float sin_phi_next = std::sin( phi_next );
                const float cos_phi_next = std::cos( phi_next );
	   
                // v0
                const float x0 = r * sin_theta * cos_phi;
                const float y0 = r * sin_theta * sin_phi;
                const float z0 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
                // v1
                const float x1 = r * sin_theta * cos_phi_next;
                const float y1 = r * sin_theta * sin_phi_next;
                const float z1 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
                // v2
                const float x2 = r * sin_theta_next * cos_phi_next;
                const float y2 = r * sin_theta_next * sin_phi_next;
                const float z2 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;

                // v2
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
                // v3
                const float x3 = r * sin_theta_next * cos_phi;
                const float y3 = r * sin_theta_next * sin_phi;
                const float z3 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
                // v0
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;

                if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
                {
                    this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
                    index++;
                    this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
                    index++;
                }
                else
                {
                    this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
                    index++;
                    this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
                    index++;
                }
            }
        }
    }

    
    // r = dim_r - 1
    {
        const float r = r_to; //grid_size.x() * ( dim_r - 1 );
        for ( size_t j = 0; j < dim_theta - 1; j++ )
        {
            const float theta = theta_from + theta_d * j;
            const float theta_next = theta + theta_d;
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            const float sin_theta_next = std::sin( theta_next );
            const float cos_theta_next = std::cos( theta_next );
            for ( size_t k = 0; k < dim_phi - 1; k++ )
            {
                const float phi = phi_from + phi_d * k;
                const float phi_next = phi + phi_d;
		    
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );
                const float sin_phi_next = std::sin( phi_next );
                const float cos_phi_next = std::cos( phi_next );
	     
	        // v3
                const float x3 = r * sin_theta_next * cos_phi;
                const float y3 = r * sin_theta_next * sin_phi;
                const float z3 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
                // v2
                const float x2 = r * sin_theta_next * cos_phi_next;
                const float y2 = r * sin_theta_next * sin_phi_next;
                const float z2 = r * cos_theta_next;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
                // v1
                const float x1 = r * sin_theta * cos_phi_next;
                const float y1 = r * sin_theta * sin_phi_next;
                const float z1 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;

                // v1
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
                // v0
                const float x0 = r * sin_theta * cos_phi;
                const float y0 = r * sin_theta * sin_phi;
                const float z0 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
                // v3
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;

                if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
                {
                    this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
                    index++;
                    this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
                    index++;
                }
                else
                {
                    this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
                    index++;
                    this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
                    index++;
                }
            }
        }
    }

    
    // theta = 0
    {
        const float theta = theta_from;
	const float sin_theta = std::sin( theta );
	const float cos_theta = std::cos( theta );
        for ( size_t k = 0; k < dim_phi - 1; k++ )
	{
	    const float phi = phi_from + phi_d * k;
	    const float phi_next = phi + phi_d;
	    for ( size_t i = 0; i < dim_r - 1; i++ )
	    {
	        const float r = r_from + r_d * i;
		const float r_next = r + r_d;
		    
		const float sin_phi = std::sin( phi );
		const float cos_phi = std::cos( phi );
		const float sin_phi_next = std::sin( phi_next );
		const float cos_phi_next = std::cos( phi_next );
		
		// v0
		const float x0 = r * sin_theta * cos_phi;
		const float y0 = r * sin_theta * sin_phi;
		const float z0 = r * cos_theta;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
		// v1
		const float x1 = r_next * sin_theta * cos_phi;
		const float y1 = r_next * sin_theta * sin_phi;
		const float z1 = r_next * cos_theta;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
		// v2
		const float x2 = r_next * sin_theta * cos_phi_next;
		const float y2 = r_next * sin_theta * sin_phi_next;
		const float z2 = r_next * cos_theta;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;

		// v2
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
		// v3
		const float x3 = r * sin_theta * cos_phi_next;
		const float y3 = r * sin_theta * sin_phi_next;
		const float z3 = r * cos_theta;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
		// v0
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
		*( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;

		if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
                {
		    this->calculate_normal(x0,y0,z0,x1,y1,z1,x2,y2,z2,normal,index);
		    index++;
		    this->calculate_normal(x2,y2,z2,x3,y3,z3,x0,y0,z0,normal,index);
		    index++;
                }
		else
                {
		    this->calculate_normal(-x0,z0,y0,-x1,z1,y1,-x2,z2,y2,normal,index);
		    index++;
		    this->calculate_normal(-x2,z2,y2,-x3,z3,y3,-x0,z0,y0,normal,index);
		    index++;
                }
            }
	}
    }

    // theta = dim_theta - 1
    {
        const float theta = theta_to; //grid_size.y() * ( dim_theta - 1 );
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );
        for ( size_t k = 0; k < dim_phi - 1; k++ )
        {
            const float phi = phi_from + phi_d * k;
            const float phi_next = phi + phi_d;
            for ( size_t i = 0; i < dim_r - 1; i++ )
            {
                const float r = r_from + r_d * i;
                const float r_next = r + r_d;
		    
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );
                const float sin_phi_next = std::sin( phi_next );
                const float cos_phi_next = std::cos( phi_next );
		    
                // v3
                const float x3 = r * sin_theta * cos_phi_next;
                const float y3 = r * sin_theta * sin_phi_next;
                const float z3 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;
                // v2
                const float x2 = r_next * sin_theta * cos_phi_next;
                const float y2 = r_next * sin_theta * sin_phi_next;
                const float z2 = r_next * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x2 : -x2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y2 : z2;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z2 : y2;
                // v1
                const float x1 = r_next * sin_theta * cos_phi;
                const float y1 = r_next * sin_theta * sin_phi;
                const float z1 = r_next * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;

                // v1
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x1 : -x1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y1 : z1;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z1 : y1;
                // v0
                const float x0 = r * sin_theta * cos_phi;
                const float y0 = r * sin_theta * sin_phi;
                const float z0 = r * cos_theta;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x0 : -x0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y0 : z0;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z0 : y0;
                // v3
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x3 : -x3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y3 : z3;
                *( coord++) = ( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z3 : y3;

                if( yvolume->gridType() == YYZVis::YinYangVolumeObject::Yin )
                {
                    this->calculate_normal(x3,y3,z3,x2,y2,z2,x1,y1,z1,normal,index);
                    index++;
                    this->calculate_normal(x1,y1,z1,x0,y0,z0,x3,y3,z3,normal,index);
                    index++;
                }
                else
                {
                    this->calculate_normal(-x3,z3,y3,-x2,z2,y2,-x1,z1,y1,normal,index);
                    index++;
                    this->calculate_normal(-x1,z1,y1,-x0,z0,y0,-x3,z3,y3,normal,index);
                    index++;
                }
            }
        }
    }

    SuperClass::setCoords( coords );
    SuperClass::setNormals( normals );
}

void ExternalFaces::calculate_yinyang_colors( const YYZVis::YinYangVolumeObject* yvolume
    )
{
    size_t dim_r = 0, dim_theta = 0, dim_phi = 0;
    const kvs::Real64 min_value = yvolume->minValue();
    const kvs::Real64 max_value = yvolume->maxValue();
    kvs::AnyValueArray value =  yvolume->values();

    kvs::UInt32 node_index[4];
    kvs::UInt32 color_level[4];
    const kvs::ColorMap cmap( BaseClass::colorMap() );

    dim_r = yvolume->dimR();
    dim_theta= yvolume->dimTheta();
    dim_phi = yvolume->dimPhi();

    float theta_from,theta_to,phi_from,phi_to,r_from,r_to;
    float r_d,theta_d,phi_d;

    r_from = yvolume->rangeR().min;
    r_to =  yvolume->rangeR().max;
    theta_from =  yvolume->rangeTheta().min;
    theta_to =  yvolume->rangeTheta().max;
    phi_from =  yvolume->rangePhi().min;
    phi_to =  yvolume->rangePhi().max;

    r_d = yvolume->rangeR().d;
    theta_d = yvolume->rangeTheta().d;
    phi_d = yvolume->rangePhi().d;

    const size_t nfaces = ((dim_r - 1) * (dim_theta - 1) * 2
                           + (dim_r - 1) * (dim_phi - 1) * 2
                           + (dim_theta - 1) * (dim_phi - 1) * 2) * 2;
    const size_t nverts = nfaces * 3;

    const size_t veclen = yvolume->veclen();

    kvs::ValueArray<kvs::UInt8> colors( 3 * nverts );
    kvs::UInt8* color = colors.data();

    // phi = 0 
    {
        const size_t k = 0;
	const size_t offset0 = k * dim_r;
        for ( size_t j = 0, offset = offset0 ; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset++ )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(),node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );


                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }

    // phi = dim_phi - 1
    {
        const size_t k = dim_phi - 1;
        const size_t offset0 = k * dim_r * dim_theta;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // r = 0
    {
        const size_t i = 0;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1 ; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t k = 0; k < dim_phi - 1; k++, offset += dim_r * dim_theta )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + ( dim_r * dim_theta );
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // r = dim_r - 1
    {
        const size_t i = dim_r - 1;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim_theta - 1; j++, offset = offset0 + j * dim_r )
        {
            for ( size_t k = 0; k < dim_phi - 1; k++, offset += ( dim_r * dim_theta ) )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + ( dim_r * dim_theta );
                node_index[2] = node_index[1] + dim_r;
                node_index[3] = node_index[0] + dim_r;
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }

    // theta = 0
    {
        const size_t j = 0;
        const size_t offset0 = j * dim_r;
        for ( size_t k = 0, offset = offset0; k < dim_phi - 1; k++, offset =  offset0 + k * ( dim_r * dim_theta ) )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + ( dim_r * dim_theta );
                node_index[3] = node_index[0] + ( dim_r * dim_theta );
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();

                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
            }
        }
    }

    // theta = dim_theta - 1
    {
        const size_t j = dim_theta;
        const size_t offset0 = j * dim_r;
        for ( size_t k = 0, offset = offset0; k < dim_phi - 1; k++, offset = offset0 + k * ( dim_r * dim_theta ) )
        {
            for ( size_t i = 0; i < dim_r - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + ( dim_r * dim_theta );
                node_index[3] = node_index[0] + ( dim_r * dim_theta );
//                this->GetColorIndices( value, min_value, max_value, cmap.resolution(), node_index, &color_level );
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
                // v2
                *( color++ ) = cmap[ color_level[2] ].r();
                *( color++ ) = cmap[ color_level[2] ].g();
                *( color++ ) = cmap[ color_level[2] ].b();
                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();

                // v1
                *( color++ ) = cmap[ color_level[1] ].r();
                *( color++ ) = cmap[ color_level[1] ].g();
                *( color++ ) = cmap[ color_level[1] ].b();
                // v0
                *( color++ ) = cmap[ color_level[0] ].r();
                *( color++ ) = cmap[ color_level[0] ].g();
                *( color++ ) = cmap[ color_level[0] ].b();
                // v3
                *( color++ ) = cmap[ color_level[3] ].r();
                *( color++ ) = cmap[ color_level[3] ].g();
                *( color++ ) = cmap[ color_level[3] ].b();
            }
        }
    }
    SuperClass::setColors( colors );
}

void ExternalFaces::calculate_colors( const YYZVis::ZhongVolumeObject* zvolume )
{
    // Parameters of the volume data.
    const kvs::Real64 min_value = zvolume->minValue();
    const kvs::Real64 max_value = zvolume->maxValue();
    const kvs::AnyValueArray& value = zvolume->values();
    const size_t veclen = zvolume->veclen();
    const size_t dim = zvolume->dim();
    const size_t nnodes_per_line = dim;
    const size_t nnodes_per_slice = dim * dim;
    const size_t nfaces = (dim - 1) * (dim - 1) * 6 * 2;
    const size_t nverts = nfaces * 3;

    const kvs::ColorMap cmap( BaseClass::colorMap() );

    kvs::ValueArray<kvs::UInt8> colors( 3 * nverts );
    kvs::UInt8* color = colors.data();

    kvs::UInt32 node_index[4];
    kvs::UInt32 color_level[4];

    // XY (Z=Zmin) plane.
    {
        const size_t k = 0;
        const size_t offset0 = k * nnodes_per_line;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
        {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v3-v2-v1
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                // v1-v0-v3
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
            }
        }
    }

    // XY (Z=Zmax) plane.
    {
        const size_t k = dim - 1;
        const size_t offset0 = k * nnodes_per_slice;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
        {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0-v1-v2
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                // v2-v3-v0
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
            }
        }
    }

    // YZ (X=Xmin) plane.
    {
        const size_t i = 0;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
        {
            for ( size_t k = 0; k < dim - 1; k++, offset += nnodes_per_slice )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + nnodes_per_slice;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0-v1-v2
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                // v2-v3-v0
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
            }
        }
    }

    // YZ (X=Xmax) plane.
    {
        const size_t i = dim - 1;
        const size_t offset0 = i;
        for ( size_t j = 0, offset = offset0; j < dim - 1; j++, offset = offset0 + j * nnodes_per_line )
        {
            for ( size_t k = 0; k < dim - 1; k++, offset += nnodes_per_slice )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + nnodes_per_slice;
                node_index[2] = node_index[1] + nnodes_per_line;
                node_index[3] = node_index[0] + nnodes_per_line;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v3-v2-v1
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                // v1-v0-v3
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
            }
        }
    }

    // XZ (Y=Ymin) plane.
    {
        const size_t j = 0;
        const size_t offset0 = j * nnodes_per_line;
        for ( size_t k = 0, offset = offset0; k < dim - 1; k++, offset =  offset0 + k * nnodes_per_slice )
        {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_slice;
                node_index[3] = node_index[0] + nnodes_per_slice;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v0-v1-v2
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                // v2-v3-v0
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
            }
        }
    }

    // XZ (Y=Ymax) plane.
    {
        const size_t j = dim - 1;
        const size_t offset0 = j * nnodes_per_line;
        for ( size_t k = 0, offset = offset0; k < dim - 1; k++, offset = offset0 + k * nnodes_per_slice )
        {
            for ( size_t i = 0; i < dim - 1; i++, offset += 1 )
            {
                node_index[0] = offset;
                node_index[1] = node_index[0] + 1;
                node_index[2] = node_index[1] + nnodes_per_slice;
                node_index[3] = node_index[0] + nnodes_per_slice;
                ::GetColorIndices<4>( value, min_value, max_value, veclen, cmap.resolution(), node_index, &color_level );
                // v3-v2-v1
                SET_COLOR( color, cmap[ color_level[3] ] );
                SET_COLOR( color, cmap[ color_level[2] ] );
                SET_COLOR( color, cmap[ color_level[1] ] );
                // v1-v0-v3
                SET_COLOR( color, cmap[ color_level[1] ] );
                SET_COLOR( color, cmap[ color_level[0] ] );
                SET_COLOR( color, cmap[ color_level[3] ] );
            }
        }
    }

    SuperClass::setColors( colors );
}

void ExternalFaces::calculate_normal( const float x0, const float y0, const float z0,
                                      const float x1, const float y1, const float z1,
				      const float x2, const float y2, const float z2, kvs::Real32* normal, size_t index )
{
    float x4=0,y4=0,z4=0;
    float x5=0,y5=0,z5=0;
    float x6=0,y6=0,z6=0;
    float er=0;

    x4 = x1 - x0;
    y4 = y1 - y0;
    z4 = z1 - z0;

    x5 = x2 - x0;
    y5 = y2 - y0;
    z5 = z2 - z0;

    x6 = y4 * z5 - y5 * z4;
    y6 = z4 * x5 - x4 * z5;
    z6 = x4 * y5 - x5 * y4;

    er = sqrt ( x6 * x6 + y6 * y6 + z6 * z6 );
    x6 = x6 / er;
    y6 = y6 / er;
    z6 = z6 / er;

    normal[ 3 * index ] = x6;
    normal[ 3 * index + 1 ] = y6;
    normal[ 3 * index + 2 ] = z6;
}

} // end of namespace local
