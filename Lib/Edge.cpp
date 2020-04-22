#include "Edge.h"
#include "YinYangVolumeObject.h"
#include "ZhongVolumeObject.h"
#include <kvs/Xform>
#include <kvs/XformControl>
#include <kvs/LineObject>


namespace YYZVis
{

namespace Edge
{

kvs::LineObject* CreateLineMeshObject( const YYZVis::YinYangVolumeObject* volume, const size_t dim_edge )
{
    const size_t dim_r = volume->dimR(); // radius
    const size_t dim_theta = volume->dimTheta(); // latitude
    const size_t dim_phi = volume->dimPhi(); // longitude

    const YYZVis::YinYangVolumeObject::Range range_r = volume->rangeR();
    const YYZVis::YinYangVolumeObject::Range range_theta = volume->rangeTheta();
    const YYZVis::YinYangVolumeObject::Range range_phi = volume->rangePhi();

    float step_r = ( dim_r - 1.0f ) / dim_edge;
    float step_theta = ( dim_theta - 1.0f ) / dim_edge;
    float step_phi = ( dim_phi - 1.0f ) / dim_edge;

    const size_t nnodes = dim_theta * ( dim_edge * 2 + 2 ) + dim_phi * ( dim_edge * 2 + 2 ) + dim_theta * ( dim_edge * 2 - 2 ) + dim_phi * ( dim_edge * 2 - 2 );
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();

    for ( int k = 0; k < dim_phi; k += dim_phi - 1 )
    {
        const float phi = range_phi.min + range_phi.d * ( k - 2 );
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );
        for ( int j = 0, nstep = 1; j < dim_r; j = step_r * nstep, nstep++ )
        {
            const float r = range_r.min + range_r.d * j;
            for ( int i = 0; i < dim_theta; i++ )
            {
                const float theta = range_theta.min + range_theta.d * ( i - 1 );
                const float sin_theta = std::sin( theta );
                const float cos_theta = std::cos( theta );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    for ( int k = 0; k < dim_theta; k += dim_theta - 1 )
    {
        const float theta = range_theta.min + range_theta.d * ( k - 1 );
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );
        for ( int j = 0, nstep = 1; j < dim_r; j = step_r * nstep, nstep++ )
        {
            const float r = range_r.min + range_r.d * j;
            for ( int i = 0; i < dim_phi; i++ )
            {
                const float phi = range_phi.min + range_phi.d * ( i - 2 );
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    for ( int k = 0; k < dim_r; k += dim_r - 1 )
    {
        const float r = range_r.min + range_r.d * k;
        for ( int j = step_phi, nstep = 2; j < dim_phi - step_phi + 1; j = step_phi * nstep, nstep++ )
        {
            const float phi = range_phi.min + range_phi.d * ( j - 2 );
            const float sin_phi = std::sin( phi );
            const float cos_phi = std::cos( phi );
            for ( int i = 0; i < dim_theta; i++ )
            {
                const float theta = range_theta.min + range_theta.d * ( i - 1 );
                const float sin_theta = std::sin( theta );
                const float cos_theta = std::cos( theta );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    for ( int k = 0; k < dim_r; k += dim_r - 1 )
    {
        const float r = range_r.min + range_r.d * k;
        for ( int j = step_theta, nstep = 2; j < dim_theta - step_theta + 1; j = step_theta * nstep, nstep++ )
        {
            const float theta = range_theta.min + range_theta.d * ( j - 1 );
            const float sin_theta = std::sin( theta );
            const float cos_theta = std::cos( theta );
            for ( int i = 0; i < dim_phi; i++ )
            {
                const float phi = range_phi.min + range_phi.d * ( i - 2 );
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    const size_t nconnections = ( dim_theta - 1 ) * ( dim_edge * 2 + 2 ) + ( dim_phi - 1 ) * ( dim_edge * 2 + 2 ) + ( dim_theta - 1 ) * ( dim_edge * 2 - 2 ) + ( dim_phi - 1 ) * ( dim_edge * 2 - 2 ) + ( dim_edge * 4 + 4 );
    kvs::ValueArray<kvs::UInt32> connections( nconnections * 2 );
    kvs::UInt32* pconnections = connections.data();

    size_t index;
    size_t index2;

    for ( int j = 0; j < dim_theta * ( dim_edge * 2 + 1 ) + 1; j += dim_theta, index = j )
    {
        for ( int i = 0; i < dim_theta - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    index2 = index;
    for ( int j = index; j < dim_phi * ( dim_edge * 2 + 1 ) + 1 + index2; j += dim_phi, index = j )
    {
        for ( int i = 0; i < dim_phi - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    index2 = index;
    for ( int j = index; j < dim_theta * ( dim_edge * 2 - 3 ) + 1 + index2; j += dim_theta, index = j )
    {
        for ( int i = 0; i < dim_theta - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    index2 = index;
    for ( int j = index; j < dim_phi * ( dim_edge * 2 - 3 ) + 1 + index2; j += dim_phi, index = j )
    {
        for ( int i = 0; i < dim_phi - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    for ( int j = 0; j < 2 * dim_theta * ( dim_edge + 1 ); j += dim_theta * ( dim_edge + 1 ), index = j )
    {
        for ( int i = 0; i < dim_edge + 1; i++ )
        {
            *(pconnections++) = step_theta * i + j;
            *(pconnections++) = dim_theta * dim_edge + step_theta * i + j;
        }
    }

    index2 = index;
    for ( int j = index; j < 2 * dim_phi * ( dim_edge + 1 ) + index2; j += dim_phi * ( dim_edge + 1 ) )
    {
        for ( int i = 0; i < dim_edge + 1; i++ )
        {
            *(pconnections++) = step_phi * i + j;
            *(pconnections++) = dim_phi * dim_edge + step_phi * i + j;
        }
    }

    kvs::ValueArray<kvs::UInt8> colors( nnodes * 3 );
    kvs::UInt8* pcolors = colors.data();
    for ( int i = 0; i < nnodes; i++ )
    {
        *(pcolors++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? 255 : 0;
        *(pcolors++) = 0;
        *(pcolors++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? 0 : 255;
    }

    kvs::LineObject* object = new kvs::LineObject();
    object->setCoords( coords );
    object->setColors( colors );
    object->setConnections( connections );
    object->setLineType( kvs::LineObject::Segment );
    object->setColorType( kvs::LineObject::VertexColor );

//  kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
//  object->multiplyXform( x );

    return object;
}

kvs::LineObject* CreateLineEdgeObject( const YYZVis::YinYangVolumeObject* volume )
{
    const size_t dim_r = volume->dimR(); // radius
    const size_t dim_theta = volume->dimTheta(); // latitude
    const size_t dim_phi = volume->dimPhi(); // longitude

    const YYZVis::YinYangVolumeObject::Range range_r = volume->rangeR();
    const YYZVis::YinYangVolumeObject::Range range_theta = volume->rangeTheta();
    const YYZVis::YinYangVolumeObject::Range range_phi = volume->rangePhi();

    const size_t nnodes = dim_theta * 4 + dim_phi * 4;
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();

    for ( int k = 0; k < dim_phi; k += dim_phi - 1 )
    {
        const float phi = range_phi.min + range_phi.d * ( k - 2 );
        const float sin_phi = std::sin( phi );
        const float cos_phi = std::cos( phi );
        for ( int j = 0; j < dim_r; j += dim_r - 1 )
        {
            const float r = range_r.min + range_r.d * j;
            for ( int i = 0; i < dim_theta; i++ )
            {
                const float theta = range_theta.min + range_theta.d * ( i - 1 );
                const float sin_theta = std::sin( theta );
                const float cos_theta = std::cos( theta );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }
    for ( int k = 0; k < dim_theta; k += dim_theta - 1 )
    {
        const float theta = range_theta.min + range_theta.d * ( k - 1 );
        const float sin_theta = std::sin( theta );
        const float cos_theta = std::cos( theta );
        for ( int j = 0; j < dim_r; j += dim_r - 1 )
        {
            const float r = range_r.min + range_r.d * j;
            for ( int i = 0; i < dim_phi; i++ )
            {
                const float phi = range_phi.min + range_phi.d * ( i - 2 );
                const float sin_phi = std::sin( phi );
                const float cos_phi = std::cos( phi );

                const float x = r * sin_theta * cos_phi;
                const float y = r * sin_theta * sin_phi;
                const float z = r * cos_theta;

                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? x : -x;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? y : z;
                *(pcoords++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? z : y;
            }
        }
    }

    const size_t nconnections = ( dim_theta - 1 ) * 4 + ( dim_phi - 1 ) * 4 + 4;
    kvs::ValueArray<kvs::UInt32> connections( nconnections * 2 );
    kvs::UInt32* pconnections = connections.data();

    size_t index;
    for ( int j = 0; j < dim_theta * 3 + 1; j += dim_theta, index = j )
    {
        for ( int i = 0; i < dim_theta - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    for ( int j = index; j < dim_phi * 3 + 1 + index; j += dim_phi )
    {
        for ( int i = 0; i < dim_phi - 1; i++ )
        {
            *(pconnections++) = i + j;
            *(pconnections++) = (i + 1) + j;
        }
    }

    *(pconnections++) = 0;
    *(pconnections++) = dim_theta;
    *(pconnections++) = dim_theta - 1;
    *(pconnections++) = dim_theta * 2 - 1;
    *(pconnections++) = dim_theta * 2;
    *(pconnections++) = dim_theta * 3;
    *(pconnections++) = dim_theta * 3 - 1;
    *(pconnections++) = dim_theta * 4 - 1;

    kvs::ValueArray<kvs::UInt8> colors( nnodes * 3 );
    kvs::UInt8* pcolors = colors.data();
    for ( int i = 0; i < nnodes; i++ )
    {
        *(pcolors++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? 255 : 0;
        *(pcolors++) = 0;
        *(pcolors++) = ( volume->gridType() == YYZVis::YinYangVolumeObject::Yin ) ? 0 : 255;
    }

    kvs::LineObject* object = new kvs::LineObject();
    object->setCoords( coords );
    object->setColors( colors );
    object->setConnections( connections );
    object->setLineType( kvs::LineObject::Segment );
    object->setColorType( kvs::LineObject::VertexColor );

//    kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
//    object->multiplyXform( x );

    return object;
}


kvs::LineObject* CreateLineEdgeObject( const YYZVis::ZhongVolumeObject* volume )
{
    const size_t dim = volume->dim();
    const float r_min = volume->rangeR().min;
    const float r_d = volume->rangeR().d;

    //ix(= iy, iz),dix(= diy, diz)
    const float r_i = r_min + r_d * 2;
    const float ix_max = r_i;
    const float ix_min = -r_i;
    const float dix = ( ix_max - ix_min ) / ( dim - 1 );

    const size_t nnodes = volume->numberOfNodes();
    kvs::ValueArray<kvs::Real32> coords( nnodes * 3 );
    kvs::Real32* pcoords = coords.data();
    for ( size_t k = 0; k < dim; k += dim - 1 )
    {
        const float z = ix_min + dix * k;
        for ( size_t j = 0; j < dim; j += dim - 1 )
        {
            const float y = ix_min + dix * j;
            for ( size_t i = 0; i < dim; i += dim - 1 )
            {
                const float x = ix_min + dix * i;
                *(pcoords++) = x;
                *(pcoords++) = y;
                *(pcoords++) = z;
            }
        }
    }

    const size_t nconnections = 12;
    kvs::ValueArray<kvs::UInt32> connections( nconnections * 2 );
    kvs::UInt32* pconnections = connections.data();

    //size_t index;
    *(pconnections++) = 0;
    *(pconnections++) = 1;
    *(pconnections++) = 0;
    *(pconnections++) = 2;
    *(pconnections++) = 0;
    *(pconnections++) = 4;
    *(pconnections++) = 1;
    *(pconnections++) = 3;
    *(pconnections++) = 1;
    *(pconnections++) = 5;
    *(pconnections++) = 2;
    *(pconnections++) = 3;
    *(pconnections++) = 2;
    *(pconnections++) = 6;
    *(pconnections++) = 3;
    *(pconnections++) = 7;
    *(pconnections++) = 4;
    *(pconnections++) = 5;
    *(pconnections++) = 4;
    *(pconnections++) = 6;
    *(pconnections++) = 5;
    *(pconnections++) = 7;
    *(pconnections++) = 6;
    *(pconnections++) = 7;

    kvs::ValueArray<kvs::UInt8> colors( nnodes * 3 );
    kvs::UInt8* pcolors = colors.data();
    for ( int i = 0; i < nnodes; i++ )
    {
        *(pcolors++) = 0;
        *(pcolors++) = 0;
        *(pcolors++) = 0;
    }

    kvs::LineObject* object = new kvs::LineObject();
    object->setCoords( coords );
    object->setColors( colors );
    object->setConnections( connections );
    object->setLineType( kvs::LineObject::Segment );
    object->setColorType( kvs::LineObject::VertexColor );

//    kvs::Xform x = kvs::Xform::Rotation( kvs::Mat3::RotationX(-135) );
//    object->multiplyXform( x );

    return object;
}

}

} // end of namespace YYZVis
