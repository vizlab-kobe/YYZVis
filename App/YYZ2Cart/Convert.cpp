#include "Convert.h"


namespace local
{

kvs::StructuredVolumeObject Convert(
    YinYangVis::YinYangVolumeObject& yin_volume,
    YinYangVis::YinYangVolumeObject& yang_volume,
    YinYangVis::ZhongVolumeObject& zhong_volume )
{
    // Vector length
    //const int veclen = 1;
    
    // Grid resolution (rx, ry, rz)
    //const kvs::Vec3u resolution( rx, ry, rz );

    // Values (physical quantity)
    //kvs::ValueArray<kvs::Real32> values( rx, * ry * rz );

    kvs::StructuredVolumeObject object;
    //object.setGridTypeToUniform();
    //object.setVeclen( veclen );
    //object.setResolution( resolution );
    //object.setValues( values );
    return object;
}

    kvs::ValueArray<kvs::Real32> set_metric( const YinYangVis::YinYangVolumeObject* object )
    {
        const size_t dim_r = object->dimR(); // radius
        const size_t dim_theta = object->dimTheta(); // latitude
        const size_t dim_phi = object->dimPhi(); // longitude
        
        const YinYangVis::YinYangVolumeObject::Range range_r = object->rangeR();
        const YinYangVis::YinYangVolumeObject::Range range_theta = object->rangeTheta();
        const YinYangVis::YinYangVolumeObject::Range range_phi = object->rangePhi();
        
        const kvs::ValueArray<kvs::Real32> ogrid_rad( dim_r );
        const kvs::ValueArray<kvs::Real32> ogrid_theta( dim_theta );
        const kvs::ValueArray<kvs::Real32> ogrid_phi( dim_phi );
        
        for ( int i = 0; i < (int)dim_r; i++ )
        {
            ogrid_rad(i) = range_r.min + range_r.d * i;
        }
        
        for ( int i = 0; i < (int)dim_theta; i++ )
        {
            ogrid_theta(i) = range_theta.min + range_theta.d * ( i - 1 );
        }
        
        for ( int i = 0; i < (int)dim_phi; i++ )
        {
            ogrid_phi(i) = range_phi.min + range_phi.d * ( i - 2 );
        }

        return 0;
    }
    
    kvs::ValueArray<kvs::Real32> ogrid__find_near_corner ( rad,theta,phi,i1,j1,k1,wr1,wt1,wp1 )
    {
        for ( int i = 0; i < (int)dim_r; i++)
        {
            if ( (rad >= ogrid_rad(i)) &&
                (rad <= ogrid_rad(i+1))     )
            {
                i1 = i;
                wr1 = (ogrid_rad(i+1)-rad)/range_r.d;
                break;
            }
        }
        
        for ( int j = 0; j < (int)dim_theta; j++)
        {
            if ( (theta >= ogrid_theta(j)) &&
                (theta <= ogrid_theta(j+1))     )
            {
                j1 = j;
                wt1 = (ogrid_theta(j+1)-theta)/range_theta.d;
                break;
            }
        }
        
        for ( int k = 0; k < (int)dim_phi; k++)
        {
            if ( (phi >= ogrid_phi(k)) &&
                (phi <= ogrid_phi(k+1))     )
            {
                k1 = k;
                wp1 = (ogrid_phi(k+1)-phi)/range_phi.d;
                break;
            }
        }
        return 0;
    }
    
} // end of namespace local
