#include "sgrid.h"
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <kvs/ValueArray>
//#include <bits/stdc++.h>
#include <math.h>

namespace local
{
    Sgrid::Sgrid( const YinYangVis::YinYangVolumeObject& object )
    {
        local::Sgrid::Sgrid__size sgrid__size;
        this->set_ogrid( object );
        this->sgrid__make();
        this->sgrid__localize();
    }
    
    void Sgrid::set_ogrid( const YinYangVis::YinYangVolumeObject& object )
    {
        OGRID__THETA_FROM = M_PI/4;
        OGRID__THETA_TO   = M_PI - M_PI/4;
        OGRID__PHI_FROM   = -(3*M_PI)/4.0;
        OGRID__PHI_TO     =  (3*M_PI)/4.0;
        
        ogrid__range_r = object.rangeR();
        ogrid__range_theta = object.rangeTheta();
        ogrid__range_phi = object.rangePhi();
    
    }
    void Sgrid::set_minmax()
    {
        SGRID__THETA_FROM = 0;
        SGRID__THETA_TO   = M_PI;
        SGRID__PHI_FROM   = -M_PI;
        SGRID__PHI_TO     =  M_PI;
        
        sgrid__rad_min = ogrid__range_r.d / 2;
        sgrid__rad_max = ogrid__range_r.max + ogrid__range_r.d;
        
       // Changed by Yoshiki Ueda on 2017.06.23
        //sgrid__tht_min = SGRID__THETA_FROM + ogrid__dtht / 2;
        //sgrid__tht_max = SGRID__THETA_TO - ogrid__dtht / 2;
        sgrid__tht_min = SGRID__THETA_FROM + ogrid__range_theta.d / 10000;
        sgrid__tht_max = SGRID__THETA_TO - ogrid__range_theta.d / 10000;
        
        sgrid__phi_min = SGRID__PHI_FROM + ogrid__range_phi.d / 2;
       // Changed by Yoshiki Ueda on 2017.06.23
        // sgrid__phi_max = SGRID__PHI_TO - ogrid__dphi / 2;
        sgrid__phi_max = SGRID__PHI_TO + ogrid__range_phi.d / 2;

    }
    
    void Sgrid::set_nrtp()
    {
        float drad, dtht, dphi;
        float rad_culling, tht_culling, phi_culling;
        
        rad_culling = 2.5; //namelist__double('DIM1_culling')
        tht_culling = 2.5; //namelist__double('DIM2_culling')
        phi_culling = 2.5; //namelist__double('DIM3_culling')
        
       /* if(rad_culling < 1.0_DP ||
           tht_culling < 1.0_DP ||
           phi_culling < 1.0_DP) then
        call ut__fatal("< sgrid > Error : Culling values must be bigger than 1.0")
        end if*/
        
        drad = rad_culling * ogrid__range_r.d;
        dtht = tht_culling * ogrid__range_theta.d;
        dphi = phi_culling * ogrid__range_phi.d;
        
        sgrid__size.nr = ( sgrid__rad_max - sgrid__rad_min ) / drad + 1;
        sgrid__size.nt = ( sgrid__tht_max - sgrid__tht_min ) / dtht + 1;
        sgrid__size.np = ( sgrid__phi_max - sgrid__phi_min ) / dphi + 1;
        
        if ( sgrid__size.nt%2 == 1 )
        {
            sgrid__size.nt = sgrid__size.nt + 1;
        }

        if ( sgrid__size.np%2 == 1 )
        {
            sgrid__size.np = sgrid__size.np + 1;
        }
      /*  call ut__message('<sgrid> sgrid__size%nr,nt,np = ', &
                         sgrid__size%nr,              &
                         sgrid__size%nt,              &
                         sgrid__size%np)*/
    }
    
    void Sgrid::set_drtp()
    {
        sgrid__drad = ( sgrid__rad_max - sgrid__rad_min ) / sgrid__size.nr - 1;
        sgrid__dtht = ( sgrid__tht_max - sgrid__tht_min ) / sgrid__size.nt - 1;
        sgrid__dphi = ( sgrid__phi_max - sgrid__phi_min ) / sgrid__size.np - 1;
    }
    
    void set_metric()
    {
        int i, j, k;
        float tht, phi;
        size_t nr = sgrid__size.nr;
        kvs::ValueArray<float> sgrid__rad( nr );
        kvs::ValueArray<float> sgrid__theta( sgrid__size.nt );
        for ( i = 0; i < sgrid__size.nr; i++ )
        {
            sgrid__rad[i] = sgrid__rad_min + sgrid__drad * i;
        }
        for ( j = 0; j < sgrid__size.nt; j++ )
        {
            sgrid__theta[j] = sgrid__tht_min + sgrid__dtht * j;
        }
        for ( k = 0; k < sgrid__size.np; k++ )
        {
            sgrid__phi(k) = sgrid__phi_min + sgrid__dphi * k;
        }
        
    }
    
    void Sgrid::sgrid__make()
    {
        this->set_minmax();
        this->set_nrtp();
        this->set_drtp();
        this->set_metric();
    }
  
    void Sgrid::sgrid__localize()
    {
        int j, k;
        float tht, phi;
        
        for ( j = 0; j < sgrid__size.nt; j++)
        {
            tht = sgrid__theta(j);
            sgrid__sintht(j) = sin(tht);
            sgrid__costht(j) = cos(tht);
        }
        
        for ( k = 0; k < sgrid__size.np; k++)
        {
            phi = sgrid__phi(k);
            sgrid__sinphi(k) = sin(phi);
            sgrid__cosphi(k) = cos(phi);
        }
       
    }
}  // end of namespace local
