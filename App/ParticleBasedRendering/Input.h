#pragma once
#include <kvs/CommandLine>
#include <kvs/TransferFunction>
#include <kvs/DivergingColorMap>
#include <string>


namespace local
{

/*===========================================================================*/
/**
 *  @brief  Input class.
 */
/*===========================================================================*/
class Input
{
private:
    kvs::CommandLine m_commandline; ///< command line parser

public:
    std::string filename_yin; ///< filename of yin volume data
    std::string filename_yang; ///< filename of yang volume data
    std::string filename_zhong; ///< filename of zhong volume data
    size_t repeats; ///< number of repetitions for PBR
    size_t dim_rad; ///< dimension in the radial direction
    size_t dim_lat; ///< dimension in the latitude direction
    size_t dim_lon; ///< dimension in longitude direction
    size_t dim_zhong; ///< dimension of the cubic zhong grid
    bool previous_method; ///< if true, previous PBR method will be used
    bool enable_shading; ///< if true, shading will be available in the PBR
    kvs::TransferFunction tfunc; ///< transfer function

    Input( int argc, char** argv )
    {
        m_commandline = kvs::CommandLine( argc, argv );
        m_commandline.addOption( "dim_rad", "Dimension in the radial direction.", 1, true );
        m_commandline.addOption( "dim_lat", "Dimension in the latitude direction.", 1, true );
        m_commandline.addOption( "dim_lon", "Dimension in the longitude direction.", 1, true );
        m_commandline.addOption( "dim_zhong", "Dimension of the cubic zhong grid.", 1, true );
        m_commandline.addOption( "yin", "Filename of yin volume data.", 1, true );
        m_commandline.addOption( "yang", "Filename of yang volume data.", 1, true );
        m_commandline.addOption( "zhong", "Filename of zhong volume data.", 1, true );
        m_commandline.addOption( "repeat", "Number of repetitions for PBVR. (defulat: 1)", 1, false );
        m_commandline.addOption( "prev", "Use previous particle generation technique." );
        m_commandline.addOption( "shading", "Enable shading." );
        m_commandline.addHelpOption();
    }

    bool parse()
    {
        if ( !m_commandline.parse() ) { return false; }

        // Required.
        filename_yin = m_commandline.optionValue<std::string>("yin");
        filename_yang = m_commandline.optionValue<std::string>("yang");
        filename_zhong = m_commandline.optionValue<std::string>("zhong");
        dim_rad = m_commandline.optionValue<size_t>("dim_rad");
        dim_lat = m_commandline.optionValue<size_t>("dim_lat");
        dim_lon = m_commandline.optionValue<size_t>("dim_lon");
        dim_zhong = m_commandline.optionValue<size_t>("dim_zhong");

        // Optional.
        repeats = m_commandline.hasOption("repeat") ? m_commandline.optionValue<size_t>("repeat") : 1;
        previous_method = m_commandline.hasOption("prev");
        enable_shading = m_commandline.hasOption("shading");
        tfunc = this->create_transfer_function(); // better to be set via a comman line option

        return true;
    }

private:
    kvs::TransferFunction create_transfer_function()
    {
        kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );
        kvs::OpacityMap omap( 256 );
        omap.addPoint(   0, 1.0 );
        omap.addPoint( 150, 0.0 );
        omap.addPoint( 160, 0.0 );
        omap.addPoint( 255, 1.0 );

        /*
          omap.addPoint( 0, 0.9 );
          omap.addPoint( 255, 0.9 );
        */

        /*
          omap.addPoint( 0, 1.0 );
          omap.addPoint( 90, 0.0 );
          omap.addPoint( 180, 0.0 );
          omap.addPoint( 255, 1.0 );
        */

        omap.create();
        return kvs::TransferFunction( cmap, omap );
    }
};

} // end of namespace local
