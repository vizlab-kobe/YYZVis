#pragma once
#include <kvs/CommandLine>
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
    std::string filename_output; ///< filename of output data
    size_t dim_rad; ///< dimension in the radial direction
    size_t dim_lat; ///< dimension in the latitude direction
    size_t dim_lon; ///< dimension in longitude direction
    size_t dim_zhong; ///< dimension of the cubic zhong grid

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
        m_commandline.addOption( "output", "Filename of output data.", 1, true );
        m_commandline.addHelpOption();
    }

    bool parse()
    {
        if ( !m_commandline.parse() ) { return false; }

        filename_yin = m_commandline.optionValue<std::string>("yin");
        filename_yang = m_commandline.optionValue<std::string>("yang");
        filename_zhong = m_commandline.optionValue<std::string>("zhong");
        filename_output = m_commandline.optionValue<std::string>("output");
        dim_rad = m_commandline.optionValue<size_t>("dim_rad");
        dim_lat = m_commandline.optionValue<size_t>("dim_lat");
        dim_lon = m_commandline.optionValue<size_t>("dim_lon");
        dim_zhong = m_commandline.optionValue<size_t>("dim_zhong");

        return true;
    }
};

} // end of namespace local
