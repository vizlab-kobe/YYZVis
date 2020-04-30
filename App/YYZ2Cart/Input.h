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
    std::string filename; ///< input filename
    std::string output; ///< output filename
    size_t dim; ///< dimension in the radial direction

    Input( int argc, char** argv ):
        filename(""),
        output("output.kvsml"),
        dim( 200 )
    {
        m_commandline = kvs::CommandLine( argc, argv );
        m_commandline.addOption( "dim", "Grid resolution of output data. (default: 200)", 1, false );
        m_commandline.addOption( "output", "Output filename. (default: output.kvsml)", 1, false );
        m_commandline.addValue( "Input filename." );
        m_commandline.addHelpOption();
    }

    bool parse()
    {
        if ( !m_commandline.parse() ) { return false; }
        if ( m_commandline.hasOption("dim") ) { dim = m_commandline.optionValue<size_t>("dim"); }
        if ( m_commandline.hasOption("output") ) { output = m_commandline.optionValue<std::string>("output"); }
        filename = m_commandline.value<std::string>();
        return true;
    }
};

} // end of namespace local
