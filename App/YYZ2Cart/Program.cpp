#include "Program.h"
#include "Input.h"
#include "Convert.h"
#include <YYZVis/Lib/YinVolumeImporter.h>
#include <YYZVis/Lib/YangVolumeImporter.h>
#include <YYZVis/Lib/ZhongVolumeImporter.h>
#include <YYZVis/Lib/UpdateMinMaxValues.h>
#include <YYZVis/Lib/UpdateMinMaxCoords.h>
#include <YYZVis/Lib/UniformGridMerger.h>
#include <kvs/Indent>


namespace local
{

int Program::exec( int argc, char** argv )
{
    // Import YYZ data.
    std::cout << "IMPORT VOLUMES ..." << std::endl;
    const std::string input_file( argv[1] );
    auto yin_volume = YYZVis::YinVolumeImporter( input_file );
    auto yng_volume = YYZVis::YangVolumeImporter( input_file );
    auto zng_volume = YYZVis::ZhongVolumeImporter( input_file );
    YYZVis::UpdateMinMaxValues( &yin_volume, &yng_volume, &zng_volume );
    YYZVis::UpdateMinMaxCoords( &yin_volume, &yng_volume, &zng_volume );

    // Dump.
    const kvs::Indent indent( 4 );
    yin_volume.print( std::cout << "YIN VOLUME DATA" << std::endl, indent );
    yng_volume.print( std::cout << "YANG VOLUME DATA" << std::endl, indent );
    zng_volume.print( std::cout << "ZHONG VOLUME DATA" << std::endl, indent );

    std::cout << "CONVERT VOLUMES ..." << std::endl;
    const bool ascii = false;
//    kvs::StructuredVolumeObject cart_volume = local::Convert( yin_volume, yng_volume, zng_volume );
    const size_t dim = 200;
    auto cart_volume = YYZVis::UniformGridMerger( &yin_volume, &yng_volume, &zng_volume, dim );
    cart_volume.print( std::cout << "STRUCTURED VOLUME DATA" << std::endl, indent );
    cart_volume.write( "output.kvsml", ascii );

    return 0;
}

} // end of namespace local
