#include "Program.h"
#include "Input.h"
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
    local::Input input( argc, argv );
    if ( !input.parse() ) { return 1; }

    // Import YYZ data.
    std::cout << "IMPORT VOLUMES ..." << std::endl;
    const std::string input_file( input.filename );
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
    const size_t dim = input.dim;
    auto cart_volume = YYZVis::UniformGridMerger( &yin_volume, &yng_volume, &zng_volume, dim );
    cart_volume.print( std::cout << "STRUCTURED VOLUME DATA" << std::endl, indent );
    cart_volume.write( input.output, ascii );

    return 0;
}

} // end of namespace local
