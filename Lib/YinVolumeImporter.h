#pragma once
#include "ImporterBase.h"
//#include "YinYangVolumeObject.h"
#include "YinVolumeObject.h"


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Importer class for yin volume object.
 */
/*===========================================================================*/
class YinVolumeImporter : public YYZVis::ImporterBase, public YYZVis::YinVolumeObject
{
    kvsModule( YYZVis::YinVolumeImporter, Importer );
    kvsModuleBaseClass( YYZVis::ImporterBase );
    kvsModuleSuperClass( YYZVis::YinVolumeObject );

public:
    YinVolumeImporter() {}
    YinVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    YinVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
