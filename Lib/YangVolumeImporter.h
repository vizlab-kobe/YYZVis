#pragma once
#include "ImporterBase.h"
//#include "YinYangVolumeObject.h"
#include "YangVolumeObject.h"


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Importer class for yang volume object.
 */
/*===========================================================================*/
//class YangVolumeImporter : public YYZVis::ImporterBase, public YYZVis::YinYangVolumeObject
class YangVolumeImporter : public YYZVis::ImporterBase, public YYZVis::YangVolumeObject
{
    kvsModule( YYZVis::YangVolumeImporter, Importer );
    kvsModuleBaseClass( YYZVis::ImporterBase );
    kvsModuleSuperClass( YYZVis::YangVolumeObject );

public:
    YangVolumeImporter() {}
    YangVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    YangVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
