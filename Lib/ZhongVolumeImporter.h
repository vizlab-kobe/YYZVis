#pragma once
#include "ImporterBase.h"
#include "ZhongVolumeObject.h"


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Importer class for zhong volume object.
 */
/*===========================================================================*/
class ZhongVolumeImporter : public YYZVis::ImporterBase, public YYZVis::ZhongVolumeObject
{
    kvsModule( YYZVis::ZhongVolumeImporter, Importer );
    kvsModuleBaseClass( YYZVis::ImporterBase );
    kvsModuleSuperClass( YYZVis::ZhongVolumeObject );

public:
    ZhongVolumeImporter() {}
    ZhongVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    ZhongVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
