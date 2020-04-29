#pragma once
#include <kvs/ImporterBase>
#include <kvs/Module>
#include "ZhongVolumeObject.h"


namespace YYZVis
{

class ZhongVolumeImporter : public kvs::ImporterBase, public YYZVis::ZhongVolumeObject
{
    kvsModule( YYZVis::ZhongVolumeImporter, Importer );
    kvsModuleBaseClass( kvs::ImporterBase );
    kvsModuleSuperClass( YYZVis::ZhongVolumeObject );

public:
    ZhongVolumeImporter() {}
    ZhongVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    ZhongVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
