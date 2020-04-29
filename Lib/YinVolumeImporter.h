#pragma once
#include <kvs/ImporterBase>
#include <kvs/Module>
#include "YinYangVolumeObject.h"


namespace YYZVis
{

class YinVolumeImporter : public kvs::ImporterBase, public YYZVis::YinYangVolumeObject
{
    kvsModule( YYZVis::YinVolumeImporter, Importer );
    kvsModuleBaseClass( kvs::ImporterBase );
    kvsModuleSuperClass( YYZVis::YinYangVolumeObject );

public:
    YinVolumeImporter() {}
    YinVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    YinVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
