#pragma once
#include <kvs/ImporterBase>
#include <kvs/Module>
#include "YinYangVolumeObject.h"


namespace YYZVis
{

class YangVolumeImporter : public kvs::ImporterBase, public YYZVis::YinYangVolumeObject
{
    kvsModule( YYZVis::YangVolumeImporter, Importer );
    kvsModuleBaseClass( kvs::ImporterBase );
    kvsModuleSuperClass( YYZVis::YinYangVolumeObject );

public:
    YangVolumeImporter() {}
    YangVolumeImporter( const kvs::FileFormatBase* file_format ) { this->exec( file_format ); }
    YangVolumeImporter( const std::string& filename );

    SuperClass* exec( const kvs::FileFormatBase* file_format );
};

} // end of namespace YYZVis
