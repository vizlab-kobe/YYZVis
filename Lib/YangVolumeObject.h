#pragma once
#include "YinYangVolumeObjectBase.h"


namespace YYZVis
{

class YangVolumeObject : public YYZVis::YinYangVolumeObjectBase
{
    kvsModule( YYZVis::YangVolumeObject, Object );
    kvsModuleBaseClass( YYZVis::YinYangVolumeObjectBase );

public:
    YangVolumeObject() { BaseClass::setGridTypeToYang(); }

private:
    void setGridType( const GridType grid_type );
    void setGridTypeToYin();
    void setGridTypeToYang();
};

} // end of namespace YYZVis
