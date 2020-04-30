#pragma once
#include "YinYangVolumeObjectBase.h"


namespace YYZVis
{

class YinVolumeObject : public YYZVis::YinYangVolumeObjectBase
{
    kvsModule( YYZVis::YinVolumeObject, Object );
    kvsModuleBaseClass( YYZVis::YinYangVolumeObjectBase );

public:
    YinVolumeObject() { BaseClass::setGridTypeToYin(); }

private:
    void setGridType( const GridType grid_type );
    void setGridTypeToYin();
    void setGridTypeToYang();
};

} // end of namespace YYZVis
