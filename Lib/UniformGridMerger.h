#pragma once
#include <kvs/StructuredVolumeObject>
#include <kvs/FilterBase>
#include <kvs/Module>
#include "YinVolumeObject.h"
#include "YangVolumeObject.h"
#include "ZhongVolumeObject.h"


namespace YYZVis
{

class UniformGridMerger : public kvs::FilterBase, public kvs::StructuredVolumeObject
{
    kvsModule( YYZVis::UniformGridMerger, Filter );
    kvsModuleBaseClass( kvs::FilterBase );
    kvsModuleSuperClass( kvs::StructuredVolumeObject );

    typedef YYZVis::YinVolumeObject YinVolumeObject;
    typedef YYZVis::YangVolumeObject YngVolumeObject;
    typedef YYZVis::ZhongVolumeObject ZngVolumeObject;

private:
    const YinVolumeObject* m_yin_volume;
    const YngVolumeObject* m_yng_volume;
    const ZngVolumeObject* m_zng_volume;
    size_t m_dim;

public:
    UniformGridMerger():
        m_yin_volume( NULL ),
        m_yng_volume( NULL ),
        m_zng_volume( NULL ),
        m_dim( 0 ) {}

    UniformGridMerger(
        const YinVolumeObject* yin_volume,
        const YngVolumeObject* yng_volume,
        const ZngVolumeObject* zng_volume,
        const size_t dim ):
        m_yin_volume( yin_volume ),
        m_yng_volume( yng_volume ),
        m_zng_volume( zng_volume ),
        m_dim( dim )
    {
        this->exec( m_yin_volume );
    }

    void setDim( const size_t dim ) { m_dim = dim; }
    void setYinVolumeObject( const YinVolumeObject* yin_volume ) { m_yin_volume = yin_volume; }
    void setYangVolumeObject( const YngVolumeObject* yng_volume ) { m_yng_volume = yng_volume; }
    void setZhongVolumeObject( const ZngVolumeObject* zng_volume ) { m_zng_volume = zng_volume; }

    SuperClass* exec( const kvs::ObjectBase* object );
};

} // end of namespace YYZVis
