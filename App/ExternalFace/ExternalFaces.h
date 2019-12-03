#pragma once
#include <kvs/Module>
#include <kvs/MapperBase>
#include <kvs/PolygonObject>
#include <kvs/VolumeObjectBase>
#include <kvs/TransferFunction>
#include <YinYangVis/Lib/YinYangVolumeObject.h>
#include <YinYangVis/Lib/ZhongVolumeObject.h>


namespace YinYangVis
{

class ExternalFaces : public kvs::MapperBase, public kvs::PolygonObject
{
    kvsModule( YinYangVis::ExternalFaces, Mapper );
    kvsModuleBaseClass( kvs::MapperBase );
    kvsModuleSuperClass( kvs::PolygonObject );

public:
    ExternalFaces();
    ExternalFaces( const kvs::VolumeObjectBase* volume );
    ExternalFaces( const kvs::VolumeObjectBase* volume, const kvs::TransferFunction& transfer_function );
    SuperClass* exec( const kvs::ObjectBase* object );

private:
    void mapping( const YinYangVis::ZhongVolumeObject* zvolume );
    void mapping( const YinYangVis::YinYangVolumeObject* yvolume );
};

} // end of namespace YinYangVis
