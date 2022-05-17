#include "ZhongVolumeImporter.h"
#include <kvs/Json>
#include <kvs/File>
#include <iostream>
#include <fstream>


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Constructs a new ZhongVolumeImporter class.
 *  @param  filename [in] filename of zhong volume data
 */
/*===========================================================================*/
ZhongVolumeImporter::ZhongVolumeImporter( const std::string& filename )
{
    kvs::Json json( filename );
    this->exec( &json );
}

/*===========================================================================*/
/**
 *  @brief  Executes the importing process.
 *  @param  file_format [in] pointer to the data file
 *  @return true if the importing process was done successfully
 */
/*===========================================================================*/
ZhongVolumeImporter::SuperClass* ZhongVolumeImporter::exec( const kvs::FileFormatBase* file_format )
{
    if ( !file_format )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input file is NULL." << std::endl;
        return NULL;
    }

    if ( !dynamic_cast<const kvs::Json*>( file_format ) )
    {
        BaseClass::setSuccess( false );
        kvsMessageError() << "Input file is not JSON file." << std::endl;
        return NULL;
    }

    kvs::Json* json = (kvs::Json*)( file_format );
    auto& object = json->rootObject();

    const size_t dim_rad( object["dim_rad"].get<double>() );
    const size_t dim_zhong( object["dim_zhong"].get<double>() );
    const std::string endian( object["endian"].get<std::string>() );

    auto& zhong_value = object["zhong_value"].get<kvs::Json::Array>();
    const std::string path_name = kvs::File( json->filename() ).pathName( true );
    const bool swap = BaseClass::needsByteSwap( endian );
    std::vector<kvs::ValueArray<kvs::Real32>> temp;
    for ( auto& f : zhong_value )
    {
//        const std::string data_file = f.get<std::string>();
//        const std::string filename =
//            ( BaseClass::isAbsolutePath( data_file ) ) ? data_file :
//            ( data_file[0] == '~' ) ? BaseClass::absolutePath( data_file ) :
//            path_name + kvs::Directory::Separator() + data_file;
        const auto filename = f.get<std::string>();
        temp.push_back( BaseClass::readBinary<kvs::Real32>( filename, 4, swap ) );
    }

    auto values = BaseClass::interleaveArrays( temp );

    SuperClass::setDimR( dim_rad );
    SuperClass::setDim( dim_zhong );
    SuperClass::setVeclen( temp.size() );
    SuperClass::calculateCoords();
    SuperClass::setValues( values );
    SuperClass::updateMinMaxCoords();
    SuperClass::updateMinMaxValues();
    return this;
}

} // end of namespace YYZVis
