#include "YinVolumeImporter.h"
#include <kvs/Json>
#include <kvs/File>
#include <iostream>
#include <fstream>


namespace YYZVis
{

/*===========================================================================*/
/**
 *  @brief  Constructs a new YinVolumeImporter class.
 *  @param  filename [in] filename of yin volume data
 */
/*===========================================================================*/
YinVolumeImporter::YinVolumeImporter( const std::string& filename )
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
YinVolumeImporter::SuperClass* YinVolumeImporter::exec( const kvs::FileFormatBase* file_format )
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
    const size_t dim_lat( object["dim_lat"].get<double>() );
    const size_t dim_lon( object["dim_lon"].get<double>() );
    const std::string endian( object["endian"].get<std::string>() );

    auto& yin_value = object["yin_value"].get<kvs::Json::Array>();
    const std::string path_name = kvs::File( json->filename() ).pathName( true );
    const bool swap = BaseClass::needsByteSwap( endian );
    std::vector<kvs::ValueArray<kvs::Real32>> temp;
    for ( auto& f : yin_value )
    {
//        const std::string data_file = f.get<std::string>();
//        const std::string filename =
//            ( BaseClass::isAbsolutePath( data_file ) ) ? data_file :
//            ( data_file[0] == '~' ) ? BaseClass::absolutePath( data_file ) :
//            path_name + kvs::Directory::Separator() + data_file;
        const std::string filename = f.get<std::string>();
        temp.push_back( BaseClass::readBinary<kvs::Real32>( filename, 4, swap ) );
    }

    auto values = BaseClass::interleaveArrays( temp );

//    SuperClass::setGridTypeToYin();
    SuperClass::setDimR( dim_rad );
    SuperClass::setDimTheta( dim_lat );
    SuperClass::setDimPhi( dim_lon );
    SuperClass::setVeclen( temp.size() );
    SuperClass::calculateCoords();
    SuperClass::setValues( values );
    SuperClass::updateMinMaxCoords();
    SuperClass::updateMinMaxValues();


    std::cout << "min = " << SuperClass::minValue() << std::endl;
    std::cout << "max = " << SuperClass::maxValue() << std::endl;

    return this;
}

} // end of namespace YYZVis
