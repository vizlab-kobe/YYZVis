#include "ZhongVolumeImporter.h"
#include <kvs/Json>
#include <kvs/File>
#include <kvs/Directory>
#include <kvs/Endian>
#include <iostream>
#include <fstream>


namespace
{

inline bool ByteSwap( const std::string& endian )
{
    if ( endian == std::string("little") )
    {
        return !kvs::Endian::IsLittle();
    }
    else if ( endian == std::string("big") )
    {
        return !kvs::Endian::IsBig();
    }
    return false;
}

inline kvs::ValueArray<float> ReadBinary( const std::string& filename, const size_t offset, const bool swap )
{
    std::ifstream ifs( filename.c_str(), std::ifstream::binary );

    ifs.seekg( 0, ifs.end );
    const size_t file_size = ifs.tellg();
    const size_t data_size = file_size - offset * 2;

    ifs.seekg( offset, ifs.beg );
    kvs::ValueArray<float> values( data_size / sizeof(float) );
    ifs.read( reinterpret_cast<char*>( values.data() ), values.byteSize() );
    if ( swap ) { kvs::Endian::Swap( values.data(), values.size() ); }

    return values;
}

inline kvs::ValueArray<float> Interleave(
    const kvs::ValueArray<float>& a,
    const kvs::ValueArray<float>& b,
    const kvs::ValueArray<float>& c )
{
    kvs::ValueArray<float> abc( a.size() * 3 );
    for ( size_t i = 0; i < a.size(); i++ )
    {
        abc[ 3 * i + 0 ] = a[i];
        abc[ 3 * i + 1 ] = b[i];
        abc[ 3 * i + 2 ] = c[i];
    }
    return abc;
}

inline bool IsAbsolutePath( const std::string& data_file )
{
    return kvs::File( data_file ).pathName( true ) == kvs::File( data_file ).pathName();
}

inline std::string AbsolutePath( const std::string& data_file )
{
    return kvs::Directory::Absolute( data_file );
}

}

namespace YYZVis
{

ZhongVolumeImporter::ZhongVolumeImporter( const std::string& filename )
{
    kvs::Json json( filename );
    this->exec( &json );
}

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
    const bool swap = ::ByteSwap( endian );
    std::vector<kvs::ValueArray<kvs::Real32>> temp;
    for ( auto& f : zhong_value )
    {
        const std::string data_file = f.get<std::string>();
        if ( ::IsAbsolutePath( data_file ) )
        {
            temp.push_back( ::ReadBinary( data_file, 4, swap ) );
        }
        else if ( data_file[0] == '~' )
        {
            const std::string filename = ::AbsolutePath( data_file );
            temp.push_back( ::ReadBinary( filename, 4, swap ) );
        }
        else
        {
            const std::string filename = path_name + kvs::File::Separator() + data_file;
            temp.push_back( ::ReadBinary( filename, 4, swap ) );
        }
    }

    auto values = ( temp.size() == 1 ) ? temp[0] : ::Interleave( temp[0], temp[1], temp[2] );

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
