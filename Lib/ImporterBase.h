#pragma once
#include <kvs/ImporterBase>
#include <kvs/Module>
#include <kvs/ValueArray>
#include <kvs/Endian>
#include <kvs/File>
#include <kvs/Directory>
#include <iostream>
#include <fstream>
#include <vector>


namespace YYZVis
{

class ImporterBase : public kvs::ImporterBase
{
    kvsModule( YYZVis::ImporterBase, Importer );
    kvsModuleBaseClass( kvs::ImporterBase );

public:
    ImporterBase() {}

protected:
    bool needsByteSwap( const std::string& endian );
    bool isAbsolutePath( const std::string& data_file );
    std::string absolutePath( const std::string& data_file );

    template <typename T>
    kvs::ValueArray<T> readBinary( const std::string& filename, const size_t offset, const bool swap );

    template <typename T>
    kvs::ValueArray<T> interleaveArrays( const std::vector<kvs::ValueArray<T>>& arrays );
};

inline bool ImporterBase::needsByteSwap( const std::string& endian )
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

inline bool ImporterBase::isAbsolutePath( const std::string& data_file )
{
    return kvs::File( data_file ).pathName( true ) == kvs::File( data_file ).pathName();
}

inline std::string ImporterBase::absolutePath( const std::string& data_file )
{
    return kvs::Directory::Absolute( data_file );
}

template <typename T>
inline kvs::ValueArray<T> ImporterBase::readBinary(
    const std::string& filename,
    const size_t offset,
    const bool swap )
{
    std::ifstream ifs( filename.c_str(), std::ifstream::binary );
    if ( !ifs )
    {
        kvsMessageError() << "Cannot open " << filename << std::endl;
        return kvs::ValueArray<T>();
    }

    ifs.seekg( 0, ifs.end );
    const size_t file_size = ifs.tellg();
    const size_t data_size = file_size - offset * 2;

    ifs.seekg( offset, ifs.beg );
    kvs::ValueArray<float> values( data_size / sizeof(float) );
    ifs.read( reinterpret_cast<char*>( values.data() ), values.byteSize() );
    if ( swap ) { kvs::Endian::Swap( values.data(), values.size() ); }

    return values;
}

template <typename T>
inline kvs::ValueArray<T> ImporterBase::interleaveArrays( const std::vector<kvs::ValueArray<T>>& arrays )
{
    const size_t narrays = arrays.size();
    if ( narrays == 1 ) { return arrays[0]; }

    const size_t array_size = arrays[0].size();
    kvs::ValueArray<T> result( array_size * narrays );
    for ( size_t i = 0; i < array_size; ++i )
    {
        const size_t offset = narrays * i;
        for ( size_t j = 0; j < narrays; ++j )
        {
            result[ offset + j ] = arrays[j][i];
        }
    }

    return result;
}

} // end of namespace YYZVis
