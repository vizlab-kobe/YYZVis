/****************************************************************************/
/**
 *  @file   HexahedralCell.h
 *  @author Naohisa Sakamoto
 */
/*----------------------------------------------------------------------------
 *
 *  Copyright (c) Visualization Laboratory, Kyoto University.
 *  All rights reserved.
 *  See http://www.viz.media.kyoto-u.ac.jp/kvs/copyright/ for details.
 *
 *  $Id: HexahedralCell.h 1792 2014-07-31 04:50:42Z naohisa.sakamoto@gmail.com $
 */
/****************************************************************************/
#pragma once

#include <kvs/Type>
#include <kvs/Vector4>
#include <kvs/Matrix44>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/CellBase>


namespace local
{

/*===========================================================================*/
/**
 *  @brief  HexahedralCell class.
 */
/*===========================================================================*/
class HexahedralCell : public kvs::CellBase
{
public:

    typedef kvs::CellBase BaseClass;

public:

    HexahedralCell( const kvs::UnstructuredVolumeObject* volume );

    void updateInterpolationFunctions( const kvs::Vec3& local ) const;
    void updateDifferentialFunctions( const kvs::Vec3& local ) const;
    const kvs::Real32 volume() const;
};

} // end of namespace local
