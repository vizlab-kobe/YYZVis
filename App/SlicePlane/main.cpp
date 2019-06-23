#include <string>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/glut/CheckBox>
#include <kvs/glut/CheckBoxGroup>
#include <kvs/Program>
#include <kvs/CommandLine>
#include <kvs/DivergingColorMap>
#include <kvs/TransferFunction>
#include <kvs/ParticleBasedRenderer>
#include <kvs/PolygonRenderer>
#include <kvs/LineRenderer>
#include <kvs/StochasticLineRenderer>
#include <kvs/StochasticPolygonRenderer>
#include <kvs/StochasticTetrahedraRenderer>
#include <kvs/StochasticRenderingCompositor>
#include <kvs/ExternalFaces>
#include <kvs/CellByCellMetropolisSampling>
#include <kvs/Isosurface>
#include <kvs/ColorMap>
#include <kvs/Vector3>
#include <kvs/Vector4>
#include <Lib/YinYangVolumeObject.h>
#include <Lib/YinYangGridSampling.h>
#include <Lib/ZhongVolumeObject.h>
#include <Lib/ZhongGridSampling.h>
#include <Lib/Edge.h>
#include <Lib/MarchingHexahedra.h>
#include <Lib/MarchingHexahedraTable.h>
#include <Lib/SlicePlane.h>

class CheckBox : public kvs::glut::CheckBox
{
private:
    std::string m_name; ///< object name
    kvs::Scene* m_scene; ///< scene

public:
    CheckBox( kvs::glut::Screen* screen, const std::string name ):
        kvs::glut::CheckBox( screen ),
        m_name( name ),
        m_scene( screen->scene() ) {}

    void stateChanged()
    {
        if ( !m_scene->hasObject( m_name ) ) { return; }
        if ( state() ) { m_scene->object( m_name )->show(); }
        else { m_scene->object( m_name )->hide(); }
    }
};

class Program : public kvs::Program
{
private:

    typedef YinYangVis::YinYangVolumeObject YinYangObject;
    typedef YinYangVis::ZhongVolumeObject ZhongObject;

    kvs::CommandLine m_command;

public:

    Program()
    {
        m_command.addOption( "dim_rad", "Dimension in the radial direction.", 1, true );
        m_command.addOption( "dim_lat", "Dimension in the latitude direction.", 1, true );
        m_command.addOption( "dim_lon", "Dimension in the longitude direction.", 1, true );
        m_command.addOption( "dim_zhong", "Dimension of the cubic zhong grid.", 1, true );
        m_command.addOption( "yin", "Filename of yin volume data.", 1, true );
        m_command.addOption( "yang", "Filename of yang volume data.", 1, true );
        m_command.addOption( "zhong", "Filename of zhong volume data.", 1, true );
        m_command.addOption( "repeat", "Number of repetitions for PBVR. (defulat: 1)", 1, false );
        m_command.addOption( "prev", "Use previous particle generation technique." );
        m_command.addOption( "shading", "Enable shading." );
        m_command.addHelpOption();
    }

    int exec( int argc, char** argv )
    {
        m_command.setArguments( argc, argv );
        if ( !m_command.parse() ) { return 1; }

        kvs::glut::Application app( argc, argv );
        kvs::glut::Screen screen( &app );
        screen.setSize( 800, 600 );
        screen.setBackgroundColor( kvs::RGBColor::White() );

        // Read volumes.
        std::cout << "READ VOLUMES ..." << std::endl;
        const std::string filename_yin = m_command.optionValue<std::string>("yin");
        YinYangObject* yin = this->import_yin( filename_yin );
        yin->setName("Yin");
        yin->print( std::cout << "YIN VOLUME DATA" << std::endl, kvs::Indent(2) );

        const std::string filename_yang = m_command.optionValue<std::string>("yang");
        YinYangObject* yang = this->import_yang( filename_yang );
        yang->setName("Yang");
        yang->print( std::cout << "YANG VOLUME DATA" << std::endl, kvs::Indent(2) );

        const std::string filename_zhong = m_command.optionValue<std::string>("zhong");
        ZhongObject* zhong = this->import_zhong( filename_zhong );
        zhong->setName("Zhong");
        zhong->print( std::cout << "ZHONG VOLUME DATA" << std::endl, kvs::Indent(2) );

        this->update_min_max( yin, yang, zhong );

        // Create objects.
		this->create_meshes( screen, yin, yang );
	        this->create_edges( screen, yin, yang, zhong );
		this->create_slice( screen, yin, yang, zhong );

        delete yin;
        delete yang;
        delete zhong;

        // Compositor.
        const size_t repeats = m_command.hasOption("repeat") ? m_command.optionValue<size_t>("repeat") : 1;

        // Check box.
        CheckBox check_yin_mesh( &screen, "YinMesh" );
        check_yin_mesh.setX( 10 );
        check_yin_mesh.setY( 10 );
        check_yin_mesh.setMargin( 10 );
        check_yin_mesh.setState( true );
        check_yin_mesh.setCaption( "Mesh (Yin)" );

        CheckBox check_yang_mesh( &screen, "YangMesh" );
        check_yang_mesh.setX( check_yin_mesh.x() );
        check_yang_mesh.setY( check_yin_mesh.y() + 20 );
        check_yang_mesh.setMargin( 10 );
        check_yang_mesh.setState( true );
        check_yang_mesh.setCaption( "Mesh (Yang)" );

        CheckBox check_yin_edge( &screen, "YinEdge" );
        check_yin_edge.setX( check_yang_mesh.x() );
        check_yin_edge.setY( check_yang_mesh.y() + 20 );
        check_yin_edge.setMargin( 10 );
        check_yin_edge.setState( true );
        check_yin_edge.setCaption( "Edge (Yin)" );

        CheckBox check_yang_edge( &screen, "YangEdge" );
        check_yang_edge.setX( check_yin_edge.x() );
        check_yang_edge.setY( check_yin_edge.y() + 20 );
        check_yang_edge.setMargin( 10 );
        check_yang_edge.setState( true );
        check_yang_edge.setCaption( "Edge (Yang)" );

        CheckBox check_zhong_edge( &screen, "ZhongEdge" );
        check_zhong_edge.setX( check_yang_edge.x() );
        check_zhong_edge.setY( check_yang_edge.y() + 20 );
        check_zhong_edge.setMargin( 10 );
        check_zhong_edge.setState( true );
        check_zhong_edge.setCaption( "Edge (Zhong)" );

	CheckBox check_yin_sliceplane( &screen, "YinSlice" );
        check_yin_sliceplane.setX( check_zhong_edge.x() );
        check_yin_sliceplane.setY( check_zhong_edge.y() + 20 );
        check_yin_sliceplane.setMargin( 10 );
        check_yin_sliceplane.setState( true );
        check_yin_sliceplane.setCaption( "Sliceplane (Yin)" );

	CheckBox check_yang_sliceplane( &screen, "YangSlice" );
        check_yang_sliceplane.setX( check_yin_sliceplane.x() );
        check_yang_sliceplane.setY( check_yin_sliceplane.y() + 20 );
        check_yang_sliceplane.setMargin( 10 );
        check_yang_sliceplane.setState( true );
        check_yang_sliceplane.setCaption( "Sliceplane (Yang)" );

       	CheckBox check_zhong_sliceplane( &screen, "ZhongSlice" );
        check_zhong_sliceplane.setX( check_yang_sliceplane.x() );
        check_zhong_sliceplane.setY( check_yang_sliceplane.y() + 20 );
        check_zhong_sliceplane.setMargin( 10 );
        check_zhong_sliceplane.setState( true );
        check_zhong_sliceplane.setCaption( "Sliceplane (zhong)" );

        kvs::glut::CheckBoxGroup group;
        group.add( &check_yin_mesh );
        group.add( &check_yang_mesh );
        group.add( &check_yin_edge );
        group.add( &check_yang_edge );
        group.add( &check_zhong_edge );
        group.add( &check_yin_sliceplane );
        group.add( &check_yang_sliceplane );
	group.add( &check_zhong_sliceplane );
	
        group.show();
        screen.show();

        kvs::Light::SetModelTwoSide( true );

        return app.run();
    }

private:

    YinYangObject* import_yin( const std::string& filename )
    {
        const size_t dim_rad = m_command.optionValue<size_t>("dim_rad");
        const size_t dim_lat = m_command.optionValue<size_t>("dim_lat");
        const size_t dim_lon = m_command.optionValue<size_t>("dim_lon");

        YinYangObject* volume = new YinYangObject();
        volume->setGridTypeToYin();
        volume->setDimR( dim_rad );
        volume->setDimTheta( dim_lat );
        volume->setDimPhi( dim_lon );
        volume->setVeclen( 1 );
        volume->calculateCoords();
        volume->readValues( filename );
        volume->updateMinMaxCoords();
        volume->updateMinMaxValues();
        return volume;
    }

    YinYangObject* import_yang( const std::string& filename )
    {
        const size_t dim_rad = m_command.optionValue<size_t>("dim_rad");
        const size_t dim_lat = m_command.optionValue<size_t>("dim_lat");
        const size_t dim_lon = m_command.optionValue<size_t>("dim_lon");

        YinYangObject* volume = new YinYangObject();
        volume->setGridTypeToYang();
        volume->setDimR( dim_rad );
        volume->setDimTheta( dim_lat );
        volume->setDimPhi( dim_lon );
        volume->setVeclen( 1 );
        volume->calculateCoords();
        volume->readValues( filename );
        volume->updateMinMaxCoords();
        volume->updateMinMaxValues();
        return volume;
    }

    ZhongObject* import_zhong( const std::string& filename )
    {
        const size_t dim_rad = m_command.optionValue<size_t>("dim_rad");
        const size_t dim_zhong = m_command.optionValue<size_t>("dim_zhong");

        ZhongObject* volume = new ZhongObject();
        volume->setDimR( dim_rad );
        volume->setDim( dim_zhong );
        volume->setVeclen( 1 );
        volume->calculateCoords();
        volume->readValues( filename );
        volume->updateMinMaxCoords();
        volume->updateMinMaxValues();
        return volume;
    }

    void update_min_max( YinYangObject* yin, YinYangObject* yang, ZhongObject* zhong )
    {
        const kvs::Real32 min_value = kvs::Math::Min( yin->minValue(), yang->minValue(), zhong->minValue() );
        const kvs::Real32 max_value = kvs::Math::Max( yin->maxValue(), yang->maxValue(), zhong->maxValue() );
        yin->setMinMaxValues( min_value, max_value );
        yang->setMinMaxValues( min_value, max_value );
        zhong->setMinMaxValues( min_value, max_value );

        const kvs::Vec3& min_coord0 = yin->minObjectCoord();
        const kvs::Vec3& min_coord1 = yang->minObjectCoord();
        const kvs::Vec3& min_coord2 = zhong->minObjectCoord();
        const kvs::Vec3& max_coord0 = yin->maxObjectCoord();
        const kvs::Vec3& max_coord1 = yang->maxObjectCoord();
        const kvs::Vec3& max_coord2 = zhong->maxObjectCoord();
        const kvs::Real32 min_x = kvs::Math::Min( min_coord0.x(), min_coord1.x(), min_coord2.x() );
        const kvs::Real32 min_y = kvs::Math::Min( min_coord0.y(), min_coord1.y(), min_coord2.y() );
        const kvs::Real32 min_z = kvs::Math::Min( min_coord0.z(), min_coord1.z(), min_coord2.z() );
        const kvs::Real32 max_x = kvs::Math::Min( max_coord0.x(), max_coord1.x(), max_coord2.x() );
        const kvs::Real32 max_y = kvs::Math::Min( max_coord0.y(), max_coord1.y(), max_coord2.y() );
        const kvs::Real32 max_z = kvs::Math::Min( max_coord0.z(), max_coord1.z(), max_coord2.z() );
        const kvs::Vec3 min_coord( min_x, min_y, min_z );
        const kvs::Vec3 max_coord( max_x, max_y, max_z );
        yin->setMinMaxObjectCoords( min_coord, max_coord );
        yin->setMinMaxExternalCoords( min_coord, max_coord );
        yang->setMinMaxObjectCoords( min_coord, max_coord );
        yang->setMinMaxExternalCoords( min_coord, max_coord );
        zhong->setMinMaxObjectCoords( min_coord, max_coord );
        zhong->setMinMaxExternalCoords( min_coord, max_coord );
    }

    void create_meshes(
        kvs::glut::Screen& screen,
        const YinYangObject* yin,
        const YinYangObject* yang )
    {
        std::cout << "CREATE MESHES ..." << std::endl;

        kvs::LineObject* yin_mesh = YinYangVis::Edge::CreateLineMeshObject( yin );
        yin_mesh->setName("YinMesh");
        yin_mesh->setColor( kvs::RGBColor::Black() );
        yin_mesh->print( std::cout << "YIN MESH DATA" << std::endl, kvs::Indent(2) );
        this->register_object( screen, yin_mesh, new kvs::LineRenderer() );

        kvs::LineObject* yang_mesh = YinYangVis::Edge::CreateLineMeshObject( yang );
        yang_mesh->setName("YangMesh");
        yang_mesh->setColor( kvs::RGBColor::Black() );
        yang_mesh->print( std::cout << "YANG MESH DATA" << std::endl, kvs::Indent(2) );
        this->register_object( screen, yang_mesh, new kvs::LineRenderer() );
    }

    void create_edges(
        kvs::glut::Screen& screen,
        const YinYangObject* yin,
        const YinYangObject* yang,
        const ZhongObject* zhong )
    {
        std::cout << "CREATE EDGES ..." << std::endl;

        kvs::LineObject* yin_edge = YinYangVis::Edge::CreateLineEdgeObject( yin );
        yin_edge->setName("YinEdge");
        yin_edge->setColor( kvs::RGBColor::Black() );
        yin_edge->setSize( 2 );
        yin_edge->print( std::cout << "YIN EDGE DATA" << std::endl, kvs::Indent(2) );
        this->register_object( screen, yin_edge, new kvs::LineRenderer() );

        kvs::LineObject* yang_edge = YinYangVis::Edge::CreateLineEdgeObject( yang );
        yang_edge->setName("YangEdge");
        yang_edge->setColor( kvs::RGBColor::Black() );
        yang_edge->setSize( 2 );
        yang_edge->print( std::cout << "YANG EDGE DATA" << std::endl, kvs::Indent(2) );
        this->register_object( screen, yang_edge, new kvs::LineRenderer() );

        kvs::LineObject* zhong_edge = YinYangVis::Edge::CreateLineEdgeObject( zhong );
        zhong_edge->setName("ZhongEdge");
        zhong_edge->setColor( kvs::RGBColor::Black() );
        zhong_edge->setSize( 2 );
        zhong_edge->print( std::cout << "ZHONG EDGE DATA" << std::endl, kvs::Indent(2) );
        this->register_object( screen, zhong_edge, new kvs::LineRenderer() );
    }

  void create_slice(
			 kvs::glut::Screen& screen,
			 const YinYangObject* yin,
			 const YinYangObject* yang,
			 const ZhongObject* zhong)
  {
    std::cout << "CREATE sliceplane ..." << std::endl;
    kvs::Vec3 point(0,0,0);
    kvs::Vec3 normal(-1,0,-3);

      kvs::UnstructuredVolumeObject* yin_volume = YinYangObject::ToUnstructuredVolumeObject( yin );                  
      kvs::PolygonObject* yin_slice = new kvs::SlicePlane( yin_volume, point, normal);
      yin_slice->setName("YinSlice");
      yin_slice->setPolygonType(kvs::PolygonObject::Triangle);
      yin_slice->print( std::cout << "YIN SLICE DATA" << std::endl, kvs::Indent(2) );
      screen.registerObject( yin_slice,new kvs::PolygonRenderer() );
                                     
      delete yin_volume;

      kvs::UnstructuredVolumeObject* yang_volume = YinYangObject::ToUnstructuredVolumeObject( yang );             
      kvs::PolygonObject* yang_slice = new kvs::SlicePlane( yang_volume, point, normal);
      yang_slice->setName("YangSlice");
      yang_slice->setPolygonType(kvs::PolygonObject::Triangle);
      yang_slice->print( std::cout << "YANG SLICE DATA" << std::endl, kvs::Indent(2) );
      screen.registerObject( yang_slice,new kvs::PolygonRenderer() );                

      delete yang_volume;

      kvs::UnstructuredVolumeObject* zhong_volume = ZhongObject::ToUnstructuredVolumeObject( zhong );     
      kvs::PolygonObject* zhong_slice = new kvs::SlicePlane( zhong_volume, point, normal);
      zhong_slice->setName("ZhongSlice");
      zhong_slice->setPolygonType(kvs::PolygonObject::Triangle);
      zhong_slice->print( std::cout << "ZHONG SLICE DATA" << std::endl, kvs::Indent(2) );
      screen.registerObject( zhong_slice,new kvs::PolygonRenderer() );                              

      delete zhong_volume;
  }
  
    void register_object( kvs::glut::Screen& screen, kvs::ObjectBase* object, kvs::RendererBase* renderer )
    {
        kvs::Xform x = kvs::Xform::Rotation(
            kvs::Mat3::RotationY( 30 ) );
        object->multiplyXform( x );
        screen.registerObject( object, renderer );
    }
};

int main( int argc, char** argv )
{
    Program p;
    return p.start( argc, argv );
}

