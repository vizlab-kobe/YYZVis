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
#include <kvs/StochasticLineRenderer>
#include <kvs/StochasticPolygonRenderer>
#include <kvs/StochasticTetrahedraRenderer>
#include <kvs/StochasticRenderingCompositor>

#include <Lib/YinYangVolumeObject.h>
#include <Lib/YinYangGridSampling.h>
#include <Lib/ZhongVolumeObject.h>
#include <Lib/ZhongGridSampling.h>
#include <Lib/Edge.h>


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
        m_command.addHelpOption();
    }

    int exec( int argc, char** argv )
    {
        m_command.setArguments( argc, argv );
        if ( !m_command.parse() ) { return 1; }

        kvs::glut::Application app( argc, argv );
        kvs::glut::Screen screen( &app );
        screen.setBackgroundColor( kvs::RGBColor::White() );

        // Read volumes.
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
        this->create_particles( screen, yin, yang, zhong );

        std::cout << screen.scene()->numberOfObjects() << std::endl;
        for ( size_t i = 0; i < screen.scene()->numberOfObjects() - 1; i++ )
        {
            std::cout << screen.scene()->object(i+1)->name() << std::endl;
        }

        delete yin;
        delete yang;
        delete zhong;

        // Compositor.
        const size_t repeats = m_command.hasOption("repeat") ? m_command.optionValue<size_t>("repeat") : 1;
        kvs::StochasticRenderingCompositor* compositor = new kvs::StochasticRenderingCompositor( screen.scene() );
        compositor->setRepetitionLevel( repeats );
        compositor->disableLODControl();
        screen.setEvent( compositor );

        // Check box.
        CheckBox check_yin_volume( &screen, "Yin" );
        check_yin_volume.setX( 10 );
        check_yin_volume.setY( 10 );
        check_yin_volume.setMargin( 10 );
        check_yin_volume.setState( true );
        check_yin_volume.setCaption( "Yin" );

        CheckBox check_yang_volume( &screen, "Yang" );
        check_yang_volume.setX( check_yin_volume.x() );
        check_yang_volume.setY( check_yin_volume.y() + 20 );
        check_yang_volume.setMargin( 10 );
        check_yang_volume.setState( true );
        check_yang_volume.setCaption( "Yang" );

        CheckBox check_zhong_volume( &screen, "Zhong" );
        check_zhong_volume.setX( check_yang_volume.x() );
        check_zhong_volume.setY( check_yang_volume.y() + 20 );
        check_zhong_volume.setMargin( 10 );
        check_zhong_volume.setState( true );
        check_zhong_volume.setCaption( "Zhong" );

        CheckBox check_yin_mesh( &screen, "YinMesh" );
        check_yin_mesh.setX( check_zhong_volume.x() );
        check_yin_mesh.setY( check_zhong_volume.y() + 20 );
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

        kvs::glut::CheckBoxGroup group;
        group.add( &check_yin_volume );
        group.add( &check_yang_volume );
        group.add( &check_zhong_volume );
        group.add( &check_yin_mesh );
        group.add( &check_yang_mesh );
        group.add( &check_yin_edge );
        group.add( &check_yang_edge );
        group.add( &check_zhong_edge );

        group.show();
        screen.show();

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

    void create_meshes( kvs::glut::Screen& screen, const YinYangObject* yin, const YinYangObject* yang )
    {
        kvs::LineObject* yin_mesh = YinYangVis::Edge::CreateLineMeshObject( yin );
        yin_mesh->setName("YinMesh");
        screen.registerObject( yin_mesh, new kvs::StochasticLineRenderer() );

        kvs::LineObject* yang_mesh = YinYangVis::Edge::CreateLineMeshObject( yang );
        yang_mesh->setName("YangMesh");
        screen.registerObject( yang_mesh, new kvs::StochasticLineRenderer() );
    }

    void create_edges( kvs::glut::Screen& screen, const YinYangObject* yin, const YinYangObject* yang, const ZhongObject* zhong )
    {
        kvs::LineObject* yin_edge = YinYangVis::Edge::CreateLineEdgeObject( yin );
        yin_edge->setName("YinEdge");
        screen.registerObject( yin_edge, new kvs::StochasticLineRenderer() );

        kvs::LineObject* yang_edge = YinYangVis::Edge::CreateLineEdgeObject( yang );
        yang_edge->setName("YangEdge");
        screen.registerObject( yang_edge, new kvs::StochasticLineRenderer() );

        kvs::LineObject* zhong_edge = YinYangVis::Edge::CreateLineEdgeObject( zhong );
        zhong_edge->setName("ZhongEdge");
        screen.registerObject( zhong_edge, new kvs::StochasticLineRenderer() );
    }

    void create_particles( kvs::glut::Screen& screen, const YinYangObject* yin, const YinYangObject* yang, const ZhongObject* zhong )
    {
        // Transfer function.
        kvs::ColorMap cmap = kvs::DivergingColorMap::CoolWarm( 256 );
        kvs::OpacityMap omap( 256 );
        omap.addPoint( 0, 1.0 );
        omap.addPoint( 90, 0.0 );
        omap.addPoint( 180, 0.0 );
        omap.addPoint( 255, 1.0 );
        omap.create();
        const kvs::TransferFunction tfunc( cmap, omap );

        // Parameters for PBVR.
        const size_t repeats = m_command.hasOption("repeat") ? m_command.optionValue<size_t>("repeat") : 1;
        const size_t subpixels = 1; // fixed to '1'
        const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repeats ) ) );
        const float step = 0.1f;

        // Generate particles.
        std::cout << "PARTICLE GENERATION..." << std::endl;

        // Yin
        kvs::Timer timer( kvs::Timer::Start );
        kvs::PointObject* yin_particle = new YinYangVis::YinYangGridSampling( yin, level, step, tfunc );
        timer.stop();
        std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;
        yin_particle->setName( yin->name() );
        yin_particle->print( std::cout << "YIN PARTICLE DATA" << std::endl, kvs::Indent(2) );

        kvs::glsl::ParticleBasedRenderer* yin_particle_renderer = new kvs::glsl::ParticleBasedRenderer();
        yin_particle_renderer->disableShading();
        screen.registerObject( yin_particle, yin_particle_renderer );

        // Ynag
        timer.start();
        kvs::PointObject* yang_particle = new YinYangVis::YinYangGridSampling( yang, level, step, tfunc );
        timer.stop();
        std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;
        yang_particle->setName( yang->name() );
        yang_particle->print( std::cout << "YANG PARTICLE DATA" << std::endl, kvs::Indent(2) );

        kvs::glsl::ParticleBasedRenderer* yang_particle_renderer = new kvs::glsl::ParticleBasedRenderer();
        yang_particle_renderer->disableShading();
        screen.registerObject( yang_particle, yang_particle_renderer );

        // Zhong
        timer.start();
        kvs::PointObject* zhong_particle = new YinYangVis::ZhongGridSampling( zhong, level, step, tfunc );
        timer.stop();
        std::cout << std::endl << "Particle generation time: " << timer.sec() << " [sec]" << std::endl;
        zhong_particle->setName( zhong->name() );
        zhong_particle->print( std::cout << "ZHONG PARTICLE DATA" << std::endl, kvs::Indent(2) );

        kvs::glsl::ParticleBasedRenderer* zhong_particle_renderer = new kvs::glsl::ParticleBasedRenderer();
        zhong_particle_renderer->disableShading();
        screen.registerObject( zhong_particle, zhong_particle_renderer );
    }
};

int main( int argc, char** argv )
{
    Program p;
    return p.start( argc, argv );
}
