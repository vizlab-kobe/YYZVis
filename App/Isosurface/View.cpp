#include "View.h"
#include <kvs/LineRenderer>
#include <kvs/PolygonRenderer>
#include <kvs/Light>
#include <kvs/Math>

#if defined( JSST2019_TEST )
#include <kvs/PaintEventListener>
namespace
{

class FPSTimer : public kvs::PaintEventListener
{
    void update()
    {
        static int counter = 0;
        static float elapse_time = 0.0f;
        const int N = 50;
        if ( counter == N )
        {
            std::cout << "FPS: " << counter / elapse_time << std::endl;
            elapse_time = 0.0f;
            counter = 0;
        }
        else
        {
            elapse_time += scene()->renderer("YinIso")->timer().sec();
            elapse_time += scene()->renderer("YangIso")->timer().sec();
            elapse_time += scene()->renderer("ZhongIso")->timer().sec();
            elapse_time += scene()->renderer("YinIso2")->timer().sec();
            elapse_time += scene()->renderer("YangIso2")->timer().sec();
            elapse_time += scene()->renderer("ZhongIso2")->timer().sec();
            elapse_time += scene()->renderer("YinIso3")->timer().sec();
            elapse_time += scene()->renderer("YangIso3")->timer().sec();
            elapse_time += scene()->renderer("ZhongIso3")->timer().sec();
            counter++;
        }
    }
};
}
#endif

namespace local
{

View::View( kvs::glut::Application* app, local::Model* model ):
    m_model( model ),
    m_screen( app )
{
    this->setup();
    this->show();
}

void View::setup()
{
    m_screen.setBackgroundColor( kvs::RGBColor::White() );
    m_screen.setSize( 800, 600 );

    this->setup_meshes();
    this->setup_edges();
    this->setup_isosurfaces();

#if defined( JSST2019_TEST )
    const kvs::Xform S = kvs::Xform::Scaling( kvs::Vec3::Constant( 1.3 ) );
    const kvs::Xform R = kvs::Xform::Rotation( kvs::Mat3::RotationX( 30 ) * kvs::Mat3::RotationY( -40 ) );
    const kvs::Xform x = S * R;
    m_screen.scene()->object("YinMesh")->multiplyXform( x );
    m_screen.scene()->object("YangMesh")->multiplyXform( x );
    m_screen.scene()->object("YinEdge")->multiplyXform( x );
    m_screen.scene()->object("YangEdge")->multiplyXform( x );
    m_screen.scene()->object("ZhongEdge")->multiplyXform( x );
    m_screen.scene()->object("YinIso")->multiplyXform( x );
    m_screen.scene()->object("YangIso")->multiplyXform( x );
    m_screen.scene()->object("ZhongIso")->multiplyXform( x );
    m_screen.scene()->object("YinIso2")->multiplyXform( x );
    m_screen.scene()->object("YangIso2")->multiplyXform( x );
    m_screen.scene()->object("ZhongIso2")->multiplyXform( x );
    m_screen.scene()->object("YinIso3")->multiplyXform( x );
    m_screen.scene()->object("YangIso3")->multiplyXform( x );
    m_screen.scene()->object("ZhongIso3")->multiplyXform( x );
    m_screen.addEvent( new FPSTimer );
#endif
}

void View::show()
{
    m_screen.show();
    kvs::Light::SetModelTwoSide( true );
}

void View::redraw()
{
    m_screen.redraw();
}

void View::setup_meshes()
{
    const kvs::Indent indent( 4 );
    const kvs::RGBColor color = kvs::RGBColor::Black();
    typedef kvs::LineRenderer Renderer;

    std::cout << "SETUP MESHES ..." << std::endl;

    // Yin
    {
        kvs::LineObject* mesh = m_model->newYinMeshes();
        mesh->setName( "YinMesh" );
        mesh->setColor( color );
        mesh->print( std::cout << "YIN MESH DATA" << std::endl, indent );
        m_screen.registerObject( mesh, new Renderer() );
    }

    // Yang
    {
        kvs::LineObject* mesh = m_model->newYangMeshes();
        mesh->setName( "YangMesh" );
        mesh->setColor( color );
        mesh->print( std::cout << "YANG MESH DATA" << std::endl, indent );
        m_screen.registerObject( mesh, new Renderer() );
    }
}

void View::setup_edges()
{
    const kvs::Indent indent( 4 );
    const kvs::RGBColor color = kvs::RGBColor::Black();
    const float size = 2.0f;
    typedef kvs::LineRenderer Renderer;

    std::cout << "SETUP EDGES ..." << std::endl;

    // Yin
    {
        kvs::LineObject* edge = m_model->newYinEdges();
        edge->setName( "YinEdge" );
        edge->setColor( color );
        edge->setSize( size );
        edge->print( std::cout << "YIN EDGE DATA" << std::endl, indent );
        m_screen.registerObject( edge, new Renderer() );
    }

    // Yang
    {
        kvs::LineObject* edge = m_model->newYangEdges();
        edge->setName( "YangEdge" );
        edge->setColor( color );
        edge->setSize( size );
        edge->print( std::cout << "YANG EDGE DATA" << std::endl, indent );
        m_screen.registerObject( edge, new Renderer() );
    }

    // Zhong
    {
        kvs::LineObject* edge = m_model->newZhongEdges();
        edge->setName( "ZhongEdge" );
        edge->setColor( color );
        edge->setSize( size );
        edge->print( std::cout << "ZHONG EDGE DATA" << std::endl, indent );
        m_screen.registerObject( edge, new Renderer() );
    }
}

void View::setup_isosurfaces()
{
    const kvs::Indent indent( 4 );
    const bool enable_shading = true;
    typedef kvs::glsl::PolygonRenderer Renderer;

    std::cout << "SETUP ISOSURFACES ..." << std::endl;

    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinIsosurfaces();
        object->setName( "YinIso" );
        object->print( std::cout << "YIN ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinIso" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangIsosurfaces();
        object->setName( "YangIso" );
        object->print( std::cout << "YANG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangIso" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongIsosurfaces();
        object->setName( "ZhongIso" );
        object->print( std::cout << "ZHONG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongIso" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

#if JSST2019_TEST
    const float min_value = m_model->constYinVolume().minValue();
    const float max_value = m_model->constYinVolume().maxValue();
    const float ratio2 = 0.6f;
    const float ratio3 = 0.4f;

    m_model->setIsovalue( kvs::Math::Mix( min_value, max_value, ratio2 ) );
    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinIsosurfaces();
        object->setName( "YinIso2" );
        object->print( std::cout << "YIN ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinIso2" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangIsosurfaces();
        object->setName( "YangIso2" );
        object->print( std::cout << "YANG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangIso2" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongIsosurfaces();
        object->setName( "ZhongIso2" );
        object->print( std::cout << "ZHONG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongIso2" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    m_model->setIsovalue( kvs::Math::Mix( min_value, max_value, ratio3 ) );
    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinIsosurfaces();
        object->setName( "YinIso3" );
        object->print( std::cout << "YIN ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinIso3" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangIsosurfaces();
        object->setName( "YangIso3" );
        object->print( std::cout << "YANG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangIso3" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongIsosurfaces();
        object->setName( "ZhongIso3" );
        object->print( std::cout << "ZHONG ISOSURFACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongIso3" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }
#endif
}

} // end of namespace local
