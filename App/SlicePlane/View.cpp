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
            elapse_time += scene()->renderer("YinSlice")->timer().sec();
            elapse_time += scene()->renderer("YangSlice")->timer().sec();
            elapse_time += scene()->renderer("ZhongSlice")->timer().sec();
            elapse_time += scene()->renderer("YinSlice2")->timer().sec();
            elapse_time += scene()->renderer("YangSlice2")->timer().sec();
            elapse_time += scene()->renderer("ZhongSlice2")->timer().sec();
            elapse_time += scene()->renderer("YinSlice3")->timer().sec();
            elapse_time += scene()->renderer("YangSlice3")->timer().sec();
            elapse_time += scene()->renderer("ZhongSlice3")->timer().sec();
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
    this->setup_slices();

#if defined( JSST2019_TEST )
    const kvs::Xform S = kvs::Xform::Scaling( kvs::Vec3::All( 1.3 ) );
    const kvs::Xform R = kvs::Xform::Rotation( kvs::Mat3::RotationX( 30 ) * kvs::Mat3::RotationY( -40 ) );
    const kvs::Xform x = S * R;
    m_screen.scene()->object("YinMesh")->multiplyXform( x );
    m_screen.scene()->object("YangMesh")->multiplyXform( x );
    m_screen.scene()->object("YinEdge")->multiplyXform( x );
    m_screen.scene()->object("YangEdge")->multiplyXform( x );
    m_screen.scene()->object("ZhongEdge")->multiplyXform( x );
    m_screen.scene()->object("YinSlice")->multiplyXform( x );
    m_screen.scene()->object("YangSlice")->multiplyXform( x );
    m_screen.scene()->object("ZhongSlice")->multiplyXform( x );
    m_screen.scene()->object("YinSlice2")->multiplyXform( x );
    m_screen.scene()->object("YangSlice2")->multiplyXform( x );
    m_screen.scene()->object("ZhongSlice2")->multiplyXform( x );
    m_screen.scene()->object("YinSlice3")->multiplyXform( x );
    m_screen.scene()->object("YangSlice3")->multiplyXform( x );
    m_screen.scene()->object("ZhongSlice3")->multiplyXform( x );
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

void View::setup_slices()
{
    const kvs::Indent indent( 4 );
    const bool enable_shading = true;
    typedef kvs::glsl::PolygonRenderer Renderer;

    std::cout << "SETUP SLICES ..." << std::endl;

    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinSlice();
        object->setName( "YinSlice" );
        object->print( std::cout << "YIN SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinSlice" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangSlice();
        object->setName( "YangSlice" );
        object->print( std::cout << "YANG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangSlice" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongSlice();
        object->setName( "ZhongSlice" );
        object->print( std::cout << "ZHONG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongSlice" );
        renderer->setEnabledShading( enable_shading );
        renderer->setPolygonOffset( -0.001 );
        m_screen.registerObject( object, renderer );
    }

#if JSST2019_TEST
    const kvs::Vec3 plane_point = m_model->planePoint();
    const kvs::Vec3 plane_normal = m_model->planeNormal();
    const kvs::Vec3 normal2( 1, 0, 0 );
    const kvs::Vec3 normal3( 0, 1, 0 );

    m_model->setPlane( plane_point, normal2 );
    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinSlice();
        object->setName( "YinSlice2" );
        object->print( std::cout << "YIN SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinSlice2" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangSlice();
        object->setName( "YangSlice2" );
        object->print( std::cout << "YANG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangSlice2" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongSlice();
        object->setName( "ZhongSlice2" );
        object->print( std::cout << "ZHONG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongSlice2" );
        renderer->setEnabledShading( enable_shading );
        renderer->setPolygonOffset( -0.001 );
        m_screen.registerObject( object, renderer );
    }

    m_model->setPlane( plane_point, normal3 );
    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinSlice();
        object->setName( "YinSlice3" );
        object->print( std::cout << "YIN SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinSlice3" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangSlice();
        object->setName( "YangSlice3" );
        object->print( std::cout << "YANG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangSlice3" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongSlice();
        object->setName( "ZhongSlice3" );
        object->print( std::cout << "ZHONG SLICE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongSlice3" );
        renderer->setEnabledShading( enable_shading );
        renderer->setPolygonOffset( -0.001 );
        m_screen.registerObject( object, renderer );
    }
#endif
}

} // end of namespace local
