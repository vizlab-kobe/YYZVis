#include "View.h"
#include <kvs/StochasticLineRenderer>
#include <kvs/ParticleBasedRenderer>
#include <kvs/Light>
#include <kvs/Math>

#if defined( JSST2019_TEST )
#include <kvs/PaintEventListener>
namespace
{

class FPSTimer : public kvs::PaintEventListener
{
private:
    local::View* m_view;

public:
    FPSTimer( local::View* view ):
        m_view( view )
    {
    }

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
            elapse_time += m_view->compositor().timer().sec();
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
    m_screen( app ),
    m_compositor( m_screen.scene() )
{
    this->setup();
    this->show();

#if defined( JSST2019_TEST )
    const kvs::Xform S = kvs::Xform::Scaling( kvs::Vec3::All( 1.3 ) );
    const kvs::Xform R = kvs::Xform::Rotation( kvs::Mat3::RotationX( 30 ) * kvs::Mat3::RotationY( -40 ) );
    const kvs::Xform x = S * R;
    m_screen.scene()->object("YinMesh")->multiplyXform( x );
    m_screen.scene()->object("YangMesh")->multiplyXform( x );
    m_screen.scene()->object("YinEdge")->multiplyXform( x );
    m_screen.scene()->object("YangEdge")->multiplyXform( x );
    m_screen.scene()->object("ZhongEdge")->multiplyXform( x );
    m_screen.scene()->object("YinVolume")->multiplyXform( x );
    m_screen.scene()->object("YangVolume")->multiplyXform( x );
    m_screen.scene()->object("ZhongVolume")->multiplyXform( x );
    m_screen.addEvent( new FPSTimer( this ) );
#endif
}

void View::setup()
{
    m_screen.setBackgroundColor( kvs::RGBColor::White() );
    m_screen.setSize( 800, 600 );

    this->setup_meshes();
    this->setup_edges();
    this->setup_particles();

    const size_t repeats = m_model->input().repeats;
    m_compositor.setRepetitionLevel( repeats );
//    m_compositor.enableLODControl();
    m_screen.setEvent( &m_compositor );
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
    typedef kvs::StochasticLineRenderer Renderer;

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
    typedef kvs::StochasticLineRenderer Renderer;

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

void View::setup_particles()
{
    const kvs::Indent indent( 4 );
    const bool enable_shading = m_model->input().enable_shading;
    typedef kvs::glsl::ParticleBasedRenderer Renderer;

    std::cout << "SETUP PARTICLES ..." << std::endl;

    // Yin
    {
        kvs::PointObject* object = m_model->newYinParticles();
        object->setName( "YinVolume" );
        object->print( std::cout << "YIN PARTICLE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YinVolume" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PointObject* object = m_model->newYangParticles();
        object->setName( "YangVolume" );
        object->print( std::cout << "YANG PARTICLE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "YangVolume" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PointObject* object = m_model->newZhongParticles();
        object->setName( "ZhongVolume" );
        object->print( std::cout << "ZHONG PARTICLE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setName( "ZhongVolume" );
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }
}

} // end of namespace local
