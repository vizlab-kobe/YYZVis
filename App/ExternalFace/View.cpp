#include "View.h"
#include <kvs/LineRenderer>
#include <kvs/PolygonRenderer>
#include <kvs/Light>


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
    this->setup_faces();
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

void View::setup_faces()
{
    const kvs::Indent indent( 4 );
    const bool enable_shading = true;
    typedef kvs::glsl::PolygonRenderer Renderer;

    std::cout << "SETUP FACES ..." << std::endl;

    // Yin
    {
        kvs::PolygonObject* object = m_model->newYinFaces();
        object->setName( "YinFace" );
        object->print( std::cout << "YIN FACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Yang
    {
        kvs::PolygonObject* object = m_model->newYangFaces();
        object->setName( "YangFace" );
        object->print( std::cout << "YANG FACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }

    // Zhong
    {
        kvs::PolygonObject* object = m_model->newZhongFaces();
        object->setName( "ZhongFace" );
        object->print( std::cout << "ZHONG FACE DATA" << std::endl, indent );

        Renderer* renderer = new Renderer();
        renderer->setEnabledShading( enable_shading );
        m_screen.registerObject( object, renderer );
    }
}

} // end of namespace local
