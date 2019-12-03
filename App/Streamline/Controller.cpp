#include "Controller.h"
#include <kvs/Time>
#include <kvs/Date>


namespace local
{

Controller::KeyPressEvent::KeyPressEvent( local::Controller* controller ):
    m_controller( controller )
{
    m_controller->view()->screen().addEvent( this );
}

void Controller::KeyPressEvent::update( kvs::KeyEvent* event )
{
    switch ( event->key() )
    {
    case kvs::Key::i:
    {
        if ( m_controller->isShown() ) { m_controller->hide(); }
        else { m_controller->show(); }
        break;
    }
    case kvs::Key::s:
    {
        kvs::Date date;
        kvs::Time time;
        std::string filename = date.today().toString("") + "_" + time.now().toString("") + ".bmp";
        kvs::ColorImage image = m_controller->view()->screen().capture();
        image.write( filename );
        std::cout << "Capture: " << filename << std::endl;
    }
    default: break;
    }
}

Controller::Controller( local::Model* model, local::View* view ):
    m_view( view ),
//    m_isosurface_checkbox( view, "Isosurace", "YinIso", "YangIso", "ZhongIso" ),
    m_mesh_checkbox( view, "Mesh", "YinMesh", "YangMesh" ),
    m_edge_checkbox( view, "Edge", "YinEdge", "YangEdge", "ZhongEdge" ),
//    m_isovalue_slider( model, view ),
    m_key_event( this )
{
//    m_isosurface_checkbox.setMargin( 10 );
//    m_isosurface_checkbox.setPosition( 0, 0 );
//    m_isosurface_checkbox.show();

    m_mesh_checkbox.setMargin( 10 );
    m_mesh_checkbox.setPosition( 0, 150 );
    m_mesh_checkbox.show();

    m_edge_checkbox.setMargin( 10 );
    m_edge_checkbox.setPosition( 0, 270 );
    m_edge_checkbox.show();

//    m_isovalue_slider.setMargin( 10 );
//    m_isovalue_slider.setPosition( m_view->screen().width() - m_isovalue_slider.width() - 10, 0 );
//    m_isovalue_slider.show();
}

void Controller::show()
{
//    m_isosurface_checkbox.show();
    m_mesh_checkbox.show();
    m_edge_checkbox.show();
//    m_isovalue_slider.show();
}

void Controller::hide()
{
//    m_isosurface_checkbox.hide();
    m_mesh_checkbox.hide();
    m_edge_checkbox.hide();
//    m_isovalue_slider.hide();
}

bool Controller::isShown()
{
    return
//        m_isosurface_checkbox.isShown() ||
        m_mesh_checkbox.isShown() ||
        m_edge_checkbox.isShown();
//        m_isovalue_slider.isShown();
}

} // end of namespace local
