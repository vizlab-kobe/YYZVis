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
    m_volume_checkbox( view, "Volume", "YinVolume", "YangVolume", "ZhongVolume" ),
    m_mesh_checkbox( view, "Mesh", "YinMesh", "YangMesh" ),
    m_edge_checkbox( view, "Edge", "YinEdge", "YangEdge", "ZhongEdge" ),
    m_key_event( this )
{
    m_volume_checkbox.setMargin( 10 );
    m_volume_checkbox.setPosition( 0, 0 );
    m_volume_checkbox.show();

    m_mesh_checkbox.setMargin( 10 );
    m_mesh_checkbox.setPosition( 0, 150 );
    m_mesh_checkbox.show();

    m_edge_checkbox.setMargin( 10 );
    m_edge_checkbox.setPosition( 0, 270 );
    m_edge_checkbox.show();
}

void Controller::show()
{
    m_volume_checkbox.show();
    m_mesh_checkbox.show();
    m_edge_checkbox.show();
}

void Controller::hide()
{
    m_volume_checkbox.hide();
    m_mesh_checkbox.hide();
    m_edge_checkbox.hide();
}

bool Controller::isShown()
{
    return
        m_volume_checkbox.isShown() ||
        m_mesh_checkbox.isShown() ||
        m_edge_checkbox.isShown();
}

} // end of namespace local
