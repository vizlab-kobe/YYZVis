#include "Controller.h"


namespace local
{

Controller::Controller( local::Model* model, local::View* view ):
    m_view( view ),
    m_isosurface_checkbox( view, "Isosurace", "YinIso", "YangIso", "ZhongIso" ),
    m_mesh_checkbox( view, "Mesh", "YinMesh", "YangMesh" ),
    m_edge_checkbox( view, "Edge", "YinEdge", "YangEdge", "ZhongEdge" ),
    m_isovalue_slider( model, view )
{
    m_isosurface_checkbox.setMargin( 10 );
    m_isosurface_checkbox.setPosition( 0, 0 );
    m_isosurface_checkbox.show();

    m_mesh_checkbox.setMargin( 10 );
    m_mesh_checkbox.setPosition( 0, 150 );
    m_mesh_checkbox.show();

    m_edge_checkbox.setMargin( 10 );
    m_edge_checkbox.setPosition( 0, 270 );
    m_edge_checkbox.show();

    m_isovalue_slider.setMargin( 10 );
    m_isovalue_slider.setPosition( m_view->screen().width() - m_isovalue_slider.width() - 10, 0 );
    m_isovalue_slider.show();
}

} // end of namespace local
