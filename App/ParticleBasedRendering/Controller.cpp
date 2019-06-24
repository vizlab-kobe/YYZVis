#include "Controller.h"


namespace local
{

Controller::Controller( local::Model* model, local::View* view ):
    m_view( view ),
    m_volume_checkbox( view, "Volume", "YinVolume", "YangVolume", "ZhongVolume" ),
    m_mesh_checkbox( view, "Mesh", "YinMesh", "YangMesh" ),
    m_edge_checkbox( view, "Edge", "YinEdge", "YangEdge", "ZhongEdge" )
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

} // end of namespace local
