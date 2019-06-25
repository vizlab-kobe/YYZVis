#include "Controller.h"


namespace local
{

Controller::Controller( local::Model* model, local::View* view ):
    m_view( view ),
    m_slice_checkbox( view, "Slice", "YinSlice", "YangSlice", "ZhongSlice" ),
    m_mesh_checkbox( view, "Mesh", "YinMesh", "YangMesh" ),
    m_edge_checkbox( view, "Edge", "YinEdge", "YangEdge", "ZhongEdge" )
{
    m_slice_checkbox.setMargin( 10 );
    m_slice_checkbox.setPosition( 0, 0 );
    m_slice_checkbox.show();

    m_mesh_checkbox.setMargin( 10 );
    m_mesh_checkbox.setPosition( 0, 150 );
    m_mesh_checkbox.show();

    m_edge_checkbox.setMargin( 10 );
    m_edge_checkbox.setPosition( 0, 270 );
    m_edge_checkbox.show();
}

} // end of namespace local
