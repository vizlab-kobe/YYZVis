#pragma once
#include "View.h"
#include "Model.h"
#include "UI.h"


namespace local
{

class Controller
{
private:
    local::View* m_view;
    local::UI::CheckBoxGroup m_isosurface_checkbox;
    local::UI::CheckBoxGroup m_mesh_checkbox;
    local::UI::CheckBoxGroup m_edge_checkbox;
    local::UI::Slider m_isovalue_slider;

public:
    Controller( local::Model* model, local::View* view );
};

} // end of namespace local
