#pragma once
#include "View.h"
#include "Model.h"
#include "UI.h"
#include <kvs/KeyPressEventListener>


namespace local
{

class Controller
{
    class KeyPressEvent : public kvs::KeyPressEventListener
    {
    private:
        local::Controller* m_controller;
    public:
        KeyPressEvent( local::Controller* controller );
        void update( kvs::KeyEvent* event );
    };

private:
    local::View* m_view;
    local::UI::CheckBoxGroup m_slice_checkbox;
    local::UI::CheckBoxGroup m_mesh_checkbox;
    local::UI::CheckBoxGroup m_edge_checkbox;
    KeyPressEvent m_key_event;

public:
    Controller( local::Model* model, local::View* view );

protected:
    local::View* view() { return m_view; }
    void show();
    void hide();
    bool isShown();
};

} // end of namespace local
