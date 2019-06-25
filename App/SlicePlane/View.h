#pragma once
#include "Model.h"
#include <kvs/glut/Application>
#include <kvs/glut/Screen>


namespace local
{

/*===========================================================================*/
/**
 *  @brief  View class.
 */
/*===========================================================================*/
class View
{
private:
    local::Model* m_model; ///< pointer to the model
    kvs::glut::Screen m_screen; ///< screen

public:
    View( kvs::glut::Application* app, local::Model* model );

    kvs::glut::Screen& screen() { return m_screen; }
    void setup();
    void show();
    void redraw();

private:
    void setup_meshes();
    void setup_edges();
    void setup_slices();
};

} // end of namespace local
