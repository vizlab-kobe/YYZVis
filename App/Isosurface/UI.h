#pragma once
#include "View.h"
#include "Model.h"
#include <kvs/CheckBox>
#include <kvs/CheckBoxGroup>
#include <kvs/Label>
#include <kvs/Slider>
#include <kvs/Font>
#include <kvs/FontMetrics>
#include <string>


namespace local
{

namespace UI
{

class CheckBox : public kvs::CheckBox
{
private:
    local::View* m_view; ///< pointer to the view
    std::string m_name; ///< object name

public:
    CheckBox( local::View* view, const std::string name ):
        kvs::CheckBox( &view->screen() ),
        m_view( view ),
        m_name( name )
    {
        setState( true );
    }

    const std::string& name() const { return m_name; }

    void stateChanged()
    {
        if ( m_view->screen().scene()->hasObject( m_name ) )
        {
            if ( state() ) { m_view->screen().scene()->object( m_name )->show(); }
            else { m_view->screen().scene()->object( m_name )->hide(); }

#if defined( JSST2019_TEST )
            if ( m_name == "YinIso" )
            {
                if ( state() ) { m_view->screen().scene()->object( "YinIso2" )->show(); }
                else { m_view->screen().scene()->object( "YinIso2" )->hide(); }

                if ( state() ) { m_view->screen().scene()->object( "YinIso3" )->show(); }
                else { m_view->screen().scene()->object( "YinIso3" )->hide(); }
            }

            if ( m_name == "YangIso" )
            {
                if ( state() ) { m_view->screen().scene()->object( "YangIso2" )->show(); }
                else { m_view->screen().scene()->object( "YangIso2" )->hide(); }

                if ( state() ) { m_view->screen().scene()->object( "YangIso3" )->show(); }
                else { m_view->screen().scene()->object( "YangIso3" )->hide(); }
            }

            if ( m_name == "ZhongIso" )
            {
                if ( state() ) { m_view->screen().scene()->object( "ZhongIso2" )->show(); }
                else { m_view->screen().scene()->object( "ZhongIso2" )->hide(); }

                if ( state() ) { m_view->screen().scene()->object( "ZhongIso3" )->show(); }
                else { m_view->screen().scene()->object( "ZhongIso3" )->hide(); }
            }
#endif
        }
    }
};

class CheckBoxGroup : public kvs::CheckBoxGroup
{
private:
    kvs::Label m_label;
    CheckBox m_yin_checkbox;
    CheckBox m_yang_checkbox;
    CheckBox m_zhong_checkbox;

public:
    CheckBoxGroup(
        local::View* view,
        const std::string& label_text,
        const std::string& yin_name,
        const std::string& yang_name,
        const std::string& zhong_name = "" ):
        kvs::CheckBoxGroup( &view->screen() ),
        m_label( &view->screen() ),
        m_yin_checkbox( view, yin_name ),
        m_yang_checkbox( view, yang_name ),
        m_zhong_checkbox( view, zhong_name )
    {
        m_label.setFont( kvs::Font( kvs::Font::Sans, kvs::Font::Bold, 24 ) );
        m_label.setText( label_text.c_str() );

        m_yin_checkbox.setFont( kvs::Font( kvs::Font::Sans, 22 ) );
        m_yin_checkbox.setCaption( "Yin" );

        m_yang_checkbox.setFont( kvs::Font( kvs::Font::Sans, 22 ) );
        m_yang_checkbox.setCaption( "Yang" );

        m_zhong_checkbox.setFont( kvs::Font( kvs::Font::Sans, 22 ) );
        m_zhong_checkbox.setCaption( "Zhong" );

        add( &m_yin_checkbox );
        add( &m_yang_checkbox );
        if ( !zhong_name.empty() ) { add( &m_zhong_checkbox ); }
    }

    void show()
    {
        kvs::CheckBoxGroup::show();
        m_label.show();
    }

    void hide()
    {
        kvs::CheckBoxGroup::hide();
        m_label.hide();
    }

    void setMargin( const int margin )
    {
        m_label.setMargin( margin );
        m_yin_checkbox.setMargin( margin );
        m_yang_checkbox.setMargin( margin );
        m_zhong_checkbox.setMargin( margin );
    }

    void screenUpdated()
    {
        const kvs::Vec2i p( this->x(), this->y() );
        m_label.setPosition( p.x(), p.y() );
        m_yin_checkbox.setPosition( p.x(), p.y() + m_label.height() - m_label.margin() );
        m_yang_checkbox.setPosition( p.x(), m_yin_checkbox.y() + m_yin_checkbox.height() - m_yin_checkbox.margin() );
        if ( !m_zhong_checkbox.name().empty() )
        {
            m_zhong_checkbox.setPosition( p.x(), m_yang_checkbox.y() + m_yang_checkbox.height() - m_yang_checkbox.margin() );
        }
    }
};

class Slider : public kvs::Slider
{
private:
    local::Model* m_model;
    local::View* m_view;
    kvs::Label m_label;

public:
    Slider( local::Model* model, local::View* view ):
        kvs::Slider( &view->screen() ),
        m_model( model ),
        m_view( view ),
        m_label( &view->screen() )
    {
        const float min_value = m_model->constYinVolume().minValue();
        const float max_value = m_model->constYinVolume().maxValue();
        kvs::Slider::setRange( min_value, max_value );
        kvs::Slider::setValue( m_model->isovalue() );
        kvs::Slider::setCaption( "" );
        kvs::Slider::setWidth( 180 );

        m_label.setFont( kvs::Font( kvs::Font::Sans, kvs::Font::Bold, 24 ) );
        m_label.setText( "Isovalue" );
    }

    void setMargin( const int margin )
    {
        kvs::Slider::setMargin( margin );
        m_label.setMargin( margin );
    }

    void show()
    {
        kvs::Slider::show();
        m_label.show();
    }

    void hide()
    {
        kvs::Slider::hide();
        m_label.hide();
    }

    void sliderReleased()
    {
        m_model->setIsovalue( this->value() );

        // Yin
        {
            kvs::PolygonObject* object = m_model->newYinIsosurfaces();
            object->setName( "YinIso" );
            m_view->screen().scene()->replaceObject( "YinIso", object );
        }

        // Yang
        {
            kvs::PolygonObject* object = m_model->newYangIsosurfaces();
            object->setName( "YangIso" );
            m_view->screen().scene()->replaceObject( "YangIso", object );
        }

        // Zhong
        {
            kvs::PolygonObject* object = m_model->newZhongIsosurfaces();
            object->setName( "ZhongIso" );
            m_view->screen().scene()->replaceObject( "ZhongIso", object );
        }
    }

    void screenUpdated()
    {
        const kvs::Vec2 p( m_view->screen().width() - this->width() - this->margin(), 0 );
        m_label.setPosition( p.x(), p.y() );
        kvs::Slider::setPosition( p.x(), p.y() + 8 );
    }
};

} // end of namespace UI

} // end of namespace local
