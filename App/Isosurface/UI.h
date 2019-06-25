#pragma once
#include "View.h"
#include <kvs/CheckBox>
#include <kvs/CheckBoxGroup>
#include <kvs/Label>
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

} // end of namespace UI

} // end of namespace local
