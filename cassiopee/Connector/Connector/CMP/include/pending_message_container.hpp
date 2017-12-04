// pending_message_container.hpp
//////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef _CMP_PENDING_MESSAGE_CONTAINER_HPP_
#define _CMP_PENDING_MESSAGE_CONTAINER_HPP_
#include <vector>
#if __cplusplus > 199711L
#include <memory>
using std::shared_ptr;
#else
#include "shared_ptr.hpp"
using CMP::shared_ptr;
#endif

namespace CMP {
    class NullType {};
    template <typename MsgBuffer, typename LocalData = NullType>
    class PendingMsgContainer {
    public:
        typedef shared_ptr<MsgBuffer> pointer_message_buffer;
        typedef shared_ptr<LocalData> pointer_local_data;
        struct pending_message {
            pointer_message_buffer message_buffer;
            pointer_local_data     local_data;
            pending_message( pointer_message_buffer pt_msg, pointer_local_data pt_data )
                : message_buffer( pt_msg ), local_data( pt_data ) {}
            MsgBuffer&       get_message_buffer( ) { return *message_buffer; }
            LocalData&       get_local_data( ) { return *local_data; }
            const MsgBuffer& get_message_buffer( ) const { return *message_buffer; }
            const LocalData& get_local_data( ) const { return *local_data; }
        };
        typedef std::vector<pending_message>       container;
        typedef typename container::iterator       iterator;
        typedef typename container::const_iterator const_iterator;

        PendingMsgContainer( ) : m_pending_message( ) {}
        PendingMsgContainer( std::size_t maxSize ) : m_pending_message( ) { m_pending_message.reserve( maxSize ); }
        ~PendingMsgContainer( ) {}

        void push_back( MsgBuffer* pt_buff, LocalData* data = NULL );
        void pop( const iterator& it );

        bool        empty( ) const { return m_pending_message.empty( ); }
        std::size_t size( ) const { return m_pending_message.size( ); }
        void        clear( ) { m_pending_message.clear( ); }

        void waitAll( );
        // If iterator is equal to end(), not complete message found
        iterator get_first_complete_message( );

        pending_message& operator[]( std::size_t i );
        const pending_message& operator[]( std::size_t i ) const;

        pending_message&       back( ) { return m_pending_message.back( ); }
        const pending_message& back( ) const { return m_pending_message.back( ); }
        MsgBuffer&             back_message_buffer( ) { return *m_pending_message.back( ).message_buffer; }
        LocalData&             back_local_data( ) { return *m_pending_message.back( ).local_data; }
        const MsgBuffer&       back_message_buffer( ) const { return *m_pending_message.back( ).message_buffer; }
        const LocalData&       back_local_data( ) const { return *m_pending_message.back( ).local_data; }

        pending_message&       front( ) { return m_pending_message.front( ); }
        const pending_message& front( ) const { return m_pending_message.front( ); }
        MsgBuffer&             front_message_buffer( ) { return *m_pending_message.front( ).message_buffer; }
        LocalData&             front_local_data( ) { return *m_pending_message.front( ).local_data; }
        const MsgBuffer&       front_message_buffer( ) const { return *m_pending_message.front( ).message_buffer; }
        const LocalData&       front_local_data( ) const { return *m_pending_message.front( ).local_data; }

        iterator       begin( ) { return m_pending_message.begin( ); }
        const_iterator begin( ) const { return cbegin( ); }
        const_iterator cbegin( ) const { return m_pending_message.begin( ); }
        iterator       end( ) { return m_pending_message.end( ); }
        const_iterator end( ) const { return cend( ); }
        const_iterator cend( ) const { return m_pending_message.end( ); }

    private:
        container m_pending_message;
    };
}
#endif
