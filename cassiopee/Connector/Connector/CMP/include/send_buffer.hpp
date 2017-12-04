// SendBuffer.hpp
#ifndef _CMP_SENDBUFFER_HPP_
#define _CMP_SENDBUFFER_HPP_
#include <mpi.h>
#include <string>
#if __cplusplus > 199711L
#include <memory>
using std::shared_ptr;
#else
#include "shared_ptr.hpp"
using CMP::shared_ptr;
#endif
#include "vector_view.hpp"

namespace CMP {
    class SendBuffer {
    public:
        class PackedData {
        public:
            typedef char*       iterator;
            typedef const char* const_iterator;

            template <typename K>
            PackedData( const K& val )
                : m_nbBytes( sizeof( K ) ), m_pt_data( (char*)new K( val ) ), m_is_owner( true ) {}

            template <typename K>
            PackedData( const std::size_t nbElts, const K* pt_data )
                : m_nbBytes( nbElts * sizeof( K ) ), m_pt_data( (char*)pt_data ), m_is_owner( pt_data == NULL ) {
                if ( m_is_owner ) m_pt_data = (char*)new K[nbElts];
            }
            template <typename InpIterator>
            PackedData( InpIterator beg, InpIterator end ) : m_nbBytes( 0 ), m_pt_data( NULL ), m_is_owner( false ) {
                char* pt_beg = (char*)&( *beg );
                char* pt_end = (char*)&( *end );
                m_nbBytes    = pt_end - pt_beg;
                m_pt_data    = pt_beg;
            }

            PackedData( const std::string& str );

            /** Copie plus dans l'esprit du déplacement */
            PackedData( const PackedData& data )
                : m_nbBytes( data.m_nbBytes ), m_pt_data( data.m_pt_data ), m_is_owner( data.m_is_owner ) {
                data.m_is_owner = false;
            }
            /** Destructeur
            */
            ~PackedData( ) {
                if ( m_is_owner ) delete m_pt_data;
            }

            PackedData& operator=( PackedData& data ) {
                if ( this != &data ) {
                    if ( m_is_owner ) delete m_pt_data;
                    m_nbBytes       = data.m_nbBytes;
                    m_pt_data       = data.m_pt_data;
                    m_is_owner      = data.m_is_owner;
                    data.m_is_owner = false;
                }
                return *this;
            }

            /** Adresse des données à sérialiser */
            template <typename K>
            K* data( ) {
                return static_cast<K*>( m_pt_data );
            }
            /** Adresse des données à sérialiser
            */
            template <typename K>
            const K* data( ) const {
                return static_cast<K*>( m_pt_data );
            }
            std::size_t    size( ) const { return m_nbBytes; }  // Taille en octet des données sérialisés
            iterator       begin( ) { return static_cast<iterator>( m_pt_data ); }
            const_iterator begin( ) const { return static_cast<const_iterator>( m_pt_data ); }
            iterator       end( ) { return static_cast<iterator>( m_pt_data ) + size( ); }
            const_iterator end( ) const { return static_cast<const_iterator>( m_pt_data ) + size( ); }

        private:
            std::size_t  m_nbBytes;   // Taille en octet des données
            char*        m_pt_data;   // Adresse
            mutable bool m_is_owner;  // Données à détruire par la suite ?
        };
        // -----------------------------------------------------------------------------
        SendBuffer( int recv_rank, int id_tag = MPI_ANY_TAG, const MPI_Comm& comm = MPI_COMM_WORLD );
        SendBuffer( const SendBuffer& s_buf );
        ~SendBuffer( );

        SendBuffer& operator=( const SendBuffer& s_buf );

        SendBuffer& operator<<( const PackedData& data );
        template <typename K>
        SendBuffer& operator<<( const K& val ) {
            *this << PackedData( val );
            return *this;
        }
        template <typename K>
        SendBuffer& operator<<( const std::vector<K>& tab ) {
            *this << PackedData( tab.begin( ), tab.end( ) );
            return *this;
        }
        template <typename K>
        SendBuffer& operator<<( const vector_view<K>& tab ) {
            *this << PackedData( tab.begin( ), tab.end( ) );
            return *this;
        }

        int  isend( );  // Envoie asynchrone du buffer, retour un entier pour signaler une erreur éventuelle
        bool test( MPI_Status* pt_status = NULL );  // Renvoie vrai si l'envoi est fini.
        int wait( MPI_Status* pt_status = NULL );   // Attend que l'envoi se finisse
        int         receiver( ) const;              // Retourne le rang du receveur
        int         tag( ) const;                   // Retourne l'identifiant du message envoyé par le buffer
        std::size_t size( ) const;                  // Retourne la taille en octet pris par le buffer d'envoi
        const void* data( ) const;                  // Renvoie l'adresse du buffer d'envoi ( pour déboguage )

    private:  // ######################### PRIVATE PART BELOW #############################################
        // PIMPL template
        class Implementation;
        shared_ptr<Implementation> m_pt_implementation;
    };
    inline std::ostream& operator<<( std::ostream& out, const SendBuffer& s_buf ) {
        out << "SendBuffer{ receiver : " << s_buf.receiver( ) << ", id tag : " << s_buf.tag( )
            << ", size of buffer : " << s_buf.size( ) << ", addr buffer : " << s_buf.data( ) << " }";
        return out;
    }
}

#endif
