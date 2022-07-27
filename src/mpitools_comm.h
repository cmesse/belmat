//
// Created by christian on 5/26/22.
//

#ifndef BELMAT_MPITOOLS_COMM_H
#define BELMAT_MPITOOLS_COMM_H

#include <vector>
#include "mpitools_datatypes.h"

namespace mpitools
{

//------------------------------------------------------------------------------

    /**
     * return the rank of this proc
     */
    inline proc_t
    comm_rank()
    {
#ifdef MPI
        int rank ;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        return rank ;
#else
        return 0 ;
#endif
    }
//------------------------------------------------------------------------------

    /**
     * return the number of procs
     */
    inline proc_t
    comm_size()
    {
#ifdef MPI
        int size ;
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        return size ;
#else
        return 1 ;
#endif
    }

//------------------------------------------------------------------------------

    /**
     * wait for all procs to arrive here
     */
     inline void
     comm_barrier()
     {
#ifdef MPI
        MPI_Barrier( MPI_COMM_WORLD );
#endif
     }

//------------------------------------------------------------------------------

    /**
     * creates a unique communication tag
     */
    inline int
    comm_tag( const proc_t source, const proc_t target, const proc_t commsize=0 )
    {
#ifdef MPI
        // get the numbe rof procs
        proc_t size = commsize == 0 ? comm_size() : commsize ;

        // the smaller proc
        proc_t low = source < target ? source : target ;

        // the bigger proc
        proc_t high = source > target ? source : target ;

        return target == low ?
            65536 * ( high * size + low ) + 1 :
            65536 * ( high * size + low ) + 32769;
#else
        return -1 ;
#endif
    }

//------------------------------------------------------------------------------

    /**
     * this function is needed to split an mpi message into parts
     */
    inline std::vector< int >
    split_message( const std::size_t messagelength )
    {
        // number of full packages
        std::size_t numfull = messagelength / 65536 ;

        std::vector< int > steps( numfull + 1 );

        for( std::size_t k=0; k<numfull; ++k )
        {
            steps[ k ] = 65536 ;
        }
        steps[ numfull ] = messagelength - numfull * 65536 ;
        return steps ;
    }

//------------------------------------------------------------------------------

    /**
     * this function is used to distribute an array to other procs
     */
    template< typename T >
    void
    distribute( std::vector< T > & data )
    {
#ifdef MPI
        // the rank of the current proc
        proc_t myRank = comm_rank() ;

        // the number of procs
        proc_t numProcs = comm_size() ;

        // the datatype for the vector that is to be distributed
        mpi_t  dataType =  mpitools::mpi_type<T>();

        // type of length
        mpi_t sizeType = mpitools::mpi_type<std::size_t >();

        if( numProcs == 1 )
        {
            // do nothing since we have only one proc
            return ;
        }

        // get length of vector
        std::size_t dataSize = data.size() ;

        if( myRank == 0 )
        {
            // Allocate memory for status/request vector
            MPI_Status*  status1  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * numProcs );
            MPI_Request* request1 = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * numProcs );

            // loop over all other procs
            for( proc_t target=1; target<numProcs; ++target)
            {
                // create communication tag
                int tag = comm_tag( myRank, target );

                // send length to target
                MPI_Isend( &dataSize,
                           1,
                           sizeType,
                           target,
                           tag,
                           MPI_COMM_WORLD,
                           &request1[ target ] );

                // wait until send is complete
                MPI_Wait( &request1[ target ], &status1[ target ] );

                // continue if there is data to be sent
                if( dataSize > 0 )
                {
                    // calculate number of individual messages to be sent
                    std::vector< int > packSize = split_message( dataSize );

                    // number of messages
                    std::size_t numPacks = packSize.size() ;

                    // Allocate status and request containers
                    MPI_Status*  status2  = ( MPI_Status*  ) alloca( numPacks * sizeof( MPI_Status  ) );
                    MPI_Request* request2 = ( MPI_Request* ) alloca( numPacks * sizeof( MPI_Request ) );

                    // get raw pointer of vector
                    T * array = data.data() ;

                    // offset in data array
                    std::size_t offset = 0 ;

                    // loop over all messages
                    for( std::size_t m=0; m<numPacks; ++m )
                    {
                        // increment comm tag
                        ++tag;

                        // send data
                        MPI_Isend( &array[ offset ],
                                   packSize[ m ],
                                   dataType,
                                   target,
                                   tag,
                                   MPI_COMM_WORLD,
                                   &request2[ m ] );

                        // wait until send is complete
                        MPI_Wait( &request2[ m ], &status2[ m ] );

                        // increment offset
                        offset += packSize[ m ];
                    }
                } // end if dataSize > 0
            } // end loop over all procs
        } // end if master proc
        else
        {
            // Allocate memory for status/request vector
            MPI_Status  status1;
            MPI_Request request1;

            // create communication tag
            int tag = comm_tag( 0, myRank );

            // create communication tag
            // receive length from target
            MPI_Irecv( &dataSize,
                       1,
                       sizeType,
                       0,
                       tag,
                       MPI_COMM_WORLD,
                       &request1 );

            // wait until receive is complete
            MPI_Wait( &request1, &status1 );

            // continue if there is data to be sent
            if( dataSize > 0 )
            {
                // calculate number of individual messages to be sent
                std::vector<int> packSize = split_message(dataSize);

                // number of messages
                std::size_t numPacks = packSize.size();

                // Allocate status and request containers
                MPI_Status *status2 = (MPI_Status *) alloca(numPacks * sizeof(MPI_Status));
                MPI_Request *request2 = (MPI_Request *) alloca(numPacks * sizeof(MPI_Request));

                // allocate standard vector
                data.resize(dataSize);

                // get raw pointer of vector
                T *array = data.data();

                // offset in data array
                std::size_t offset = 0;

                // loop over all messages
                for (std::size_t m = 0; m < numPacks; ++m)
                {
                    // increment comm tag
                    ++tag;

                    // receive data
                    MPI_Irecv( &array[ offset ],
                               packSize[ m ],
                               dataType,
                               0,
                               tag,
                               MPI_COMM_WORLD,
                               &request2[ m ] );

                    // wait until send is complete
                    MPI_Wait( &request2[ m ], &status2[ m ] );

                    // increment offset
                    offset += packSize[ m ];
                } // end loop over all messages
            } // end if dataSize > 0
            else
            {
                // populate empty array
                data = {} ;
            }
        }
#endif
    }

//------------------------------------------------------------------------------
}

#endif //BELMAT_MPITOOLS_COMM_H
