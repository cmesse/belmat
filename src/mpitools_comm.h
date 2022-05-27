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
        MPI_Comm_size( MPI_COMM_WORLD, size );
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
            8192 * ( high * size + low ) + 1 :
            8192 * ( high * size + low ) + 4097;
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

//------------------------------------------------------------------------------
}

#endif //BELMAT_MPITOOLS_COMM_H
