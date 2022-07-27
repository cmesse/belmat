//
// Created by christian on 5/25/22.
//

#ifndef BELMAT_MPITOOLS_DATATYPES_H
#define BELMAT_MPITOOLS_DATATYPES_H

#include <complex>

#ifdef MPI
#include <mpi.h>
#endif

namespace mpitools
{
#ifdef MPI
    typedef MPI_Datatype mpi_t ;
#else
    // define fake type if there is no MPI
    typedef int mpi_t;
#endif
    typedef int proc_t ;

//------------------------------------------------------------------------------

    /**
     *
     * /brief                 returns an MPI enum defining the
     *                        data type that is to be communicated.
     *
     * see also https://support.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
     */
    template<typename T> mpi_t
    mpi_type()
    {
        return T::mpi_type();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef MPI
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< int > ()
    {
        return MPI_INT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< long int > ()
    {
        return MPI_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< unsigned int > ()
    {
        return MPI_UNSIGNED;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< long unsigned int > ()
    {
        return MPI_UNSIGNED_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< double > ()
    {
        return MPI_DOUBLE ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline mpi_t
    mpi_type< long double > ()
    {
        return MPI_LONG_DOUBLE ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// MPI_CXX_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_DOUBLE_COMPLEX
    template<> inline mpi_t
    mpi_type< std::complex< double > > ()
    {
        return MPI_CXX_DOUBLE_COMPLEX ;
    }
#endif

// MPI_CXX_LONG_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    template<> inline mpi_t
    mpi_type< std::complex< long double > > ()
    {
        return MPI_CXX_LONG_DOUBLE_COMPLEX ;
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif // MPI
//------------------------------------------------------------------------------
}

#endif //BELMAT_MPITOOLS_DATATYPES_H
