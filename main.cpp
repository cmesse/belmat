#include <iostream>

#ifdef PETSC
#include <type_traits>  // for type checking is_same
#include <petscsys.h> // for PetscInt and PetscReal
#endif

#include "mpitools_comm.h"

#include "material.h"
namespace mpitools
{
//------------------------------------------------------------------------------

    int
    init( int    *argc,
              char ***argv )
    {
#ifdef PETSC
        return PetscInitialize( argc, argv, NULL, NULL );
#elif MPI
        return  MPI_Init( argc, argv );
#else
        return 0 ;
#endif
    }

//------------------------------------------------------------------------------

    int
    finanize()
    {
#ifdef PETSC
        MPI_Barrier( MPI_COMM_WORLD );
        return PetscFinalize();
#elif MPI
        MPI_Barrier( MPI_COMM_WORLD );
        return MPI_Finalize();
#endif
        return 0 ;
    }
//------------------------------------------------------------------------
}

int main( int    argc,
         char * argv[] )
{
    // initialize MPI
    mpitools::init( &argc, &argv );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /*
     * Material selection.
     */

    // we assume copper rrr=50 as default
    std::string database = "materials.hdf5" ;
    std::string label = "copper" ;
    uint rrr = 0 ;

    // otherwise, we take the values from the console
    if( argc > 1 )
    {
        label = std::string( argv[ 1 ] );
    }
    if( argc > 2 )
    {
        rrr = std::stoi( std::string( argv[ 2 ] ) );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // create the material
    belmat::material mat( database, label, rrr );

    if( mpitools::comm_rank() == 0 )
    {
        double T = 0.0 ;
        std::cout << mat.label() << std::endl ;
        std::cout << "    density @298.15K: " << mat.density() << std::endl ;
        while( T < 300 )
        {
            std::cout << "    " << T ;

            if( mat.has_cp() )
            {
                std::cout << " " << mat.cp( T ) ;
            }
            if( mat.has_k() )
            {
                std::cout << " " << mat.k( T ) ;
            }
            if( mat.has_rho() )
            {
                std::cout << " " << mat.rho( T ) ;
            }
            std::cout << std::endl ;
            T += 5.0 ;
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    return mpitools::finanize() ;
}
