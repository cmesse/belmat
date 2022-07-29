#include <iostream>

#ifdef PETSC
#include <type_traits>  // for type checking is_same
#include <petscsys.h> // for PetscInt and PetscReal
#endif

#include "mpitools_comm.h"

#include "database.h"
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

    mpitools::init( &argc, &argv );

    if( mpitools::comm_rank() == 0 ) {
        std::cout << "Database Hello World!" << std::endl;
    }

    // this is the constructor for the database
    //belmat::database currentmap( "database.hdf5", "SultanDelta" );
    belmat::database heat( "materials.hdf5", "HastelloyCp" );
    //belmat::database lambda( "materials.hdf5", "CopperK_RRR300" );

    // the temperature in K
    double T = 20.0 ; // 7 ;

    // magnetic flux density in T
    double B = 15.0 ;; 6.11 ;

    // angle in rad
    double theta = 1.5708 ; //0.5236 ;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // now we compute the analytic value of the dummy function
    // this function will act as our blackbox.
    // now we must do the following things
    //
    // 1: find out how to link the libbelmat.a to SparseLizard
    // 2: find out how to translate this blackbox function into
    //    a SparseLizard expression
    //double jc = currentmap.compute( T, B, theta ) ;
    //double lambda_value = lambda.compute( T, B );
    double cp_value = heat.compute( T );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // print the result
    //std::cout << "proc " << mpitools::comm_rank() << ": T=" << T << " B=" << B << " theta=" << theta << " jc=" << jc << std::endl ;


    std::cout << "proc " << mpitools::comm_rank() << ": T=" << T << " B=" << B << " cp=" << cp_value << std::endl ;

    return mpitools::finanize() ;
}
