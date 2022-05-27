#include <iostream>
#include "database.h"


int main() {
    std::cout << "Database Hello World!" << std::endl;

    // this is the constructor for the database
    // we don't need to pass anything yet
    // later, we will pass the filename to the HDF5 as a string
    // and the parameter from MPI_COMM_WORLD
    belmat::database database ;

    database.load_from_file( "database.hdf5", "SultanDelta" );

    // this dummy function just calls the analytic function so that we get some
    // numbers. Nico and I need to discuss if we want to introduce a material
    // subclass so that we can store different datasets in one database

    // the temperature in K
    double T = 7 ;

    // magnetic flux density in T
    double B = 6.11 ;

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
    double jc = database.jc( T, B, theta ) ;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // print the result

    std::cout << "demo: T=" << T << " B=" << B << " theta=" << theta << " jc=" << jc << std::endl ;

    return 0;
}
