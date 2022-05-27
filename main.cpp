#include <iostream>
#include "database.h"


int main() {
    std::cout << "Database Hello World!" << std::endl;

    // this is the constructor for the database
    belmat::database currentmap( "database.hdf5", "SultanDelta" );

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
    double jc = currentmap.compute( T, B, theta ) ;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // print the result

    std::cout << "demo: T=" << T << " B=" << B << " theta=" << theta << " jc=" << jc << std::endl ;

    return 0;
}
