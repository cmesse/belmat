//
// Created by christian on 4/28/22.
//
#include <cmath>
#include "database.h"

namespace belmat
{
    database::database()
    {
        // do nothing for now
        // later, a lot of stuff will happen here
    }

//------------------------------------------------------------------------
    double
    database::jc( const double T, const double B, const double theta ) const
    {
        double Jc0 = 3.5e10/46.5 ;
        double alpha=1.2 ;
        double beta=1 ;
        double gamma=8.02 ;
        double Bc0=240e-3;

        double eps= std::sqrt(std::pow(std::sin(theta),2)/gamma
                    +std::pow(std::cos(theta),2) );

        return Jc0/std::pow(1+eps*std::pow(B/Bc0,beta),alpha);
    }
//------------------------------------------------------------------------
}