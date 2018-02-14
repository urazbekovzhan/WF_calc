//
//  functions.hpp
//  WF_calc
//
//  Created by Bakytzhan on 28.01.2018.
//  Copyright © 2018 Bakytzhan. All rights reserved.
//

#ifndef FUNCTIONS_hpp
#define FUNCTIONS_hpp

#include <gsl/gsl_matrix.h>
#include "CLASSES.hpp"

double INTEGRAL_FERMDIST(int l, double a);   // Integral[ x^(2+2l)Exp(-a x^2 )  , x:0->Inf]

double COERRELATION_FUNC (double x, double y);

double COERRELATION_FUNC_OF2DSET ( double x, double y);

double A (WaveFunction& I, double j7, double j8, int i, int j);

double MATTER_DIST_SUBALPHA(double R, double aij, double bij, int lx, int ly);

double MatterDistSubNucleon(double y, double aij, double bij, int lx, int ly);


#endif /* functions_hpp */
