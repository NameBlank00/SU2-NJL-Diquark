#ifndef __TEST_H
#define __TEST_H

#include "TMQGP/IntGSL.h"
#include "TMQGP/Integrator.h"

double integrate_x(double a, double b){
	return b*b/2 - a*a/2;
}

IntGSL<std::function<double(double)>> integ_k;

double integrate_x_num(double a, double b);

#endif