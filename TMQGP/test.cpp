#include "test.h"

double integrate_x_num(double a, double b){
    funct f = [&](double x) -> double {
        return x;  // Using fixed y and z values
    };

    double res, err;

    integ_k.integrate(&f, a, b, res, err);

    return res;
}