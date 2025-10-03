#include "TMInterpolator_P.h"



XInterpolator_P::XInterpolator_P(Interpolator2D *ReX00, Interpolator2D *ReX01, Interpolator2D *ReX10, Interpolator2D *ReX11, Interpolator2D *ImX00, Interpolator2D *ImX01, Interpolator2D *ImX10, Interpolator2D *ImX11)
{
    this->ReX00 = ReX00;
    this->ReX01 = ReX01;
    this->ReX10 = ReX10;
    this->ReX11 = ReX11;

    this->ImX00 = ImX00;
    this->ImX01 = ImX01;
    this->ImX10 = ImX10;
    this->ImX11 = ImX11;
}

std::complex<double> XInterpolator_P::operator()(int i, int j, double P, double E)
{
    std::complex<double> res;

    double P_eval = P;
    double E_eval = E;

    if (this->old_P) {
        // If old_P is true, we use the old kinematics
        P_eval = 0;
        double sign_E = (E < 0) ? -1. : 1.;
        double Ecm2 = E * E - P * P;
        if (Ecm2 < 0) {
            return std::complex<double>(0, 0); // Avoid negative square root
        }
        E_eval = sign_E * sqrt(fabs(Ecm2));
    }

    if (i == 0 && j == 0){
        res = std::complex<double>(ReX00->operator()(P_eval, E_eval), ImX00->operator()(P_eval, E_eval));
    }
    else if (i == 0 && j == 1){
        res = std::complex<double>(ReX01->operator()(P_eval, E_eval), ImX01->operator()(P_eval, E_eval));
    }
    else if (i == 1 && j == 0){
        res = std::complex<double>(ReX10->operator()(P_eval, E_eval), ImX10->operator()(P_eval, E_eval));
    }
    else if (i == 1 && j == 1){
        res = std::complex<double>(ReX11->operator()(P_eval, E_eval), ImX11->operator()(P_eval, E_eval));
    }
    else{
        throw std::runtime_error("i,j not supported");
    }
    return res;
}

