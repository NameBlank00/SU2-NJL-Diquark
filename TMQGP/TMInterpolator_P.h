#ifndef _TMINTERPOLATOR_P_H
#define _TMINTERPOLATOR_P_H

#include <complex>
#include "Interpolator.h"

class XInterpolator_P{
	public:
		XInterpolator_P(){};

		XInterpolator_P(Interpolator2D * ReX00, Interpolator2D * ReX01, 
					Interpolator2D * ReX10, Interpolator2D * ReX11,
					Interpolator2D * ImX00, Interpolator2D * ImX01, 
					Interpolator2D * ImX10, Interpolator2D * ImX11);

		Interpolator2D * ReX00;
		Interpolator2D * ReX01;
		Interpolator2D * ReX10;
		Interpolator2D * ReX11;

		Interpolator2D * ImX00;
		Interpolator2D * ImX01;
		Interpolator2D * ImX10;
		Interpolator2D * ImX11;

        
		~XInterpolator_P(){};

		bool old_P = false;
		std::complex<double> operator()(int i, int j, double P, double E);
};

#endif