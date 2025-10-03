#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "gsl/gsl_spline.h"
#include "gsl/gsl_spline2d.h"
#include <map>
#include <string>
#include <vector>
#include <complex>
// #include <mlinterp.hpp>
// #include <SPLINTER/bspline.h>
// #include <SPLINTER/bsplinebasis.h>
// #include <SPLINTER/datatable.h>


using namespace std;

/**
 * @class Interpolator
 * @brief A class for 1D interpolation using GSL splines.
 * 
 * This class provides a convenient wrapper around the GSL library's spline 
 * interpolation functionality. It allows you to create an interpolation object
 * from a set of data points and then evaluate the interpolated function 
 * at arbitrary points within the data range.
 */
class Interpolator {
public:
	/**
	 * @brief Constructor for the Interpolator class.
	 * 
	 * @param x Pointer to an array containing the x-coordinates of the data points.
	 * @param dimX The number of elements in the x array.
	 * @param y Pointer to an array containing the y-coordinates of the data points.
	 * @param dimY The number of elements in the y array.
	 * @param kind A string specifying the type of interpolation to use (e.g., "linear", "cubic").
	 */
	Interpolator(double * x, int dimX, double * y, int dimY, string kind);

	/**
	 * @brief Default constructor for the Interpolator class.
	 */
	Interpolator();

	/**
	 * @brief Destructor for the Interpolator class.
	 */
	~Interpolator();

	/// @brief Pointer to the GSL spline object.
	gsl_spline * interp;
	/// @brief Pointer to the GSL interpolation accelerator.
	gsl_interp_accel * accel;
	/// @brief Map of supported interpolation methods and their corresponding GSL types.
	map<string, const gsl_interp_type *> methods;

	/// @brief Vector storing the x-coordinates of the data points.
	vector<double> x;
	/// @brief Vector storing the y-coordinates of the data points.
	vector<double> data;
	/// @brief String storing the type of interpolation used.
	string kind;

	/**
	 * @brief Overloaded function call operator for evaluating the interpolated function.
	 * 
	 * @param x The point at which to evaluate the interpolated function.
	 * @return The interpolated value at the given point.
	 */
	double operator()(double x);

	/**
	 * @brief Calculates the derivative of the interpolated function at a given point.
	 * 
	 * @param x The point at which to calculate the derivative.
	 * @return The derivative of the interpolated function at the given point.
	 */
	double D(double x) {
		return gsl_spline_eval_deriv(interp, x, accel);
	}

	// template <class Archive> void  serialize(Archive & archive);
};


/**
 * @class Interpolator2D
 * @brief A base class for 2D interpolation.
 * 
 * This class provides a basic framework for 2D interpolation. 
 * It stores the interpolation data and provides a virtual function 
 * call operator for evaluating the interpolated function. 
 * Derived classes can implement different interpolation methods.
 */
class Interpolator2D{
	public:
		/**
		 * @brief Default constructor for the Interpolator2D class.
		 */
		Interpolator2D(){};
		/**
		 * @brief Constructor for the Interpolator2D class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param z2 Pointer to a 2D array containing the z-coordinates of the data points.
		 * @param dimZ1 The number of rows in the z2 array.
		 * @param dimZ2 The number of columns in the z2 array.
		 */
		Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);

		/**
		 * @brief Destructor for the Interpolator2D class.
		 */
		~Interpolator2D();

		/// @brief Pointer to the GSL 2D spline object.
		gsl_spline2d * interp;
		/// @brief Pointer to the GSL interpolation accelerator for the x-dimension.
		gsl_interp_accel * accX;
		/// @brief Pointer to the GSL interpolation accelerator for the y-dimension.
		gsl_interp_accel * accY;
		/// @brief Flag for debugging output.
		bool debug;

		/// @brief Vector storing the x-coordinates of the data points.
		vector<double> x;
		/// @brief Vector storing the y-coordinates of the data points.
		vector<double> y;
		/// @brief Vector storing the z-coordinates of the data points.
		vector<double> z;

		/**
		 * @brief Virtual function call operator for evaluating the interpolated function.
		 * 
		 * @param x The x-coordinate at which to evaluate the interpolated function.
		 * @param y The y-coordinate at which to evaluate the interpolated function.
		 * @return The interpolated value at the given point (x, y).
		 */
		virtual double operator()(double x, double y);
};


class Interpolator2DLog : public Interpolator2D{
	public:
		Interpolator2DLog(){};
		Interpolator2DLog(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);

		double operator()(double x, double y) override;
		~Interpolator2DLog();
};
// class Interpola




/**
 * @class Interpolator2D_cubic
 * @brief A class for 2D cubic interpolation.
 * 
 * This class inherits from Interpolator2D and implements 
 * bicubic interpolation using GSL splines.
 */
class Interpolator2D_cubic : public Interpolator2D{
	public:
		/**
		 * @brief Default constructor for the Interpolator2D_cubic class.
		 */
		Interpolator2D_cubic(){};
		/**
		 * @brief Constructor for the Interpolator2D_cubic class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param z2 Pointer to a 2D array containing the z-coordinates of the data points.
		 * @param dimZ1 The number of rows in the z2 array.
		 * @param dimZ2 The number of columns in the z2 array.
		 */
		Interpolator2D_cubic(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);
};


/**
 * @class InterDenom2D
 * @brief A class for 2D interpolation of a complex function represented as a fraction.
 * 
 * This class inherits from Interpolator2D and provides functionality for 
 * interpolating both the real and imaginary parts of a complex function 
 * that is represented as a fraction, where the denominator is the same 
 * for both the real and imaginary parts.
 */
class InterDenom2D : public Interpolator2D{
	public:
		/**
		 * @brief Default constructor for the InterDenom2D class.
		 */
		InterDenom2D(){};
		/**
		 * @brief Constructor for the InterDenom2D class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 */
		InterDenom2D(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, 
			double * ImZ2, int dimZ3, int dimZ4, string what);

		/**
		 * @brief Destructor for the InterDenom2D class.
		 */
		~InterDenom2D();

		/// @brief Pointer to the GSL interpolation accelerator for the x-dimension (imaginary part).
		gsl_interp_accel * accImX;
		/// @brief Pointer to the GSL interpolation accelerator for the y-dimension (imaginary part).
		gsl_interp_accel * accImY;

		/// @brief Pointer to the GSL 2D spline object for the imaginary part.
		gsl_spline2d * iIm;

		/// @brief Vector storing the imaginary part of the numerator.
		vector<double> z2;
		/// @brief String indicating whether to return the real or imaginary part.
		string what;

		/**
		 * @brief Evaluates the real part of the interpolated function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The real part of the interpolated function at (x, y).
		 */
		virtual double real(double x, double y);
		/**
		 * @brief Evaluates the imaginary part of the interpolated function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The imaginary part of the interpolated function at (x, y).
		 */
		virtual double imag(double x, double y);

		/**
		 * @brief Overridden function call operator for evaluating the interpolated function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The real or imaginary part of the interpolated function at (x, y), depending on the value of 'what'.
		 */
		double operator()(double x, double y) override;
};

// ########## Write the integrations here temporarily


/**
 * @class PoleInterpolator
 * @brief A class for 2D interpolation of a complex function with a pole, represented as a fraction.
 *
 * This class inherits from InterDenom2D and adds functionality for interpolating 
 * the position and width of a pole in the complex plane. This is useful for 
 * representing propagators with a spectral function that has a peak (the pole).
 */
class PoleInterpolator: public InterDenom2D{
	public:
		/**
		 * @brief Default constructor for the PoleInterpolator class.
		 */
		PoleInterpolator(){};
		/**
		 * @brief Constructor for the PoleInterpolator class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param q Pointer to an array of momentum values for the pole and width.
		 * @param dimQ The number of elements in the q array.
		 * @param pole Pointer to an array of pole positions for each momentum value.
		 * @param dimPole The number of elements in the pole array.
		 * @param width Pointer to an array of pole widths for each momentum value.
		 * @param dimWidth The number of elements in the width array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 */
		PoleInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what);

		/// @brief Interpolator object for the pole position as a function of momentum.
		Interpolator * iPole;
		/// @brief Interpolator object for the pole width as a function of momentum.
		Interpolator * iWidth;

		/// @brief Vector storing the momentum values for the pole and width.
		vector<double> q;
		/// @brief Vector storing the pole positions for each momentum value.
		vector<double> pole;
		/// @brief Vector storing the pole widths for each momentum value.
		vector<double> width;

		/// @brief String indicating whether to return the real or imaginary part.
		string what;

		/**
		 * @brief Destructor for the PoleInterpolator class.
		 */
		~PoleInterpolator();
};

class PoleInterpolatorLog: public PoleInterpolator{
	public:
		PoleInterpolatorLog(){};
		PoleInterpolatorLog(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what);

		double real(double x, double y) override;
		double imag(double x, double y) override;
		~PoleInterpolatorLog();

		double operator()(double x, double y) override;
};



class PoleInterpolatorQ: public InterDenom2D{
	public:
		/**
		 * @brief Default constructor for the PoleInterpolatorQ class.
		 */
		PoleInterpolatorQ(){};
		/**
		 * @brief Constructor for the PoleInterpolatorQ class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param e Pointer to an array of energy values for the pole and width.
		 * @param dimE The number of elements in the e array.
		 * @param pole Pointer to an array of pole positions for each energy value.
		 * @param dimPole The number of elements in the pole array.
		 * @param width Pointer to an array of pole widths for each energy value.
		 * @param dimWidth The number of elements in the width array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 */
		PoleInterpolatorQ(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * e, int dimE, double * pole, int dimPole, double * width, int dimWidth, string what);

		/// @brief Interpolator object for the pole position as a function of energy.
		Interpolator * iPole;
		/// @brief Interpolator object for the pole width as a function of energy.
		Interpolator * iWidth;

		/// @brief Vector storing the energy values for the pole and width.
		vector<double> q;
		/// @brief Vector storing the pole positions for each energy value.
		vector<double> pole;
		/// @brief Vector storing the pole widths for each energy value.
		vector<double> width;
		/// @brief String indicating whether to return the real or imaginary part.
		string what;
		
		/**
		 * @brief Destructor for the PoleInterpolatorQ class.
		 */
		~PoleInterpolatorQ();
};



// Interpolating Sigma and return G on evaluation
/**
 * @class GFInterpolator
 * @brief A class for 2D interpolation of a complex Green's function with a pole.
 * 
 * This class inherits from PoleInterpolator and provides functionality for 
 * interpolating the real and imaginary parts of a complex Green's function, 
 * taking into account the presence of a pole. It assumes the Green's function 
 * is represented as a fraction, where the denominator includes the pole structure.
 */
class GFInterpolator: public PoleInterpolator{
	public:
		/**
		 * @brief Default constructor for the GFInterpolator class.
		 */
		GFInterpolator(){};
		/**
		 * @brief Constructor for the GFInterpolator class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param q Pointer to an array of momentum values for the pole and width.
		 * @param dimQ The number of elements in the q array.
		 * @param pole Pointer to an array of pole positions for each momentum value.
		 * @param dimPole The number of elements in the pole array.
		 * @param width Pointer to an array of pole widths for each momentum value.
		 * @param dimWidth The number of elements in the width array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 * @param m The mass of the particle.
		 * @param mu The chemical potential.
		 */
		GFInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		/**
		 * @brief Destructor for the GFInterpolator class.
		 */
		~GFInterpolator();

		/// @brief The mass of the particle.
		double m;
		/// @brief The chemical potential.
		double mu;

		/**
		 * @brief Evaluates the real part of the interpolated Green's function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The real part of the interpolated Green's function at (x, y).
		 */
		double real(double x, double y) override;
		/**
		 * @brief Evaluates the imaginary part of the interpolated Green's function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The imaginary part of the interpolated Green's function at (x, y).
		 */
		double imag(double x, double y) override;
		/**
		 * @brief Overridden function call operator for evaluating the interpolated Green's function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The real or imaginary part of the interpolated Green's function at (x, y), depending on the value of 'what'.
		 */
		double operator()(double x, double y) override;
};

class GFInterpolatorLog: public GFInterpolator{
	public:
		GFInterpolatorLog(){};
		GFInterpolatorLog(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		double real(double x, double y) override;
		double imag(double x, double y) override;
		~GFInterpolatorLog();

		// double operator()(double x, double y) override;
};

class GFInterpolatorLogMu6: public GFInterpolatorLog{
	public:
		GFInterpolatorLogMu6(){};
		GFInterpolatorLogMu6(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		double operator()(double x, double y) override;
		~GFInterpolatorLogMu6();

		// double operator()(double x, double y) override;
};


/**
 * @class GFInterpolator_NoPole
 * @brief A class for 2D interpolation of a complex Green's function without a pole.
 * 
 * This class inherits from InterDenom2D and provides functionality for 
 * interpolating the real and imaginary parts of a complex Green's function, 
 * without taking into account the presence of a pole. It assumes the Green's function 
 * is represented as a fraction, where the denominator does not include the pole structure.
 */
class GFInterpolator_NoPole: public InterDenom2D{
	public:
		/**
		 * @brief Default constructor for the GFInterpolator_NoPole class.
		 */
		GFInterpolator_NoPole(){};
		/**
		 * @brief Constructor for the GFInterpolator_NoPole class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 * @param m The mass of the particle.
		 * @param mu The chemical potential.
		 */
		GFInterpolator_NoPole(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			string what, double m, double mu);

		/**
		 * @brief Destructor for the GFInterpolator_NoPole class.
		 */
		~GFInterpolator_NoPole(){};

		/// @brief The mass of the particle.
		double m;
		/// @brief The chemical potential.
		double mu;
		/**
		 * @brief Overridden function call operator for evaluating the interpolated Green's function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The real or imaginary part of the interpolated Green's function at (x, y), depending on the value of 'what'.
		 */
		double operator()(double x, double y) override;
		double real(double x, double y) override;
		double imag(double x, double y) override;
};



/**
 * @class RhoInterpolator
 * @brief A class for 2D interpolation of the spectral function.
 * 
 * This class inherits from GFInterpolator and provides functionality for 
 * interpolating the spectral function, which is related to the imaginary part 
 * of the Green's function.
 */
class RhoInterpolator: public GFInterpolator{
	public:
		/**
		 * @brief Default constructor for the RhoInterpolator class.
		 */
		RhoInterpolator(){};
		/**
		 * @brief Constructor for the RhoInterpolator class.
		 * 
		 * @param x Pointer to an array containing the x-coordinates of the data points.
		 * @param dimX The number of elements in the x array.
		 * @param y Pointer to an array containing the y-coordinates of the data points.
		 * @param dimY The number of elements in the y array.
		 * @param ReZ2 Pointer to a 2D array containing the real part of the numerator.
		 * @param dimZ1 The number of rows in the ReZ2 array.
		 * @param dimZ2 The number of columns in the ReZ2 array.
		 * @param ImZ2 Pointer to a 2D array containing the imaginary part of the numerator.
		 * @param dimZ3 The number of rows in the ImZ2 array.
		 * @param dimZ4 The number of columns in the ImZ2 array.
		 * @param q Pointer to an array of momentum values for the pole and width.
		 * @param dimQ The number of elements in the q array.
		 * @param pole Pointer to an array of pole positions for each momentum value.
		 * @param dimPole The number of elements in the pole array.
		 * @param width Pointer to an array of pole widths for each momentum value.
		 * @param dimWidth The number of elements in the width array.
		 * @param what A string indicating whether to return the real or imaginary part ("real" or "imag").
		 * @param m The mass of the particle.
		 * @param mu The chemical potential.
		 */
		RhoInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		/**
		 * @brief Overridden function call operator for evaluating the interpolated spectral function.
		 * 
		 * @param x The x-coordinate at which to evaluate the function.
		 * @param y The y-coordinate at which to evaluate the function.
		 * @return The interpolated value of the spectral function at (x, y).
		 */
		double operator()(double x, double y) override;
};

class RhoInterpolatorLog : public GFInterpolatorLog{
	public:
		RhoInterpolatorLog(){};
		RhoInterpolatorLog(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		double operator()(double x, double y) override;
		~RhoInterpolatorLog(){};
};

class RhoInterpolatorLogMu6 : public RhoInterpolatorLog{
	public:
		RhoInterpolatorLogMu6(){};
		RhoInterpolatorLogMu6(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu);

		double operator()(double x, double y) override;
		~RhoInterpolatorLogMu6();
};

class XInterpolator{
	public:
		XInterpolator(){};

		XInterpolator(Interpolator * ReX00, Interpolator * ReX01, 
					Interpolator * ReX10, Interpolator * ReX11,
					Interpolator * ImX00, Interpolator * ImX01, 
					Interpolator * ImX10, Interpolator * ImX11);

		Interpolator * ReX00;
		Interpolator * ReX01;
		Interpolator * ReX10;
		Interpolator * ReX11;

		Interpolator * ImX00;
		Interpolator * ImX01;
		Interpolator * ImX10;
		Interpolator * ImX11;


		~XInterpolator(){};
		std::complex<double> operator()(int i, int j, double x);
};


#endif