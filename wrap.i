%module TMQGP

%include <std_map.i>
%include <std_string.i>
%include <std_unordered_map.i>
%include <std_vector.i>

%pythoncode %{
    import numpy as np
%}

%{
#define SWIG_FILE_WITH_INIT
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_API
%}
%include "numpy.i"
%init %{
import_array();
%}

%feature("autodoc", 1);

typedef std::string String;
%template(vec_type) std::vector<double>;
%template(dval) std::map<std::string, double>;
%template(content_type) std::map<std::string, std::map<std::string, double>>;
%template(str_vector) std::vector<std::string>;
%template(suppress_type) std::map<std::string, int>;
%template(eos_type) std::map<std::string, std::vector<double>>;
%template(mr_type) std::vector<std::map<std::string, double>>;
%template(integfunc_type) std::map<std::string, std::function<double (double, const double[])>>;


//%template(outfunc_type) std::map<std::string, std::vector<double>>;


%naturalvar Solver::variables;
%naturalvar Model::active_fields;
%naturalvar Model::charged_fields;
%naturalvar TOVsolver::iEoS;
%naturalvar TOVsolver::variables;

%{

//#include "stacktrace/stack_exception.hpp"
//#include "stacktrace/call_stack.hpp"
#include <functional>
#include <map>

#include "TMQGP/Interpolator.h"
#include "TMQGP/TMInterpolator_P.h"
// #include "TMQGP/TMtotColorProvider.h"
#include "TMQGP/test.h"
#include "TMQGP/Meson_Polarization.h"
#include "TMQGP/Meson_PhaseShift_and_Pressure.h"

%}


%include "std_function.i";
%include "std_complex.i"

%std_function(Functor, double, double);

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* x, int dimX)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* y, int dimY)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* z2, int dimZ1, int dimZ2)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* ReZ2, int dimZ1, int dimZ2)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* ImZ2, int dimZ3, int dimZ4)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* q, int dimQ)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* e, int dimE)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* pole, int dimPole)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* width, int dimWidth)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* v, int dimV)};
%include "TMQGP/Interpolator.h"
%include "TMQGP/TMInterpolator_P.h"
%include "TMQGP/test.h"
%include "TMQGP/Meson_Polarization.h"
%include "TMQGP/Meson_PhaseShift_and_Pressure.h"




%extend Interpolator {
%pythoncode {
    def __getstate__(self):
        x = list(self.x)
        data = list(self.data)
        return x, data, self.kind

    def __setstate__(self, state):
        x, data, kind = state
        self.__init__(np.array(x), np.array(data), kind)
}
}

%extend Interpolator2D {
%pythoncode {
    def __getstate__(self):
        x = list(self.x)
        y = list(self.y)
        z = list(self.z)
        return x, y, z

    def __setstate__(self, state):
        _x, _y, _z = state
        x = np.array(_x)
        y = np.array(_y)
        z = np.array(_z).reshape(len(y), len(x))
        self.__init__(x, y, z)
}
}


%extend InterDenom2D{
    %pythoncode{
        def __getstate__(self):
            x = list(self.x)
            y = list(self.y)
            z = list(self.z)
            z2 = list(self.z2)
            what = self.what
            return x, y, z, z2, what

        def __setstate__(self, state):
            _x, _y, _z, _z2, _what = state
            x = np.array(_x)
            y = np.array(_y)
            z = np.array(_z).reshape(len(y), len(x))
            z2 = np.array(_z2).reshape(len(y), len(x))
            what = _what
            self.__init__(x, y, z, z2, what)
    }
}


%extend RhoInterpolator {
%pythoncode {
    def __getstate__(self):
        x = list(self.x)
        y = list(self.y)
        z = list(self.z)
        z2 = list(self.z2)
        m = self.m
        mu = self.mu
        q = list(self.q)
        pole = list(self.pole)
        width = list(self.width)
        what = self.what

        return x, y, z, z2, m, mu, q, pole, width, what

    def __setstate__(self, state):
        _x, _y, _z, _z2, _m, _mu, _q, _pole, _width, _what = state
        x = np.array(_x)
        y = np.array(_y)
        z = np.array(_z).reshape(len(y), len(x))
        z2 = np.array(_z2).reshape(len(y), len(x))
        q = np.array(_q)
        pole = np.array(_pole)
        width = np.array(_width)
        m = _m
        mu = _mu
        what = _what
        self.__init__(x, y, z, z2, q, pole, width, what, m, mu)
}
}

%extend GFInterpolator {
%pythoncode {
    def __getstate__(self):
        x = list(self.x)
        y = list(self.y)
        z = list(self.z)
        z2 = list(self.z2)
        m = self.m
        mu = self.mu
        q = list(self.q)
        pole = list(self.pole)
        width = list(self.width)
        what = self.what

        return x, y, z, z2, m, mu, q, pole, width, what

    def __setstate__(self, state):
        _x, _y, _z, _z2, _m, _mu, _q, _pole, _width, _what = state
        x = np.array(_x)
        y = np.array(_y)
        z = np.array(_z).reshape(len(y), len(x))
        z2 = np.array(_z2).reshape(len(y), len(x))
        q = np.array(_q)
        pole = np.array(_pole)
        width = np.array(_width)
        m = _m
        mu = _mu
        what = _what

        self.__init__(x, y, z, z2, q, pole, width, what, m, mu)
    }
}


%extend PoleInterpolator{
%pythoncode{
    def __getstate__(self):
        x = list(self.x)
        y = list(self.y)
        z = list(self.z)
        z2 = list(self.z2)
        q = list(self.q)
        pole = list(self.pole)
        width = list(self.width)
        what = self.what

        return x, y, z, z2, q, pole, width, what

    def __setstate__(self, state):
        _x, _y, _z, _z2, _q, _pole, _width, _what = state
        x = np.array(_x)
        y = np.array(_y)
        z = np.array(_z).reshape(len(y), len(x))
        z2 = np.array(_z2).reshape(len(y), len(x))
        q = np.array(_q)
        pole = np.array(_pole)
        width = np.array(_width)
        what = _what

        self.__init__(x, y, z, z2, q, pole, width, what)
}
}

// %template(TMlist_type) std::vector<TMInterpolator *>;
%template(pair_dd) std::pair<double, double>;
// %naturalvar TM_tot_Interpolator::TMs;


