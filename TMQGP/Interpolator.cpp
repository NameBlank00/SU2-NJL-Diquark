#include "Interpolator.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <cmath>
// #include<libInterpolate/Interpolate.hpp>



using namespace std;


Interpolator::Interpolator() {

}

Interpolator::Interpolator(double * x, int dimX, double * y, int dimY, string kind) {
	methods = { { "linear", gsl_interp_linear },{ "cubic", gsl_interp_cspline },
	{ "steffen", gsl_interp_steffen } };
	if (dimX != dimY) {
		throw;
	}
	this->kind = kind;
	if (methods.find(kind) != methods.end()) {
		interp = gsl_spline_alloc(methods[kind], dimX);
		accel = gsl_interp_accel_alloc();
		gsl_spline_init(interp, x, y, dimX);
	}
	else {
		throw;
	}

	this->x.clear();
	this->data.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
		this->data.push_back(y[i]);
	}

}

Interpolator::~Interpolator() {
	gsl_spline_free(interp);
	gsl_interp_accel_free(accel);
}

double Interpolator::operator()(double x) {
// 	cout << interp->x[0] << " ; " << interp->x[interp->size - 1] << endl;
// 	cout << interp->y[0] << " ; " << interp->y[interp->size - 1] << endl;
	double res = gsl_spline_eval(interp, x, accel);
	if (gsl_isnan(res)){
		return 0;
	}
	
	return res;
}


Interpolator2D::Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2){
	interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, z2, dimX, dimY);

	this->x.clear();
	this->y.clear();
	this->z.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(z2[i]);
	}

	debug = 0;
}

Interpolator2D::~Interpolator2D()
{
	gsl_spline2d_free(this->interp);
	gsl_interp_accel_free(this->accX);
	gsl_interp_accel_free(this->accY);


}

double Interpolator2D::operator()(double x, double y){
	// gsl_set_error_handler_off();
	double res;
	// try
	// {
	// 	res = gsl_spline2d_eval(interp, x, y, accX, accY);
	// }
	// catch(const std::exception& e)
	// {
	// 	std::cerr << e.what() << '\n';
	// 	res = 0;
	// }
	try{
		res = gsl_spline2d_eval(interp, x, y, accX, accY);
	}
	catch (exception){
		return 0.;
	};
	if (gsl_isnan(res)){
		return 0;
	}

	if (debug){
		printf("Interp2D: %.3e %.3e %.3e", x, y, res);
	}
	return res;
}

Interpolator2D_cubic::Interpolator2D_cubic(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2)
{
	interp = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, z2, dimX, dimY);

	this->x.clear();
	this->y.clear();
	this->z.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(z2[i]);
	}

	debug = 0;
}

InterDenom2D::InterDenom2D(double *x, int dimX, double *y, int dimY, double *ReZ2, 
			int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, string what){
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2, dimX, dimY);

	this->what = what;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2[i]);
	}
	debug = 0;
}

InterDenom2D::~InterDenom2D()
{
	// gsl_spline2d_free(this->interp);
	// gsl_interp_accel_free(this->accX);
	// gsl_interp_accel_free(this->accY);
	gsl_spline2d_free(this->iIm);
	gsl_interp_accel_free(this->accImX);
	gsl_interp_accel_free(this->accImY);
}

double InterDenom2D::real(double x, double y){
	double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	double im = gsl_spline2d_eval(iIm, x, y, accImX, accImY);
	double res = re / (re*re + im*im);
	
	if (gsl_isnan(res)){
		return 0;
	}
	return res;

}

double InterDenom2D::imag(double x, double y){
	double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	// cout << "re = " << re << endl;
	double im = gsl_spline2d_eval(iIm, x, y, accImX, accImY);
	// cout << "im = " << im << endl;
	double res = - im / (re*re + im*im);
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

double InterDenom2D::operator()(double x, double y){
	// cout << "this = " << this << endl;
	// cout << "what = " << this->what << endl;
	if (this->what == "real"){
		return this->real(x, y);
	}
	if (this->what == "imag"){
		return this->imag(x, y);
	}
	else{
		return -1;
	}
}

PoleInterpolator::PoleInterpolator(double *x, int dimX, double *y, int dimY,
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4,
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what) : 
			InterDenom2D(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, what)
			{
	iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}

PoleInterpolator::~PoleInterpolator(){
	delete iPole;
	delete iWidth;
}

PoleInterpolatorQ::PoleInterpolatorQ(double *x, int dimX, double *y, int dimY,
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4,
			double * e, int dimE, double * pole, int dimPole, double * width, int dimWidth, string what) : 
			InterDenom2D(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, what)
			{
	iPole = new Interpolator(e, dimE, pole, dimPole, "cubic");
	iWidth = new Interpolator(e, dimE, width, dimWidth, "cubic");
}

PoleInterpolatorQ::~PoleInterpolatorQ(){
	delete iPole;
	delete iWidth;
}

GFInterpolator::GFInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu)
{
	this->m = m;
	this->what = what;
	this->mu = mu;
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2, dimX, dimY);

	this->what = what;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2[i]);
	}

	for (int i = 0; i < dimQ; i++){
		this->q.push_back(q[i]);
	}

	for (int i = 0; i < dimPole; i++){
		this->pole.push_back(pole[i]);
	}

	for (int i = 0; i < dimWidth; i++){
		this->width.push_back(width[i]);
	}

	debug = 0;
	iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}

GFInterpolator::~GFInterpolator()
{
	
}

double GFInterpolator::real(double x, double y)
{
	return gsl_spline2d_eval(interp, x, y, accX, accY);
}

double GFInterpolator::imag(double x, double y)
{
	return gsl_spline2d_eval(iIm, x, y, accImX, accImY);
}

double GFInterpolator::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	double res;
	// complex<double> G = 1./(y - eps - complex<double>(this->real(x, y + this->mu), this->imag(x, y + this->mu)) 
	// 		+ this->mu);
	complex<double> G = 1./(y - eps - complex<double>(this->real(x, y), this->imag(x, y)) 
			+ this->mu);
			
	if (this->what == "real"){
		res = G.real();
	}
	else if (this->what == "imag"){
		res = G.imag();
	}
	else{
		res = -1;
	}
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}



RhoInterpolator::RhoInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m, double mu) : 
			GFInterpolator(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, q, dimQ, pole, dimPole, width, dimWidth, what, m, mu)
{
	
}

double RhoInterpolator::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	// double re = this->real(x, y + this->mu);
	double re = this->real(x, y);
	if (gsl_isnan(re)){
		re = 0;
	}

	// double im = this->imag(x, y + this->mu);
	double im = this->imag(x, y);
	if (gsl_isnan(im)){
		im = 0;
	}

	double res;
	complex<double> G = 1./(y - eps - complex<double>(re, im) + this->mu);
	
	return -G.imag() / M_PI;
}

// double TMR2Interpolator::operator()(double x, double y)
// {
//     return 0.0;
// }

// Interpolator2D_ML::Interpolator2D_ML(double *x, int dimX, double *y, int dimY, double *z2, int dimZ1, int dimZ2)
// {
// 	this->x.clear();
// 	this->y.clear();
// 	this->z.clear();

// 	for (int i = 0; i < dimX; i++){
// 		this->x.push_back(x[i]);
// 	}

// 	for (int i = 0; i < dimY; i++){
// 		this->y.push_back(y[i]);
// 	}

// 	for (int i = 0; i < dimY*dimX; i++){
// 		this->z.push_back(z2[i]);
// 	}


// 	debug = 0;
// }

GFInterpolator_NoPole::GFInterpolator_NoPole(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			string what, double m, double mu) :
			InterDenom2D(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, what)
{
	this->m = m;
	this->what = what;
	this->mu = mu;
}

double GFInterpolator_NoPole::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	double res;
	complex<double> G = 1./(y - eps - complex<double>(this->real(x, y), this->imag(x, y)) 
			+ this->mu);
			
	if (this->what == "real"){
		res = G.real();
	}
	else if (this->what == "imag"){
		res = G.imag();
	}
	else{
		res = -1;
	}
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

double GFInterpolator_NoPole::real(double x, double y){
	return gsl_spline2d_eval(interp, x, y, accX, accY);
}

double GFInterpolator_NoPole::imag(double x, double y){
	return gsl_spline2d_eval(iIm, x, y, accImX, accImY);
}





XInterpolator::XInterpolator(Interpolator * ReX00, Interpolator * ReX01, 
										Interpolator * ReX10, Interpolator * ReX11,
										Interpolator * ImX00, Interpolator * ImX01, 
										Interpolator * ImX10, Interpolator * ImX11)
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

std::complex<double> XInterpolator::operator()(int i, int j, double x)
{
    std::complex<double> res;
	if (i == 0 && j == 0){
		res = std::complex<double>(ReX00->operator()(x), ImX00->operator()(x));
	}
	else if (i == 0 && j == 1){
		res = std::complex<double>(ReX01->operator()(x), ImX01->operator()(x));
	}
	else if (i == 1 && j == 0){
		res = std::complex<double>(ReX10->operator()(x), ImX10->operator()(x));
	}
	else if (i == 1 && j == 1){
		res = std::complex<double>(ReX11->operator()(x), ImX11->operator()(x));
	}
	else{
		throw std::runtime_error("i,j not supported");
	}
		
	return res;
}


double xlog(double x){
	// return x;
	// return ((x > 0) - (x < 0)) * log(1 + fabs(x));
	return log(-x);
}

double ylog(double y){
	// return y;//
	// return ((y > 0) - (y < 0)) * (exp(fabs(y)) - 1);
	return -exp(y);
}

double xlog1(double x){
	return ((x > 0) - (x < 0)) * log(1 + fabs(x));
}

double ylog1(double y){
	return ((y > 0) - (y < 0)) * (exp(fabs(y)) - 1);
}


PoleInterpolatorLog::PoleInterpolatorLog(double *x, int dimX, double *y, int dimY, double *ReZ2, int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, double *q, int dimQ, double *pole, int dimPole, double *width, int dimWidth, string what)
{
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();

	double * ImZ2_log = new double[dimY*dimX];

	for (int i = 0; i < dimY*dimX; i++){
		// if (ImZ2[i] >= 0){
		// 	ImZ2[i] = -1e-50;
		// }
		ImZ2_log[i] = xlog1(ImZ2[i]);
	}

	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2_log, dimX, dimY);

	this->what = what;
	// cout << "In PoleInterpolatorLog: this = " << this << endl;
	// cout << "In PoleInterpolatorLog: what = " << this->what << endl;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2_log[i]);
	}

	debug = 0;

	for (int i = 0; i < dimQ; i++){
		this->q.push_back(q[i]);
	}

	for (int i = 0; i < dimPole; i++){
		this->pole.push_back(pole[i]);
	}

	for (int i = 0; i < dimWidth; i++){
		this->width.push_back(width[i]);
	}

	this->iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	this->iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}

double PoleInterpolatorLog::real(double x, double y)
{
    double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	double im = ylog1(gsl_spline2d_eval(iIm, x, y, accImX, accImY));
	double res = re / (re*re + im*im);
	
	if (gsl_isnan(res)){
		return 0;
	}
	return res;

}

double PoleInterpolatorLog::imag(double x, double y)
{
    double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	// cout << "re = " << re << endl;
	double im = ylog1(gsl_spline2d_eval(iIm, x, y, accImX, accImY));
	// cout << "im = " << im << endl;
	double res = - im / (re*re + im*im);
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

PoleInterpolatorLog::~PoleInterpolatorLog()
{

}

double PoleInterpolatorLog::operator()(double x, double y)
{
	// cout << "Log (): this = " << this << endl;
	// cout << "Log (): what = " << this->what << endl;
	if (this->what == "real"){
		return this->real(x, y);
	}
	if (this->what == "imag"){
		return this->imag(x, y);
	}
	else{
		return -1;
	}
}

Interpolator2DLog::Interpolator2DLog(double *x, int dimX, double *y, int dimY, double *z2, int dimZ1, int dimZ2)
{
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();

	double * ImZ2_log = new double[dimY*dimX];

	for (int i = 0; i < dimY*dimX; i++){
		if (z2[i] >= 0){
			z2[i] = -1e-50;
		}
		ImZ2_log[i] = xlog(z2[i]);
	}

	gsl_spline2d_init(interp, x, y, ImZ2_log, dimX, dimY);

	this->x.clear();
	this->y.clear();
	this->z.clear();


	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ImZ2_log[i]);
	}

	debug = 0;
}

double Interpolator2DLog::operator()(double x, double y)
{
    double res = ylog(gsl_spline2d_eval(interp, x, y, accX, accY));
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

Interpolator2DLog::~Interpolator2DLog()
{
	// gsl_spline2d_free(this->interp);
	// gsl_interp_accel_free(this->accX);
	// gsl_interp_accel_free(this->accY);
}

GFInterpolatorLog::GFInterpolatorLog(double *x, int dimX, double *y, int dimY, double *ReZ2, int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, double *q, int dimQ, double *pole, int dimPole, double *width, int dimWidth, string what, double m, double mu)
{
	gsl_set_error_handler_off();
	
	this->m = m;
	this->what = what;
	this->mu = mu;
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();

	double * ImZ2_log = new double[dimY*dimX];

	for (int i = 0; i < dimY*dimX; i++){
		if (ImZ2[i] >= 0){
			ImZ2[i] = -1e-50;
		}
		ImZ2_log[i] = xlog(ImZ2[i]);
	}

	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2_log, dimX, dimY);

	this->what = what;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2_log[i]);
	}

	for (int i = 0; i < dimQ; i++){
		this->q.push_back(q[i]);
	}

	for (int i = 0; i < dimPole; i++){
		this->pole.push_back(pole[i]);
	}

	for (int i = 0; i < dimWidth; i++){
		this->width.push_back(width[i]);
	}

	debug = 0;
	iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}


double GFInterpolatorLog::real(double x, double y)
{
	return gsl_spline2d_eval(interp, x, y, accX, accY);
}

double GFInterpolatorLog::imag(double x, double y)
{
	return ylog(gsl_spline2d_eval(iIm, x, y, accImX, accImY));
}

GFInterpolatorLog::~GFInterpolatorLog()
{
	// delete iPole;
	// delete iWidth;
}

RhoInterpolatorLog::RhoInterpolatorLog(double *x, int dimX, double *y, int dimY, double *ReZ2, int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, double *q, int dimQ, double *pole, int dimPole, double *width, int dimWidth, string what, double m, double mu)
: GFInterpolatorLog(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, q, dimQ, pole, dimPole, width, dimWidth, what, m, mu)
{
	
}

double RhoInterpolatorLog::operator()(double x, double y)
{
    double eps = sqrt(x*x + this->m*this->m);
	double re = this->real(x, y);
	if (gsl_isnan(re)){
		re = 0;
	}

	double im = this->imag(x, y);
	if (gsl_isnan(im)){
		im = 0;
	}

	double res;
	complex<double> G = 1./(y - eps - complex<double>(re, im) + this->mu);
	
	return -G.imag() / M_PI;
}

GFInterpolatorLogMu6::GFInterpolatorLogMu6(double *x, int dimX, double *y, int dimY, double *ReZ2, int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, double *q, int dimQ, double *pole, int dimPole, double *width, int dimWidth, string what, double m, double mu)
: GFInterpolatorLog(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, q, dimQ, pole, dimPole, width, dimWidth, what, m, mu)
{
	
}

double GFInterpolatorLogMu6::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	double re = this->real(x, y + this->mu);
	if (gsl_isnan(re)){
		re = 0;
	}

	double im = this->imag(x, y + this->mu);
	if (gsl_isnan(im)){
		im = 0;
	}

	double res;
	complex<double> G = 1./(y - eps - complex<double>(re, im) + this->mu);
	
	return -G.imag() / M_PI;
}

GFInterpolatorLogMu6::~GFInterpolatorLogMu6()
{
}

RhoInterpolatorLogMu6::RhoInterpolatorLogMu6(double *x, int dimX, double *y, int dimY, double *ReZ2, int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, double *q, int dimQ, double *pole, int dimPole, double *width, int dimWidth, string what, double m, double mu)
: RhoInterpolatorLog(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, q, dimQ, pole, dimPole, width, dimWidth, what, m, mu)
{
	
}

double RhoInterpolatorLogMu6::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	double re = this->real(x, y + this->mu);
	if (gsl_isnan(re)){
		re = 0;
	}

	double im = this->imag(x, y + this->mu);
	if (gsl_isnan(im)){
		im = 0;
	}

	double res;
	complex<double> G = 1./(y - eps - complex<double>(re, im) + this->mu);
	
	return -G.imag() / M_PI;
}

RhoInterpolatorLogMu6::~RhoInterpolatorLogMu6()
{
}

// Interpolator2D_SPLINTER::Interpolator2D_SPLINTER(double *x, int dimX, double *y, int dimY, double *z2, int dimZ1, int dimZ2)
// {

// 	DataTable data;
// 	DenseVector p(2);

// 	for (int i = 0; i < dimX; i++){
// 		for (int j = 0; j < dimY; j++){
// 			p(0) = x[i];
// 			p(1) = y[j];
// 			double f = z2[i*dimY + j];
// 			data.addSample(p, f); 
// 		}
// 	}

// 	BSpline spl = SPLINTER::BSpline::Builder(data).degree(1).build();
// }

// double Interpolator2D_SPLINTER::operator()(double x, double y)
// {
// 	DenseVector p(2);
// 	p(0) = x;
// 	p(1) = y;
//     return spline->eval(p);
// }
