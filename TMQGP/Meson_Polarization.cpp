#include "Meson_Polarization.h"

#include <functional>

#include <cmath>
#include <limits>

#include <gsl/gsl_errno.h>

double heaviside(double x) {
    return (x >= 0.0) ? 1.0 : 0.0;
}




// Fermi-Dirac distribution function
double n(double E, double T) {
    if (T < 1e-12) return (E <= 0) ? 1.0 : 0.0;
    double x = E / T;
    if (x > 700) return 0.0;
    if (x < -700) return 1.0;
    if (x >= 0.0) {
        double ex = std::exp(x);
        return 1.0 / (1.0 + ex);
    } else {
        double exm = std::exp(-x);
        return exm / (1.0 + exm);
    }
}







double ImPi_Pion(double w, double q, double T, double p, double m, double mu_star, double Delta, double thermal_cutoff, double Nf, double Lambda_CutOff) {

    const double mu_star_eff = p * mu_star;  // use local const

    double s = w * w - q * q;

    double P = Nf * 3 / (8.0 * pi);

    // Set thermal cutoff based on thermal_cutoff value (0 for inf, 1 for 651)
    double L = (thermal_cutoff == 0) ? std::numeric_limits<double>::infinity() : 651.0;
    
    // Modify this value based on the value of the thermal cutoff
    const double lower_cutoff = s - 4.0 * m * m;
    const double upper_cutoff = 4.0 * (Lambda_CutOff * Lambda_CutOff + m * m) - s;
    const double upper_thermal_cutoff = 4.0 * (L*L + m * m) - s;

    double Jp = 0.0;

    if (lower_cutoff >= 0) {
        if (q == 0.0) {
            // Avoid repeated sqrt calculation
            double root_term = std::sqrt(std::max(0.0, 1.0 - 4.0 * m * m / s));
            Jp = root_term * (1 - n( (w - 2*mu_star) / 2.0, T) - n( (w + 2*mu_star) / 2.0, T) ) * heaviside(upper_thermal_cutoff);
        } else {
            double ep = w / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / s);
            double em = w / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / s);
            if (T != 0) {
                Jp = 0.5 * T * std::log((n(-em - mu_star, T) * n(em - mu_star, T)) / (n(-ep - mu_star, T) * n(ep - mu_star, T))) * heaviside(upper_thermal_cutoff);
                Jp += 0.5 * T * std::log((n(-em + mu_star, T) * n(em + mu_star, T)) / (n(-ep + mu_star, T) * n(ep + mu_star, T))) * heaviside(upper_thermal_cutoff);
            //Jp = T * std::log( (1 + cosh((w+q)/(2*T)))/(1 + cosh((w-q)/(2*T))) )*heaviside(upper_cutoff);
            } else {
                Jp = (std::abs(ep - mu_star) - std::abs(em - mu_star)) * heaviside(upper_thermal_cutoff);
            }
        }
    } 
    


    

    
    if (q != 0) {
        P *= 1 / q;

        if (s < 0) {
        double ep = w / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / s);
        double em = w / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / s);
        double a = n(em + mu_star, T) * n(em - mu_star, T);
        double b = n(-ep + mu_star, T) * n(-ep - mu_star, T);

        if (abs(a-b) < 1e-7) {
            Jp += 0;
        } else{
            Jp += 0.5*( T * log( a/b ) );
        }


        // double a = n(em + mu_star, T);
        // double b = n(em - mu_star, T);
        // double c = n(-ep + mu_star, T);
        // double d = n(-ep - mu_star, T);

        // if (abs(a*b - c*d) < 1e-10) {
        //     Jp += 0;
        // } else{
        //     Jp += ( T * log( a*b/(c*d) ) );
        // }
    } else{
        Jp += 0;
    }
    
    }

    return P * (s * Jp);
}





double ImPi_Diquark(double w, double q, double T, double p, double m, double mu_star, double Delta, double thermal_cutoff, double Nf, double Lambda_CutOff) {

    const double mu_star_eff = p * mu_star;  // use local const

    double wp = w + 2.0 * mu_star_eff;
    double sp = wp * wp - q * q;

    double P = Nf / (4.0 * pi);

    // Set thermal cutoff based on thermal_cutoff value (0 for inf, 1 for 651)
    double L = (thermal_cutoff == 0) ? std::numeric_limits<double>::infinity() : 651.0;
    
    // Modify this value based on the value of the thermal cutoff
    const double lower_cutoff = sp - 4.0 * m * m;
    const double upper_cutoff = 4.0 * (Lambda_CutOff * Lambda_CutOff + m * m) - sp;
    const double upper_thermal_cutoff = 4.0 * (L*L + m * m) - sp;

    double Jp = 0.0;

    if (lower_cutoff >= 0) {
        if (q == 0.0) {
            // Avoid repeated sqrt calculation
            double root_term = std::sqrt(std::max(0.0, 1.0 - 4.0 * m * m / sp));
            Jp = root_term * (1 - 2.0 * n(w / 2.0, T)) * heaviside(upper_thermal_cutoff);
        } else {
            double ep = wp / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
            double em = wp / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
            if (T != 0) {
                Jp = T * std::log((n(-em + mu_star, T) * n(em - mu_star, T)) / (n(-ep + mu_star, T) * n(ep - mu_star, T))) * heaviside(upper_thermal_cutoff);
            //Jp = T * std::log( (1 + cosh((w+q)/(2*T)))/(1 + cosh((w-q)/(2*T))) )*heaviside(upper_cutoff);
            } else {
                Jp = (std::abs(ep - mu_star) - std::abs(em - mu_star)) * heaviside(upper_thermal_cutoff);
            }
        }
    } 
    


    

    
    if (q != 0) {
        P *= 1 / q;

        if (sp < 0) {
        double ep = wp / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
        double em = wp / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
        double a = n(em - mu_star, T);
        double b = n(-ep + mu_star, T);

        if (abs(a-b) < 1e-7) {
            Jp += 0;
        } else{
            Jp += ( T * log( a/b ) );
        }


        // double a = n(em + mu_star, T);
        // double b = n(em - mu_star, T);
        // double c = n(-ep + mu_star, T);
        // double d = n(-ep - mu_star, T);

        // if (abs(a*b - c*d) < 1e-10) {
        //     Jp += 0;
        // } else{
        //     Jp += ( T * log( a*b/(c*d) ) );
        // }
    } else{
        Jp += 0;
    }
    
    }

    return P * (sp * Jp);
}



struct RegularFunctionParams {
    double q;
    double T;
    double p;
    double m;
    double mu_star;
    double delta;
    double Nf;
    double Lambda_CutOff;
};

// Cauchy integration (unchanged)
template <typename Func>
double integrate_cauchy_lambda(Func* f, double lower, double upper, double c, gsl_integration_workspace* workspace) {
    struct Wrapper {
        static double call(double x, void* p) {
            return (*static_cast<Func*>(p))(x);
        }
    };

    gsl_function F;
    F.function = &Wrapper::call;
    F.params   = f;

    double result = 0.0, error = 0.0;
    gsl_set_error_handler_off();
    int status = gsl_integration_qawc(&F, lower, upper, c,
                                      1e-9, 1e-9, 100000, workspace, &result, &error);
    return result;
}






// RePi_Diquark function
double RePi_Pion(double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, double Nf, double Lambda_CutOff) {
    double sum = 0.0;
    double upper_cut_2 = 2e8 + q*q;
    double upper_cut_1 = 2e4;
    double Lambda_4 = 4 * (Lambda_CutOff * Lambda_CutOff + m * m) + q*q;

    #pragma omp parallel num_threads(5)
    {
        gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100000);
        double local_sum = 0.0;

        if (thermal_cutoff == 1){
            
            auto func1 = [&](double x) {
                return ImPi_Pion(x, q, T, 1, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum += integrate_cauchy_lambda(&func1, -sqrt(Lambda_4), sqrt(Lambda_4), w, workspace);


        } else {
                        
            auto func1 = [&](double x) {
                return ImPi_Pion(x, q, 0, 0, m0, mu_star0, delta0, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum += integrate_cauchy_lambda(&func1, -sqrt(Lambda_4), sqrt(Lambda_4), w, workspace);


                        // Lambda #2 (i == 1)
            auto func2 = [&](double x) {
                return ImPi_Pion(x, q, 0.0, 0.0, m0, mu_star0, delta0, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum -= integrate_cauchy_lambda(&func2, -upper_cut_1, upper_cut_1, w, workspace);

            // Lambda #3 (i == 2)
            auto func3 = [&](double x) {
                return ImPi_Pion(x, q, T, 1.0, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
            };

            local_sum += integrate_cauchy_lambda(&func3, -upper_cut_1, upper_cut_1, w, workspace);

        }

            #pragma omp atomic
            
            sum += local_sum;

            gsl_integration_workspace_free(workspace);
    }

    return (1.0 / pi) * sum;
}






// RePi_Diquark function
double RePi_Diquark(double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, double Nf, double Lambda_CutOff) {
    double sum = 0.0;
    double upper_cut_2 = 2e8 + q*q;
    double upper_cut_1 = 2e4;
    double Lambda_4 = 4 * (Lambda_CutOff * Lambda_CutOff + m * m) + q*q;

    #pragma omp parallel num_threads(5)
    {
        gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100000);
        double local_sum = 0.0;

        if (thermal_cutoff == 1){
            
            auto func1 = [&](double x) {
                return ImPi_Diquark(x, q, T, 1, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum += integrate_cauchy_lambda(&func1, -sqrt(Lambda_4) - 2*mu_star, sqrt(Lambda_4) - 2*mu_star, w, workspace);


        } else {
                        
            auto func1 = [&](double x) {
                return ImPi_Diquark(x, q, 0, 0, m0, mu_star0, delta0, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum += integrate_cauchy_lambda(&func1, -sqrt(Lambda_4), sqrt(Lambda_4), w + 2*mu_star, workspace);


                        // Lambda #2 (i == 1)
            auto func2 = [&](double x) {
                return ImPi_Diquark(x, q, 0.0, 0.0, m0, mu_star0, delta0, thermal_cutoff, Nf, Lambda_CutOff);
            };
            local_sum -= integrate_cauchy_lambda(&func2, -upper_cut_1, upper_cut_1, w + 2*mu_star, workspace);

            // Lambda #3 (i == 2)
            auto func3 = [&](double x) {
                return ImPi_Diquark(x - 2*mu_star, q, T, 1.0, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
            };

            local_sum += integrate_cauchy_lambda(&func3, -upper_cut_1, upper_cut_1, w + 2*mu_star, workspace);

        }

            #pragma omp atomic
            
            sum += local_sum;

            gsl_integration_workspace_free(workspace);
    }

    return (1.0 / pi) * sum;
}







// RePi_Diquark function
double RePi_AntiDiquark(double w, double q, double T, double mu, double m, double mu_star, double delta, double m0, double mu_star0, double delta0, double thermal_cutoff, double Nf, double Lambda_CutOff) {
    double sum = 0.0;
    double upper_cut_2 = 2e8;
    double upper_cut_1 = 2e4;
    double Lambda_4 = 4 * (Lambda_CutOff * Lambda_CutOff + m0 * m0) + q*q;


    #pragma omp parallel num_threads(8)
    {
        gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100000);
        double local_sum = 0.0;

        // Lambda #2 (i == 1)
        auto func1 = [&](double x) {
            double safe_x = std::max(x, 1e-12);
            return ImPi_Diquark(sqrt(safe_x), q, 0.0, 0.0, m0, mu_star0, delta0, thermal_cutoff, Nf, Lambda_CutOff);
        };
        local_sum -= integrate_cauchy_lambda(&func1, Lambda_4, upper_cut_2, w * w, workspace);

        // Lambda #3 (i == 2)
        auto func2 = [&](double x) {
            return ImPi_Diquark(x, q, T, 1.0, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
        };

        local_sum += integrate_cauchy_lambda(&func2, 2 * mu_star, upper_cut_1, w, workspace);

        // Lambda #4 (i == 3)
        auto func3 = [&](double x) {
            return ImPi_Diquark(x, q, T, -1.0, m, mu_star, delta, thermal_cutoff, Nf, Lambda_CutOff);
        };
        local_sum += integrate_cauchy_lambda(&func3, -2 * mu_star, upper_cut_1, -w, workspace);

        #pragma omp atomic
        sum += local_sum;

        gsl_integration_workspace_free(workspace);
    }

    return (1.0 / M_PI) * sum;
}





















double ImPi_AntiDiquark(double w, double q, double T, double p, double m, double mu_star, double Delta, double thermal_cutoff, double Nf, double Lambda_CutOff) {

    const double mu_star_eff = p * mu_star;  // use local const

    double wp = w - 2.0 * mu_star_eff;
    double sp = wp * wp - q * q;

    double P = Nf / (4.0 * pi);

    // Set thermal cutoff based on thermal_cutoff value (0 for inf, 1 for 651)
    double L = (thermal_cutoff == 0) ? std::numeric_limits<double>::infinity() : 651.0;
    
    // Modify this value based on the value of the thermal cutoff
    const double lower_cutoff = sp - 4.0 * m * m;
    const double upper_cutoff = 4.0 * (Lambda_CutOff * Lambda_CutOff + m * m) - sp;
    const double upper_thermal_cutoff = 4.0 * (L*L + m * m) - sp;

    double Jp = 0.0;

    if (lower_cutoff >= 0) {
        if (q == 0.0) {
            // Avoid repeated sqrt calculation
            double root_term = std::sqrt(std::max(0.0, 1.0 - 4.0 * m * m / sp));
            Jp = root_term * (1 - 2.0 * n(w / 2.0, T)) * heaviside(upper_thermal_cutoff);
        } else {
            double ep = wp / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
            double em = wp / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
            if (T != 0) {
                Jp = T * std::log((n(-em - mu_star, T) * n(em + mu_star, T)) / (n(-ep - mu_star, T) * n(ep + mu_star, T))) * heaviside(upper_thermal_cutoff);
            //Jp = T * std::log( (1 + cosh((w+q)/(2*T)))/(1 + cosh((w-q)/(2*T))) )*heaviside(upper_cutoff);
            } else {
                Jp = (std::abs(ep - mu_star) - std::abs(em - mu_star)) * heaviside(upper_thermal_cutoff);
            }
        }
    } 
    


    

    
    if (q != 0) {
        P *= 1 / q;

        if (sp < 0) {
        double ep = wp / 2 + (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
        double em = wp / 2 - (q / 2) * std::sqrt(1 - (4 * m * m) / sp);
        double a = n(em + mu_star, T);
        double b = n(-ep - mu_star, T);

        if (abs(a-b) < 1e-7) {
            Jp += 0;
        } else{
            Jp += ( T * log( a/b ) );
        }


        // double a = n(em + mu_star, T);
        // double b = n(em - mu_star, T);
        // double c = n(-ep + mu_star, T);
        // double d = n(-ep - mu_star, T);

        // if (abs(a*b - c*d) < 1e-10) {
        //     Jp += 0;
        // } else{
        //     Jp += ( T * log( a*b/(c*d) ) );
        // }
    } else{
        Jp += 0;
    }
    
    }

    return P * (sp * Jp);
}



// // Integrate with Cauchy principal value around c
// template <typename Func>
// double integrate_cauchy_lambda(Func* f, double lower, double upper, double c,
//                                gsl_integration_workspace* workspace)
// {
//     struct Wrapper {
//         static double call(double x, void* p) {
//             return (*static_cast<Func*>(p))(x);
//         }
//     };

//     gsl_function F;
//     F.function = &Wrapper::call;
//     F.params   = f;

//     double result = 0.0, error = 0.0;
//     gsl_set_error_handler_off();
//     int status = gsl_integration_qawc(&F, lower, upper, c,
//                                       1e-9, 1e-9, 100000, workspace, &result, &error);
//     // Consider logging error on failure
//     return result;
// }




// double RePi_Diquark(double w, double q, double T, double mu, double m,
//                     double mu_star, double delta, double m0,
//                     double mu_star0, double delta0,
//                     const std::string& ND_Integral_Bound,
//                     double Nf, double Lambda_CutOff) 
// {
//     double sum = 0.0;

//     constexpr double upper_cut_2 = 2e8;
//     constexpr double upper_cut_1 = 2e4;
//     const double Lambda_4 = 4.0 * (Lambda_CutOff * Lambda_CutOff + m0 * m0);

//     const size_t workspace_size = 5000;
//     const int num_threads = 8;

//     #pragma omp parallel num_threads(num_threads) reduction(+:sum)
//     {
//         gsl_integration_workspace* ws0 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws1 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws2 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws3 = gsl_integration_workspace_alloc(workspace_size);

//         double local_sum = 0.0;

//         // Copy parameters into local variables for safe capture
//         const double local_w = w;
//         const double local_q = q;
//         const double local_T = T;
//         const double local_m0 = m0;
//         const double local_mu_star0 = mu_star0;
//         const double local_delta0 = delta0;
//         const double local_m = m;
//         const double local_mu_star = mu_star;
//         const double local_delta = delta;
//         const double local_Nf = Nf;
//         const double local_Lambda_CutOff = Lambda_CutOff;
//         const std::string local_ND_Integral_Bound = ND_Integral_Bound;

//         // Lambda #1
//         auto func0 = [=](double x) {
//             double safe_x = std::max(x, 1e-12);
//             return ImPi_Diquark(std::sqrt(safe_x), local_q, 0.0, 0.0, local_m0, local_mu_star0, local_delta0,
//                                 local_ND_Integral_Bound, local_Nf, local_Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(&func0, 1e-6, Lambda_4, local_w * local_w, ws0);

//         // Lambda #2
//         auto func1 = [=](double x) {
//             double safe_x = std::max(x, 1e-12);
//             return ImPi_Diquark(std::sqrt(safe_x), local_q, 0.0, 0.0, local_m0, local_mu_star0, local_delta0,
//                                 local_ND_Integral_Bound, local_Nf, local_Lambda_CutOff);
//         };
//         local_sum -= integrate_cauchy_lambda(&func1, 1e-6, upper_cut_2, local_w * local_w, ws1);

//         // Lambda #3
//         auto func2 = [=](double x) {
//             return ImPi_Diquark(x, local_q, local_T, 1.0, local_m, local_mu_star, local_delta,
//                                 local_ND_Integral_Bound, local_Nf, local_Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(&func2, -2.0 * local_mu_star, upper_cut_1, local_w, ws2);

//         // Lambda #4
//         auto func3 = [=](double x) {
//             return ImPi_Diquark(x, local_q, local_T, -1.0, local_m, local_mu_star, local_delta,
//                                 local_ND_Integral_Bound, local_Nf, local_Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(&func3, 2.0 * local_mu_star, upper_cut_1, -local_w, ws3);

//         gsl_integration_workspace_free(ws0);
//         gsl_integration_workspace_free(ws1);
//         gsl_integration_workspace_free(ws2);
//         gsl_integration_workspace_free(ws3);

//         sum += local_sum;
//     }

//     return (1.0 / pi) * sum;
// }



// struct RegularFunctionParams {
//     double q;
//     double T;
//     double p;
//     double m;
//     double mu_star;
//     double delta;
//     double thermal_cutoff;
//     double Nf;
//     double Lambda_CutOff;
// };





// template <typename Func>
// double integrate_cauchy_lambda(Func* f, double lower, double upper, double c,
//                                gsl_integration_workspace* workspace)
// {
//     struct Wrapper {
//         static double call(double x, void* p) {
//             return (*static_cast<Func*>(p))(x);
//         }
//     };

//     gsl_function F;
//     F.function = &Wrapper::call;
//     F.params   = f;

//     double result = 0.0, error = 0.0;
//     gsl_set_error_handler_off();
//     int status = gsl_integration_qawc(&F, lower, upper, c,
//                                       1e-9, 1e-9, 100000, workspace, &result, &error);
//     // if (status != GSL_SUCCESS) {
//     //     std::cerr << "GSL qawc error: " << gsl_strerror(status) << std::endl;
//     //     return 0.0;
//     // }
//     return result;
// }

// double RePi_Diquark(double w, double q, double T, double mu, double m,
//                     double mu_star, double delta, double m0,
//                     double mu_star0, double delta0,
//                     double thermal_cutoff,
//                     double Nf, double Lambda_CutOff) 
// {
//     double sum = 0.0;
//     double upper_cut_2 = 2e8;
//     double upper_cut_1 = 2e4;
//     double Lambda_4 = 4 * (Lambda_CutOff * Lambda_CutOff + m0 * m0);

//     #pragma omp parallel num_threads(8)
//     {
//         gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100000);
//         double local_sum = 0.0;


//         // Lambda #2 (i == 1)
//         auto func1 = [&](double x) {
//             double safe_x = std::max(x, 1e-12);
//             return ImPi_Diquark( sqrt(safe_x), q, 0.0, 0.0, m0, mu_star0, delta0,
//                                 thermal_cutoff, Nf, Lambda_CutOff);
//         };
//         local_sum -= integrate_cauchy_lambda(&func1, Lambda_4, upper_cut_2, w * w, workspace);



//         // Lambda #3 (i == 2)
//         auto func2 = [&](double x) {
//             return ImPi_Diquark(x, q, T, 1.0, m, mu_star, delta,
//                                 thermal_cutoff, Nf, Lambda_CutOff);
//         };


//         local_sum += integrate_cauchy_lambda(&func2, -2 * mu_star, upper_cut_1, w, workspace);



//         // Lambda #4 (i == 3)
//         auto func3 = [&](double x) {
//             return ImPi_Diquark(x, q, T, -1.0, m, mu_star, delta,
//                                 thermal_cutoff, Nf, Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(&func3, 2 * mu_star, upper_cut_1, -w, workspace);

//         #pragma omp atomic
//         sum += local_sum;

//         gsl_integration_workspace_free(workspace);
//     }

//     return (1.0 / pi) * sum;
// }










// double integrate_cauchy_lambda(std::function<double(double)> f,
//                                double lower, double upper, double c,
//                                gsl_integration_workspace* workspace)
// {
//     struct Wrapper {
//         static double call(double x, void* params) {
//             auto& func = *static_cast<std::function<double(double)>*>(params);
//             return func(x);
//         }
//     };

//     gsl_function F;
//     F.function = &Wrapper::call;
//     F.params   = &f;

//     double result = 0.0, error = 0.0;
//     gsl_set_error_handler_off();
//     gsl_integration_qawc(&F, lower, upper, c, 1e-9, 1e-9, workspace->limit, workspace, &result, &error);
//     return result;
// }

// // RePi_Diquark function
// double RePi_Diquark(double w, double q, double T, double mu, double m,
//                     double mu_star, double delta, double m0,
//                     double mu_star0, double delta0,
//                     const std::string& ND_Integral_Bound,
//                     double Nf, double Lambda_CutOff)
// {
//     double sum = 0.0;
//     constexpr double upper_cut_2 = 2e8;
//     constexpr double upper_cut_1 = 2e4;
//     const double Lambda_4 = 4.0 * (Lambda_CutOff * Lambda_CutOff + m0 * m0);
//     const size_t workspace_size = 20000;

//     #pragma omp parallel num_threads(10) reduction(+:sum)
//     {
//         // Allocate four separate workspaces per thread
//         gsl_integration_workspace* ws0 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws1 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws2 = gsl_integration_workspace_alloc(workspace_size);
//         gsl_integration_workspace* ws3 = gsl_integration_workspace_alloc(workspace_size);

//         double local_sum = 0.0;

//         // Lambda #1
//         auto func0 = [&](double x) -> double {
//             double safe_x = std::max(x, 1e-12);
//             return ImPi_Diquark(std::sqrt(safe_x), q, 0.0, 0.0, m0, mu_star0, delta0,
//                                 ND_Integral_Bound, Nf, Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(std::function<double(double)>(func0),
//                                              1e-6, Lambda_4, w*w, ws0);

//         // Lambda #2
//         auto func1 = [&](double x) -> double {
//             double safe_x = std::max(x, 1e-12);
//             return ImPi_Diquark(std::sqrt(safe_x), q, 0.0, 0.0, m0, mu_star0, delta0,
//                                 ND_Integral_Bound, Nf, Lambda_CutOff);
//         };
//         local_sum -= integrate_cauchy_lambda(std::function<double(double)>(func1),
//                                              1e-6, upper_cut_2, w*w, ws1);

//         // Lambda #3
//         auto func2 = [&](double x) -> double {
//             return ImPi_Diquark(x, q, T, 1.0, m, mu_star, delta,
//                                 ND_Integral_Bound, Nf, Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(std::function<double(double)>(func2),
//                                              -2.0*mu_star, upper_cut_1, w, ws2);

//         // Lambda #4
//         auto func3 = [&](double x) -> double {
//             return ImPi_Diquark(x, q, T, -1.0, m, mu_star, delta,
//                                 ND_Integral_Bound, Nf, Lambda_CutOff);
//         };
//         local_sum += integrate_cauchy_lambda(std::function<double(double)>(func3),
//                                              2.0*mu_star, upper_cut_1, -w, ws3);

//         // Free workspaces
//         gsl_integration_workspace_free(ws0);
//         gsl_integration_workspace_free(ws1);
//         gsl_integration_workspace_free(ws2);
//         gsl_integration_workspace_free(ws3);

//         sum += local_sum;
//     }

//     return (1.0 / pi) * sum;
// }
