#include "RKF78.hpp"
#include <iostream>

using namespace RKF78;

ODEFunction<std::array<double,1>, std::array<double,1>> f = [](const double t, const std::array<double,1> &y, std::array<double,1> &dy, const std::array<double,1> &params)
{
  dy[0] = -params[0] * y[0];
};

ODEFunction<std::array<double,2>, std::array<double,0>> f2 = [](const double t, const std::array<double,2> &y, std::array<double,2> &dy, const std::array<double,0> &params)
{
  dy[0] = -2*t*y[0]*log(y[1]);
  dy[1] = 2*t*y[1]*log(y[0]);
 
};

int main(int argc, char *argv[])
{
    std::array<double,2> y0 = {exp(1.0),1.0};
    std::array<double,0> params = {};
    double t0 = 0.0;
    double tf = 5.0;

    Options options;
    options.h0 =0.1;
    options.atol = 1e-16;
    options.rtol = 1e-16;

    Results<std::array<double,2>> results = integrate(t0, tf, y0, f2, params, options);
    std::cout << "Final value: " << results.yf[0] << ", " << results.yf[1] << std::endl;
    std::cout << "Number of steps: " << results.accepted << std::endl;
    std::cout << "Number of rejected steps: " << results.rejected << std::endl;
    std::cout << "Number of function evaluations: " << results.fevals << std::endl;

    //Exact solution at t=tf
    double exact_y0 = exp(cos(pow(tf,2)));
    double exact_y1 = exp(sin(pow(tf,2)));

    // Compute the absolute error
    double error_y0 = std::abs(results.yf[0] - exact_y0);
    double error_y1 = std::abs(results.yf[1] - exact_y1);

    //Compute the relative error
    double rel_error_y0 = error_y0 / std::abs(exact_y0);
    double rel_error_y1 = error_y1 / std::abs(exact_y1);

    // Print the results
    std::cout << "Exact solution: " << exact_y0 << ", " << exact_y1 << std::endl;
    std::cout << "Absolute error: " << error_y0 << ", " << error_y1 << std::endl;
    std::cout << "Relative error: " << rel_error_y0 << ", " << rel_error_y1 << std::endl;
  


  return 0;
}
