#include "RKF78.hpp"
#include <iostream>
#include <vector>
#include <cmath>

using namespace RKF78;

//params stored in a std::array
void f1(const double t, const std::array<double,2> &y, std::array<double,2> &dy, const std::array<double,2> params)
{
  dy[0] = params[0]*t*y[0]*log(y[1]);
  dy[1] = params[1]*t*y[1]*log(y[0]);
};

//pointer to the parameters array
void f2(const double t, const std::array<double,2> &y, std::array<double,2> &dy, const double* params)
{
  dy[0] = params[0]*t*y[0]*log(y[1]);
  dy[1] = params[1]*t*y[1]*log(y[0]);
};

//no params at all
void f3(const double t, const std::array<double,2> &y, std::array<double,2> &dy)
{
  dy[0] = -2.0*t*y[0]*log(y[1]);
  dy[1] = 2.0*t*y[1]*log(y[0]);
};

int main(int argc, char *argv[])
{
    std::array<double,2> y0 = {exp(1.0),1.0};
    double params_ptr[2] = {-2.0,2.0};
    std::array<double,2> params_array = {-2.0,2.0};
    double t0 = 0.0;
    double tf = 5.0;

    Options options;
    options.h0 =0.1;
    options.atol = 1e-16;
    options.rtol = 1e-16;

  auto results = integrate(t0, tf, y0, f1, params_array, options);
  auto results2 = integrate(t0, tf, y0, f2, params_ptr, options);
  auto results3 = integrate(t0, tf, y0, f3, options);

  //ouput the final state for each integration
  std::cout << "Results from f1 with params as std::array: " << std::endl;
  std::cout << "Final value: " << results.yf[0] << ", " << results.yf[1] << std::endl;
  std::cout << "Results from f2 with params as pointer: " << std::endl;
  std::cout << "Final value: " << results2.yf[0] << ", " << results2.yf[1] << std::endl;
  std::cout << "Results from f3 with no params: " << std::endl;
  std::cout << "Final value: " << results3.yf[0] << ", " << results3.yf[1] << std::endl;


    std::cout << "Accepted step samples: " << results.t.size() << std::endl;
    if (results.t.size() > 0)
    {
      const auto &y_first = results.y_vector(0);
      const auto &y_last  = results.y_vector(results.t.size()-1);
      std::cout << "First step value: t=" << results.t.front()
        << " y=(" << y_first[0] << ", " << y_first[1] << ")" << std::endl;
      std::cout << "Last step value: t=" << results.t.back()
        << " y=(" << y_last[0] << ", " << y_last[1] << ")" << std::endl;
    }
    std::cout << "Number of steps: " << results.accepted << std::endl;
    std::cout << "Number of rejected steps: " << results.rejected << std::endl;
    std::cout << "Number of function evaluations: " << results.fevals << std::endl;



    // check the accepted step values against the exact solution
    for (std::size_t i = 0; i < results.t.size(); i+=100)
    {
      double t_step = results.t[i];
      const std::vector<double>& y_step = results.y_vector(i);
      std::vector<double> y_exact = {exp(cos(t_step*t_step)),exp(sin(t_step*t_step))};

      std::cout << "t=" << t_step << " y_step=(" << y_step[0] << ", " << y_step[1] << ")"
                << " y_exact=(" << y_exact[0] << ", " << y_exact[1] << ")" << std::endl;

      double error0 = std::abs(y_step[0] - y_exact[0])/std::abs(y_exact[0]);
      double error1 = std::abs(y_step[1] - y_exact[1])/std::abs(y_exact[1]);

      std::cout << "Error in y[0]: " << error0 << std::endl;
      std::cout << "Error in y[1]: " << error1 << std::endl;


    
    }
  return 0;
}
