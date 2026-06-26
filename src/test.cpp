#include "RKF78.hpp"
#include <iostream>
#include <vector>
#include <cmath>

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

  auto results = integrate(t0, tf, y0, f2, params, options);
  std::cout << "Accepted step samples: " << results.t.size() << std::endl;
  if (!results.y.empty())
  {
    std::cout << "First step value: t=" << results.t.front()
          << " y=(" << results.y.front()[0] << ", " << results.y.front()[1] << ")" << std::endl;
    std::cout << "Last step value: t=" << results.t.back()
          << " y=(" << results.y.back()[0] << ", " << results.y.back()[1] << ")" << std::endl;
  }
    std::cout << "Final value: " << results.yf[0] << ", " << results.yf[1] << std::endl;
    std::cout << "Number of steps: " << results.accepted << std::endl;
    std::cout << "Number of rejected steps: " << results.rejected << std::endl;
    std::cout << "Number of function evaluations: " << results.fevals << std::endl;



    // // check the accepted step values against the exact solution
    // for (std::size_t i = 0; i < results.t.size(); ++i)
    // {
    //   double t_step = results.t[i];
    //   std::vector<double> y_step = results.y[i];
    //   std::vector<double> y_exact = {exp(cos(t_step*t_step)),exp(sin(t_step*t_step))};

    //   std::cout << "t=" << t_step << " y_step=(" << y_step[0] << ", " << y_step[1] << ")"
    //             << " y_exact=(" << y_exact[0] << ", " << y_exact[1] << ")" << std::endl;

    //   double error0 = std::abs(y_step[0] - y_exact[0])/std::abs(y_exact[0]);
    //   double error1 = std::abs(y_step[1] - y_exact[1])/std::abs(y_exact[1]);

    //   std::cout << "Error in y[0]: " << error0 << std::endl;
    //   std::cout << "Error in y[1]: " << error1 << std::endl;


    
    // }

  


  return 0;
}
