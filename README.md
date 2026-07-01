# Runge-Kutta-Fehlberg7(8)
A header-only c++ implementation of the 7/8th order Runge-Kutta-Fehlberg method. In this program the 8th order solution is used to compute the next step. Error is estimated using the difference between the 7th and 8th order solutions.

Implemented as a template class, this library allows for the numerical solution of ordinary differential equations (ODEs) of the form dy/dt = f(t, y), where y is a vector of dependent variables and t is the independent variable.

Input state can be provided as any container that supports container.size() and container.data() methods, such as std::vector, std::array, or Eigen::VectorXd. 

## Quick Start

1) Include the header in your C++ file:

```cpp
#include "RKF78.hpp"
```

2) Define your ODE dy/dt = f(t, y):

```cpp
// Example without params: y' = -0.5 y
void f_np(double, const std::array<double,1>& y, std::array<double,1>& dy) {
	dy[0] = -0.5 * y[0];
}
```

3) Integrate (no options):

```cpp
std::array<double,1> y0{1.0};
double t0 = 0.0, tf = 5.0;

RKF78::Results res = RKF78::integrate(t0, tf, y0, f_np);
```

4) Use results:

```cpp
// Accepted sample times and states
res.t;        // std::vector<double>
res.y;        // std::vector<std::vector<double>>
res.yf;       // final state (std::vector<double>)
res.accepted; // number of accepted steps
res.rejected; // number of rejected steps
res.fevals;   // function evaluations
```

### Build and run the included example

```bash
mkdir -p build && cd build
cmake ..
cmake --build . --config Release
./RKF78_example
```

### Documentation

To generate the API reference (HTML) with Doxygen:

```bash
cd build
cmake ..
cmake --build . --target docs
open ../docs/build/html/index.html  # macOS convenience
```

Online docs (GitHub Pages):

- This repository is configured to publish the Doxygen HTML to GitHub Pages via CI.
- After pushing to your default branch and enabling Pages (Source: GitHub Actions) in the repository settings, the site will be available at a URL like:

	https://YOUR-GITHUB-USERNAME.github.io/Runge-Kutta-Fehlberg78/

	Replace YOUR-GITHUB-USERNAME with your actual username or org name.

Notes:
- Header-only: just add `include/` to your include path and `#include "RKF78.hpp"`.
- State/parameter types can be `std::array`, `std::vector`, or any container with `.size()` and random access '[]'.
- The solver uses the 8th-order solution; error is estimated from the 7th/8th order difference.
 - Call forms:
	 - No params: `RKF78::integrate(t0, tf, y0, f_np)` or with options `..., f_np, opt)`.
	 - With params: `RKF78::integrate(t0, tf, y0, f, params)` or with options `..., f, params, opt)`.

	### References

	1. E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control," NASA Technical Report R-287, 1968. NTRS: https://ntrs.nasa.gov/citations/19680027281

### With params and options

```cpp
// ODE with parameters p = [k1, k2]
void f_with_params(double t, const std::array<double,1>& y,
                   std::array<double,1>& dy, const std::array<double,1>& p)
{
  dy[0] = -p[0] * y[0];
}
std::array<double,1> y0{1.0}; 
std::array<double,1> params{0.5};
RKF78::Options opt; opt.h0=0.1; // Initial step size
opt.atol=1e-14; // Absolute tolerance
opt.rtol=1e-14; // Relative tolerance
RKF78::Results res2 = RKF78::integrate(0.0, 5.0, y0, f_with_params, params, opt); 
```

### Scalar (double) usage

```cpp
// No params, no options (state is a scalar double)
void f_scalar(double /*t*/, const double& y, double& dy) {
	dy = -0.5 * y;
}
double y0s = 1.0;
RKF78::Results rs = RKF78::integrate(0.0, 5.0, y0s, f_scalar);
```

```cpp
// With params (and options)
struct P1 { double k; };
void f_scalar_p(double /*t*/, const double& y, double& dy, const P1& p) {
	dy = -p.k * y;
}
double y0p = 1.0; P1 p{0.4};
RKF78::Options opt2; opt2.h0=0.1; opt2.atol=1e-12; opt2.rtol=1e-12;
RKF78::Results rsp = RKF78::integrate(0.0, 5.0, y0p, f_scalar_p, p, opt2);
```