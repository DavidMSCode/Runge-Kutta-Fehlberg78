# Quick Start: RKF78 ODE Solver

## Quick Start

1. **Include the header:**
   ```cpp
   #include "RKF78.hpp"
   ```

2. **Define your ODE function (no params)** dy/dt = f(t, y):
   ```cpp
     // Free function without params (non-capturing)
     void f_np(double t, const std::array<double,2>& y, std::array<double,2>& dy)
     {
        dy[0] = -2*t*y[0]*log(y[1]);
        dy[1] =  2*t*y[1]*log(y[0]);
     }
   ```

3. **Integrate (no options):**
   ```cpp
   std::array<double,2> y0 = {exp(1.0), 1.0};  // Initial state
   double t0 = 0.0, tf = 5.0;                    // Time span
   
   RKF78::Results results = RKF78::integrate(t0, tf, y0, f_np);
   ```

5. **Access results:**
   ```cpp
   results.t          // Vector of accepted step times
   results.y          // Vector of state vectors at each step
   results.yf         // Final state vector
   results.accepted   // Number of accepted steps
   results.rejected   // Number of rejected steps
   results.fevals     // Total function evaluations
   ```

## Build & Run Example
```bash
cd build
cmake ..
make
./RKF78_example
```

## Notes
- State vectors can be `std::array`, `std::vector`, or any container with `.size()` and `.data()`
 - Overloads:
    - No params: `integrate(t0, tf, y0, f)` or `integrate(t0, tf, y0, f, options)`
    - With params: `integrate(t0, tf, y0, f_with_params, params)` or `... , params, options)`

---

## Alternative ways and options/params example

You can also use these forms as long as the signature matches:

```cpp
// No params:    void f(double t, const State& y, State& dy)
// With params:  void f(double t, const State& y, State& dy, const Params& p)
// With options no params: void f(double t, const State& y, State& dy, const Options& opt)
// With options and params: void f(double t, const State& y, State& dy, const Params& p, const Options& opt)
```

- Function with parameters and options with stricter error control:
```cpp
void f_with_params(double t, const std::array<double,2>& y,
                        std::array<double,2>& dy, const std::array<double,2>& p)
{
  dy[0] = p[0]*t*y[0]*log(y[1]);
  dy[1] =  p[1]*t*y[1]*log(y[0]);
};
std::array<double,2> y0 = {exp(1.0), 1.0};
std::array<double,2> params = {-2.0, 2.0};
RKF78::Options options; 
options.h0=0.1;     // Initial step size
options.atol=1e-14; // Absolute tolerance
options.rtol=1e-14; // Relative tolerance
RKF78::Results res_lp = RKF78::integrate(0.0, 5.0, y0, f_with_params, params, options);
```
