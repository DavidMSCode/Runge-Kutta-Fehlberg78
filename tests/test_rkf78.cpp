#include "RKF78.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>

using RKF78::Options;
using RKF78::Results;

static bool approx(double a, double b, double rtol=1e-9, double atol=1e-12) {
  double diff = std::abs(a-b);
  double scale = atol + rtol*std::max(std::abs(a), std::abs(b));
  return diff <= scale;
}

static int test_scalar_no_params_no_options() {
  // y' = -k y, k=0.5
  auto f = [](double /*t*/, const double& y, double& dy) {
    dy = -0.5 * y;
  };
  double y0 = 1.0;
  double t0 = 0.0, tf = 2.0;
  Results res = RKF78::integrate(t0, tf, y0, f);
  double exact = std::exp(-0.5*(tf - t0));
  if (!approx(res.yf[0], exact)) {
    std::cerr << "FAIL scalar no-params no-options: got=" << res.yf[0]
              << " exact=" << exact << std::endl;
    return 1;
  }
  if (res.accepted <= 0 || res.rejected < 0) return 2;
  return 0;
}

static int test_scalar_no_params_with_options() {
  auto f = [](double /*t*/, const double& y, double& dy) { dy = -1.0 * y; };
  double y0 = 2.0; double t0=0.0, tf=1.0;
  Options opt; opt.h0=0.05; opt.atol=1e-12; opt.rtol=1e-12;
  Results res = RKF78::integrate(t0, tf, y0, f, opt);
  double exact = y0*std::exp(-1.0*(tf-t0));
  if (!approx(res.yf[0], exact)) {
    std::cerr << "FAIL scalar no-params with-options: got=" << res.yf[0]
              << " exact=" << exact << std::endl;
    return 1;
  }
  return 0;
}

struct Params1 { double k; };

static int test_scalar_with_params_no_options() {
  auto f = [](double /*t*/, const double& y, double& dy, const Params1& p) {
    dy = -p.k * y;
  };
  double y0 = 1.5; double t0=0.0, tf=1.2; Params1 P{0.3};
  Results res = RKF78::integrate(t0, tf, y0, f, P);
  double exact = y0*std::exp(-P.k*(tf-t0));
  if (!approx(res.yf[0], exact)) {
    std::cerr << "FAIL scalar with-params no-options: got=" << res.yf[0]
              << " exact=" << exact << std::endl;
    return 1;
  }
  return 0;
}

static int test_scalar_with_params_with_options() {
  auto f = [](double /*t*/, const double& y, double& dy, const Params1& p) {
    dy = -p.k * y;
  };
  double y0 = 1.5; double t0=0.0, tf=1.2; Params1 P{0.7};
  Options opt; opt.h0=0.02; opt.atol=1e-12; opt.rtol=1e-12;
  Results res = RKF78::integrate(t0, tf, y0, f, P, opt);
  double exact = y0*std::exp(-P.k*(tf-t0));
  if (!approx(res.yf[0], exact)) {
    std::cerr << "FAIL scalar with-params with-options: got=" << res.yf[0]
              << " exact=" << exact << std::endl;
    return 1;
  }
  return 0;
}

static int test_array2_no_params_no_options() {
  auto f = [](double /*t*/, const std::array<double,2>& y, std::array<double,2>& dy) {
    dy[0] = -0.5 * y[0];
    dy[1] = -1.0 * y[1];
  };
  std::array<double,2> y0{1.0, 2.0}; double t0=0.0, tf=1.0;
  Results res = RKF78::integrate(t0, tf, y0, f);
  double e0 = y0[0]*std::exp(-0.5*(tf-t0));
  double e1 = y0[1]*std::exp(-1.0*(tf-t0));
  if (!approx(res.yf[0], e0) || !approx(res.yf[1], e1)) {
    std::cerr << "FAIL array2 no-params no-options: got=("<<res.yf[0]<<","<<res.yf[1]
              <<") exact=("<<e0<<","<<e1<<")"<<std::endl;
    return 1;
  }
  return 0;
}

struct Params2 { double k1, k2; };

static int test_array2_with_params_with_options() {
  auto f = [](double /*t*/, const std::array<double,2>& y, std::array<double,2>& dy, const Params2& p) {
    dy[0] = -p.k1 * y[0];
    dy[1] = -p.k2 * y[1];
  };
  std::array<double,2> y0{1.0, 3.0}; double t0=0.0, tf=1.5; Params2 P{0.4, 0.7};
  Options opt; opt.h0=0.03; opt.atol=1e-12; opt.rtol=1e-12;
  Results res = RKF78::integrate(t0, tf, y0, f, P, opt);
  double e0 = y0[0]*std::exp(-P.k1*(tf-t0));
  double e1 = y0[1]*std::exp(-P.k2*(tf-t0));
  if (!approx(res.yf[0], e0) || !approx(res.yf[1], e1)) {
    std::cerr << "FAIL array2 with-params with-options: got=("<<res.yf[0]<<","<<res.yf[1]
              <<") exact=("<<e0<<","<<e1<<")"<<std::endl;
    return 1;
  }
  return 0;
}

int main() {
  int rc = 0;
  rc |= test_scalar_no_params_no_options();
  rc |= test_scalar_no_params_with_options();
  rc |= test_scalar_with_params_no_options();
  rc |= test_scalar_with_params_with_options();
  rc |= test_array2_no_params_no_options();
  rc |= test_array2_with_params_with_options();
  if (rc == 0) {
    std::cout << "All RKF78 tests passed" << std::endl;
  }
  return rc;
}
