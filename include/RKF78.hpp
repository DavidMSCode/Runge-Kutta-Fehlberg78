/*
 * Created on Wed Jun 24 2026
 *
 * Copyright (c) 2026 David Stanley
 */

#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <cstddef>
#include <utility>
#include <type_traits>
#include <variant>
#include <functional>
#include <stdexcept>

#pragma once // Include guard to prevent multiple inclusions of this header file

/** \file RKF78.hpp
 *  \brief Header-only Runge–Kutta–Fehlberg 7(8) ODE solver.
 *
 *  Provides adaptive step-size integration for first-order ODE systems
 *  using the Fehlberg 7/8 tableau. Supports scalar states (double) and
 *  vector-like states (e.g., std::array, std::vector) that expose size()
 *  and operator[]. Error control uses the embedded 7/8 pair.
 */
namespace RKF78
{
  /** \defgroup rkf78_api RKF78 API
   *  Public API for the Runge–Kutta–Fehlberg 7(8) solver.
   *  \brief Primary types and functions exposed to users.
   */
  /** \brief Helper utilities to support both container states and scalar doubles. */
  inline std::size_t rkf_size(const double &)
  {
    return 1;
  }

  template <typename S>
  inline std::size_t rkf_size(const S &s)
  {
    return static_cast<std::size_t>(s.size());
  }

  /** \brief Access element j of a scalar (always returns the scalar itself). */
  inline double &elem(double &s, std::size_t /*i*/)
  {
    return s; // scalar has a single element
  }

  /** \brief Const access to element j of a scalar (always returns the scalar itself). */
  inline const double &elem(const double &s, std::size_t /*i*/)
  {
    return s; // scalar has a single element
  }

  template <typename S>
  /** \brief Access element i of a vector-like state. */
  inline auto elem(S &s, std::size_t i) -> decltype(s[i])
  {
    return s[i];
  }

  template <typename S>
  /** \brief Const access to element i of a vector-like state. */
  inline auto elem(const S &s, std::size_t i) -> decltype(s[i])
  {
    return s[i];
  }

  /** \brief Integration results and statistics.
   *  \ingroup rkf78_api
   */
  class Results
  {
  public:
    /** Accepted step time points */
    std::vector<double> t;
    /** Number of rejected steps during integration */
    int rejected = 0;
    /** Number of accepted steps during integration */
    int accepted = 0;
    /** Number of function evaluations performed */
    int fevals = 0;
    /** Final state vector after integration (flattened) */
    std::vector<double> yf;

    Results() = default;

    /** Mark whether this run used a scalar state. */
    void set_scalar(bool v) { scalar_ = v; }
    /** True if the integration was performed with a scalar state (double). */
    bool scalar() const { return scalar_; }

    /** Number of accepted samples (same as t.size()). */
    std::size_t size() const { return y_.size(); }

    /** Internal: append a sampled state (flattened). */
    void add_state(const std::vector<double>& yv) { y_.push_back(yv); }

    /** Access a sampled state as variant: double for scalar, vector ref for multi-dim. */
    auto y(std::size_t i) const -> std::variant<double, std::reference_wrapper<const std::vector<double>>>
    {
      if (scalar_) return y_.at(i).at(0);
      return std::cref(y_.at(i));
    }

    /** Convenience: access vector state (throws if scalar run).
     *  \throws std::logic_error if the integration was performed with a scalar state.
     */
    const std::vector<double>& y_vector(std::size_t i) const
    {
      if (scalar_) throw std::logic_error("Results::y_vector called for scalar run");
      return y_.at(i);
    }

    /** Convenience: access scalar state (throws if vector run).
     *  \throws std::logic_error if the integration was performed with a vector state.
     */
    double y_scalar(std::size_t i) const
    {
      if (!scalar_) throw std::logic_error("Results::y_scalar called for vector run");
      return y_.at(i).at(0);
    }

    /** Proxy view to support bracket indexing: results[i] */
    class SampleView {
    public:
      SampleView(const std::vector<std::vector<double>>* yptr, std::size_t idx, bool scalar)
        : yptr_(yptr), idx_(idx), scalar_(scalar) {}
      /** Implicit conversion to double for scalar runs.
       *  \throws std::logic_error if called on a vector run.
       */
      operator double() const {
        if (!scalar_) throw std::logic_error("Results::operator[] to double on vector run");
        return (*yptr_)[idx_][0];
      }
      /** Implicit conversion to const std::vector<double>& for vector runs.
       *  \throws std::logic_error if called on a scalar run.
       */
      operator const std::vector<double>&() const {
        if (scalar_) throw std::logic_error("Results::operator[] to vector on scalar run");
        return (*yptr_)[idx_];
      }
      /** Index into the state (works for both scalar (expect j==0) and vector states).
       *  \throws std::out_of_range if j != 0 for scalar runs, or j >= size() for vector runs.
       */
      const double& operator[](std::size_t j) const {
        if (scalar_) {
          if (j != 0) throw std::out_of_range("Results::SampleView scalar index must be 0");
          return (*yptr_)[idx_][0];
        }
        return (*yptr_)[idx_][j];
      }
      std::size_t size() const { return scalar_ ? 1 : (*yptr_)[idx_].size(); }
      /** Access the underlying vector state for vector runs.
       *  \throws std::logic_error if called on a scalar run.
       */
      const std::vector<double>& as_vector() const {
        if (scalar_) throw std::logic_error("Results::SampleView::as_vector on scalar run");
        return (*yptr_)[idx_];
      }
    private:
      const std::vector<std::vector<double>>* yptr_;
      std::size_t idx_;
      bool scalar_;
    };

    /** Bracket access: results[i]
     *  - For scalar runs: convertible to double (e.g., double v = results[i])
     *  - For vector runs: use results[i][j] or results[i].as_vector()
     *  \note Conversions and indexing may throw as documented in SampleView methods.
     */
    SampleView operator[](std::size_t i) const { return SampleView(&y_, i, scalar_); }

  private:
    /** Accepted step state values (each state flattened to std::vector<double>) */
    std::vector<std::vector<double>> y_;
    /** True if the integration was performed with a scalar state (double). */
    bool scalar_ = false;
  };

  /** \brief Solver options (with defaults).
   *  \ingroup rkf78_api
   */
  struct Options
  {
    double h0=0.1;     ///< Initial step size
    double atol=1e-12; ///< Absolute tolerance
    double rtol=1e-12; ///< Relative tolerance
  };

  // Default options used by overloads that omit Options
  inline Options default_options()
  {
    return Options{0.1, 1e-12, 1e-12};
  }

  /** \brief Number of RKF78 stages. */
  constexpr int stages = 13;

  /** \brief Stage time fractions a_k for the Fehlberg 7/8 tableau.
   *
   *  Coefficients follow the classical Fehlberg 7(8) method from:
   *  E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control",
   *  NASA Technical Report R-287, 1968. NTRS: https://ntrs.nasa.gov/citations/19680027281
   */
  static constexpr std::array<double, stages> ak = {
      0.0,
      2.0 / 27.0,
      1.0 / 9.0,
      1.0 / 6.0,
      5.0 / 12.0,
      1.0 / 2.0,
      5.0 / 6.0,
      1.0 / 6.0,
      2.0 / 3.0,
      1.0 / 3.0,
      1.0,
      0.0,
      1.0};

  /** \brief Coefficients beta_{k,lambda} of the lower-triangular RK matrix.
   *
   *  Source: E. Fehlberg, NASA-TR-R-287 (1968); see NTRS citation above.
   */
  static constexpr std::array<std::array<double, stages>, stages> bk =
      {{{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{2.0 / 27.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{1.0 / 36.0,
          1.0 / 12.0,
          0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{1.0 / 24.0,
          0.0,
          1.0 / 8.0,
          0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{5.0 / 12.0,
          0.0,
          -25.0 / 16.0,
          25.0 / 16.0,
          0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{1.0 / 20.0,
          0.0,
          0.0,
          1.0 / 4.0,
          1.0 / 5.0,
          0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{-25.0 / 108.0,
          0.0,
          0.0,
          125.0 / 108.0,
          -65.0 / 27.0,
          125.0 / 54.0,
          0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{31.0 / 300.0,
          0.0,
          0.0,
          0.0,
          61.0 / 225.0,
          -2.0 / 9.0,
          13.0 / 900.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

        {{2.0,
          0.0,
          0.0,
          -53.0 / 6.0,
          704.0 / 45.0,
          -107.0 / 9.0,
          67.0 / 90.0,
          3.0,
          0.0, 0.0, 0.0, 0.0, 0.0}},

        {{-91.0 / 108.0,
          0.0,
          0.0,
          23.0 / 108.0,
          -976.0 / 135.0,
          311.0 / 54.0,
          -19.0 / 60.0,
          17.0 / 6.0,
          -1.0 / 12.0,
          0.0, 0.0, 0.0, 0.0}},

        {{2383.0 / 4100.0,
          0.0,
          0.0,
          -341.0 / 164.0,
          4496.0 / 1025.0,
          -301.0 / 82.0,
          2133.0 / 4100.0,
          45.0 / 82.0,
          45.0 / 164.0,
          18.0 / 41.0,
          0.0, 0.0, 0.0}},

        {{3.0 / 205.0,
          0.0,
          0.0,
          0.0,
          0.0,
          -6.0 / 41.0,
          -3.0 / 205.0,
          -3.0 / 41.0,
          3.0 / 41.0,
          6.0 / 41.0,
          0.0,
          0.0, 0.0}},

        {{-1777.0 / 4100.0,
          0.0,
          0.0,
          -341.0 / 164.0,
          4496.0 / 1025.0,
          -289.0 / 82.0,
          2193.0 / 4100.0,
          51.0 / 82.0,
          33.0 / 164.0,
          12.0 / 41.0,
          0.0,
          1.0,
          0.0}}}};

  /** \brief Convert a vector-like state to std::vector<double>. */
  template <typename State>
  std::vector<double> state_to_vector(const State &state)
  {
    return std::vector<double>(state.begin(), state.end());
  }

  /** \overload */
  inline std::vector<double> state_to_vector(const double &state)
  {
    return std::vector<double>{state};
  }

  /** \brief Adaptive RKF78 integration (with parameters).
   *  \ingroup rkf78_api
   *  \tparam State  State container type (double or vector-like with size() and operator[])
   *  \tparam Ode    Callable f(t, y, dy, params)
   *  \tparam Params Parameter type passed to f
   *  \param t0     Initial time
   *  \param tf     Final time
   *  \param y0     Initial state
   *  \param f      ODE right-hand side callable
   *  \param params User parameters forwarded to f
   *  \param options Solver tolerances and initial step
   *  \return Results containing trajectory samples and statistics
   *
   *  \par Requirements
   *  - State is either double or a vector-like type with .size() and operator[] returning references.
   *  - Ode must be callable as `void f(double t, const State& y, State& dy, const Params& p)`.
   *
   *  \par Example (with parameters)
   *  \code{.cpp}
   *  struct P { double k; } p{0.5};
   *  auto f = [](double /*t*/, const std::array<double,1>& y,
   *              std::array<double,1>& dy, const P& p){ dy[0] = -p.k * y[0]; };
   *  std::array<double,1> y0{1.0};
   *  auto res = RKF78::integrate(0.0, 5.0, y0, f, p, RKF78::default_options());
   *  \endcode
   */
  template <typename State, typename Ode, typename Params>
  Results integrate(const double t0, const double tf, State y0, Ode &&f, Params &&params, Options options)
  {


    double &atol = options.atol;
    double &rtol = options.rtol;
    double &h = options.h0; // initial step size

    // Initialize the deriv variables for the Runge-Kutta-Fehlberg 7(8) method
    std::array<State, stages> k; // Array to hold the derivatives at each stage

    // Initialize the solution storage arrays (need to be resizeable)
    int dims = static_cast<int>(rkf_size(y0));              // Get the size of the state (1 for scalar)
    bool is_scalar = (dims == 1);
    int prop_dir = (tf > t0) ? 1 : -1; // Determine the propagation direction based on the time interval

    // construct the results struct to return
    Results results;
    results.t.reserve(64);
    // Reserve internal storage via number of samples (states)
    results.t.push_back(t0);
    results.set_scalar(is_scalar);
    results.add_state(state_to_vector(y0));

    int &fevals = results.fevals;     // Initialize the function evaluation counter
    int &istep = results.accepted;    // Initialize the step counter
    int &rejected = results.rejected; // Initialize the rejected step counter
    double t = t0;                    // Initialize the current time to the initial time
    int istat = 0;                    // Start integration
    // istat=1 // means integration is complete
    // istat=-1 // means error occurred
    // istat=0 // means integration is in progress
    while (istat == 0)
    {
      bool bad_stage = false; // Flag to indicate if a bad stage has occurred (e.g., due to NaN or Inf values)
      // compute the derivatives at each stage using the RKF78 method
      for (int i = 0; i < stages; ++i)
      {
        State y_temp = y0;
        for (int j = 0; j < i; ++j)
        {
          for (int d = 0; d < dims; ++d)
          {
            elem(y_temp, d) += h * bk[i][j] * elem(k[j], d);
          }
        }
        f(t + ak[i] * h, y_temp, k[i], params);
        fevals++; // Increment the function evaluation counter for each stage
        // check for NaN or Inf values in the computed derivatives
        for (int d = 0; d < dims; ++d)
        {
          if (!std::isfinite(elem(k[i], d)))
          {
            bad_stage = true;
            break;
          }
        }

        if (bad_stage)
        {
          break;
        }

      } // end computing stage derivatives
      if (bad_stage)
      {
        // If a bad stage occurred, reduce the step size and retry
        h *= 0.5;   // Reduce the step size by half
        rejected++; // Increment the rejected step counter
        continue;   // Retry the current step with the new step size
      }

      // Compute the 8th order solution for the step
      State y8 = y0;
      for (int d = 0; d < dims; ++d)
      {
        elem(y8, d) += h * (34.0 / 105.0 * elem(k[5], d) + 9.0 / 35.0 * elem(k[6], d) + 9.0 / 35.0 * elem(k[7], d) + 9.0 / 280.0 * elem(k[8], d) + 9.0 / 280.0 * elem(k[9], d) + 41.0 / 840.0 * elem(k[11], d) + 41.0 / 840.0 * elem(k[12], d));
      } // end computing 8th order solution

      // Compute the truncation error estimate for the step
      State error_estimate;
      for (int d = 0; d < dims; ++d)
      {
        elem(error_estimate, d) = h * 41.0 / 840.0 * (elem(k[0], d) + elem(k[10], d) - elem(k[11], d) - elem(k[12], d));
      }

      double err2 = 0.0; // error scaled by the tolerance

      for (int d = 0; d < dims; ++d)
      {
        double scale = atol + rtol * std::max(std::abs(elem(y0, d)), std::abs(elem(y8, d)));
        double ratio = elem(error_estimate, d) / scale;
        err2 += ratio * ratio;
      }

      double error_norm = std::sqrt(err2 / dims);
      double h_new; // New step size for the next iteration
      if (error_norm == 0.0)
      {
        h_new = 2.0 * h; // If the error is zero, increase the step size
      }
      else
      {
        h_new = 0.9 * h * std::pow(1.0 / error_norm, 1.0 / 8.0);
      }
      if (error_norm <= 1.0) // Check if the error is within the tolerance
      {
        // Accept the step
        y0 = y8; // Update the state vector to the new value
        t += h;  // Update the time to the new value
        h = h_new;
        istep++; // Increment the step counter
        results.t.push_back(t);
        results.add_state(state_to_vector(y0));

        if ((prop_dir > 0 && t >= tf) || (prop_dir < 0 && t <= tf)) // Check if the final time has been reached
        {
          istat = 1; // Integration complete
        }
      }
      else
      {
        // Reject the step and reduce the step size
        h = h_new;  // Update the step size for the next iteration
        rejected++; // Increment the rejected step counter
      }

      // check if h will overshoot the final time
      if ((prop_dir > 0 && t + h > tf) || (prop_dir < 0 && t + h < tf))
      {
        h = tf - t; // Adjust the step size to reach the final time exactly
      }

    } // end while loop

    // Store the final state vector in the results struct
    results.yf = state_to_vector(y0); // Set the final state vector after integration

    return results; // Return the results of the integration
  }

  /** \brief Adaptive RKF78 integration for ODEs without parameters.
   *  \ingroup rkf78_api
   *  Accepts any callable f(t, y, dy) and adapts it to the main integrator.
   *  \tparam State  State container type (double or vector-like)
   *  \tparam Ode    Callable f(t, y, dy)
   *  \param t0     Initial time
   *  \param tf     Final time
   *  \param y0     Initial state
   *  \param f      ODE right-hand side callable
   *  \param options Solver tolerances and initial step
   *  \return Results containing trajectory samples and statistics
   *
   *  \par Requirements
   *  - State is either double or a vector-like type with .size() and operator[] returning references.
   *  - Ode must be callable as `void f(double t, const State& y, State& dy)`.
   *
   *  \par Example (scalar)
   *  \code{.cpp}
   *  auto f = [](double /*t*/, const double& y, double& dy){ dy = -0.5 * y; };
   *  double y0 = 1.0;
   *  auto res = RKF78::integrate(0.0, 5.0, y0, f, RKF78::default_options());
   *  double y_end = res.yf[0];
   *  \endcode
   */
  template <typename State, typename Ode>
  Results integrate(const double t0, const double tf, State y0, Ode &&f, Options options)
  {
    auto adapter = [&](const double t, const State &y, State &dy, const std::nullptr_t &)
    {
      f(t, y, dy);
    };
    return integrate<State, decltype(adapter), std::nullptr_t>(
      t0, tf, std::move(y0), std::move(adapter), nullptr, options);
  }

  /** \brief Convenience overload: with parameters, using default Options.
   *  \ingroup rkf78_api
   */
  template <typename State, typename Ode, typename Params>
  Results integrate(const double t0, const double tf, State y0, Ode &&f, Params &&params)
  {
    return integrate(t0, tf, std::move(y0), std::forward<Ode>(f), std::forward<Params>(params), default_options());
  }

  /** \brief Convenience overload: no parameters, using default Options.
   *  \ingroup rkf78_api
   */
  template <typename State, typename Ode>
  Results integrate(const double t0, const double tf, State y0, Ode &&f)
  {
    return integrate(t0, tf, std::move(y0), std::forward<Ode>(f), default_options());
  }
}