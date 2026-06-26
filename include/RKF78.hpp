/*
 * Created on Wed Jun 24 2026
 *
 * Copyright (c) 2026 David Stanley
 */

#include <array>
#include <vector>
#include <cmath>
#include <span>

#pragma once // Include guard to prevent multiple inclusions of this header file

namespace RKF78
{

  template <typename State>
  struct Results
  {
    double *y;        // Pointer to the array of computed values
    double *t;        // Pointer to the array of time points
    int rejected = 0; // Number of steps taken during integration
    int accepted = 0; // Number of rejected steps during integration
    int fevals = 0;   // Number of function evaluations performed
    State yf;         // Final state vector after integration
  };

  struct Options
  {
    double h0;   // initial step size
    double atol; // absolute tolerance
    double rtol; // relative tolerance
    std::span<const double> tout;
  };

  constexpr int stages = 13;

  // a_k: stage time fractions
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

  // beta_{k,lambda}: lower-triangular RK matrix
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

  // function type for first order ODE system
  template <typename State, typename Params>
  using ODEFunction = void (*)(const double t, const State &y, State &dy, const Params &params);

  // integrate template function
  template <typename State, typename Params>
  Results<State> integrate(const double t0, const double tf, State y0, ODEFunction<State, Params> f, const Params &params, Options options)
  {
    double &atol = options.atol;
    double &rtol = options.rtol;
    double &h = options.h0; // initial step size

    // Initialize the deriv variables for the Runge-Kutta-Fehlberg 7(8) method
    std::array<State, stages> k; // Array to hold the derivatives at each stage

    // Initialize the solution storage arrays (need to be resizeable)
    std::vector<double> t_values; // Vector to hold the time points
    std::vector<State> y_values;  // Vector to hold the computed values

    int dims = y0.size();              // Get the size of the state vector
    int prop_dir = (tf > t0) ? 1 : -1; // Determine the propagation direction based on the time interval

    // construct the results struct to return
    Results<State> results;

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
            y_temp[d] += h * bk[i][j] * k[j][d];
          }
        }
        f(t + ak[i] * h, y_temp, k[i], params);
        fevals++; // Increment the function evaluation counter for each stage
        // check for NaN or Inf values in the computed derivatives
        for (int d = 0; d < dims; ++d)
        {
          if (!std::isfinite(k[i][d]))
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
        y8[d] += h * ( 34.0 / 105.0 * k[5][d] + 9.0 / 35.0 * k[6][d] + 9.0 / 35.0 * k[7][d] + 9.0 / 280.0 * k[8][d] + 9.0 / 280.0 * k[9][d] + 41.0 / 840.0 * k[11][d] +  41.0 / 840.0 * k[12][d]);
      } // end computing 8th order solution

      // Compute the truncation error estimate for the step
      State error_estimate;
      for (int d = 0; d < dims; ++d)
      {
        error_estimate[d] = h * 41.0 / 840.0 * (k[0][d] + k[10][d] - k[11][d] - k[12][d]);
      }

      double err2 = 0.0; // error scaled by the tolerance

      for (int d = 0; d < dims; ++d)
      {
        double scale = atol + rtol * std::max(std::abs(y0[d]), std::abs(y8[d]));
        double ratio = error_estimate[d] / scale;
        err2 += ratio * ratio;
      }

      double error_norm = std::sqrt(err2 / dims);
      double h_new; // New step size for the next iteration
      if (error_norm == 0.0)
      {
        h_new = 2.0 * h; // If the error is zero, increase the step size
      }
      else{
      h_new = 0.9 * h * std::pow(1.0 / error_norm, 1.0 / 8.0);
      }
      if (error_norm <= 1.0) // Check if the error is within the tolerance
      {
        // Accept the step
        y0 = y8; // Update the state vector to the new value
        t += h;  // Update the time to the new value
        h = h_new;
        istep++; // Increment the step counter

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
    results.yf = y0; // Set the final state vector after integration

    return results; // Return the results of the integration
  }
}