#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "eigen.h"

double design_vol(double m_target, double rho){
  return m_target / rho;
}

double compute_vol(double*params, double width){ 
  double r1 = params[0];
  double r2 = params[0] + width;
  double e = params[2];
  double l = params[3];
  double w = r2-r1;

  double volume = w*(M_PI*(r2*r2-r1*r1)/2 + e*w + 2*w*(fabs(l-r2)));

  return volume;
}

// Implement the bisection method
double bisection(double a, double b, double tol, double * params, double v_target)
{
    double c = 0.0;
    while (fabs(b-a) > tol)
    {
        c = (a + b) / 2.0;
        if ((compute_vol(params, c) - v_target)  == 0.0) // If we find the root exactly
        {
            return c;
        }
        else if ((compute_vol(params, a) - v_target)*(compute_vol(params, c) - v_target) < 0) // Root is in [a, c]
        {
            b = c;
        }
        else // Root is in [c, b]
        {
            a = c;
        }
    }
    return c;
}

double find_param(double *params, int to_find, double upper_bound, double lower_bound, double tol, double ** u, int * n_nodes, double v_target, double * freq, double freq_target)
{
  // f(low) < freq_target < f(up)
  double width;

  params[to_find] = lower_bound;
  width = bisection(0, 1, 1e-6, params, v_target);
  params[1] = params[0] + width; 

  double f_low = compute_freq(params, NULL, NULL);

  params[to_find] = upper_bound;
  width = bisection(0, 1, 1e-6, params, v_target);
  params[1] = params[0] + width; 

  double f_up = compute_freq(params, NULL, NULL);

  double c = (lower_bound + upper_bound) / 2;
  params[to_find] = c;
  width = bisection(0, 1, 1e-6, params, v_target);
  params[1] = params[0] + width; 

  double f_c = compute_freq(params, NULL, NULL);

  if (f_up<freq_target  || f_low>freq_target){
    printf("freq target not inside given interval");
    return 0.0;
  }

  int nbr_ite = 0;

  while (fabs(f_c - freq_target) > tol)
  {
    printf("error: %.3e\n", f_c - freq_target);
    if (f_c - freq_target < 0)
    {
      f_low = f_c;
      lower_bound = c;
    }
    else
    {
      f_up = f_c;
      upper_bound = c;
    }

    c = (lower_bound + upper_bound) / 2;
    printf("c = %.3e\n", c);
    params[to_find] = c;
    width = bisection(0, 1, 1e-6, params, v_target);
    params[1] = params[0] + width; 

    f_c = compute_freq(params, NULL, NULL);
    nbr_ite++;
  }

  f_c = compute_freq(params, u, n_nodes);
  printf("Optimization done!\n");

  double rho = 3000;
  double volume = compute_vol(params, width);
  double mass = volume * rho;

  printf("Tuning Fork mass: %.3e\n", mass);
  printf("Number of iterations needed: %d\n",nbr_ite);
  printf("Distance from target frequency: %.3eHz\n", (f_c - freq_target));

  *freq = f_c;

  return c;
}