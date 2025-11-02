#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

double design_vol(double m_target, double rho);
double find_param (double* params, int to_find, double upper_bound, double lower_bound, double tol, double ** u, int * n_nodes, double v_target, double * freq, double freq_target);

#endif
