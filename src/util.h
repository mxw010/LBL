//
//  util.h
//  
//
//  Created by Wang, Meng on 4/30/18.
//

#ifndef util_h
#define util_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

#endif /* util_h */

double calc_a(double *freq, int *per_xz, int x_length, double D);

double gen_double_exp(double mean, double SD);

void dirichlet(double *param, int dim, double *gen_sample, int n);

double sum(double *x, int n);

double find_min(double *arr, int n);

double calc_den_post(double *beta, double *freq, int **uniq_x_mat, int x_length, double D, int n_cases, int tot_hap);

double update_lambda(double *beta, double a, double b, int x_length);

void getWts(int *size, int *ID, double *wts, double *num_prob);