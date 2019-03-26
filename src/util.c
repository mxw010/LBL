//
//  util.c
//
//
//  Created by Wang, Meng on 4/30/18.
//  utility functions for LBL package:
// calcuate the frequency a(h)
// generate double exponential random numbers
// generate dirichlet random numbers
// function for sum and min.
// function updating lambda in MCMC.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "util.h"

//update the value of lambda;
double update_lambda(double *beta, double a, double b, int x_length) {

    double lambda, beta_abs[x_length];
    int i;

    for (i=0; i<x_length; ++i)
        beta_abs[i] = fabs(beta[i]);

    GetRNGstate();
    //Rprintf("slaps\n");
    lambda=rgamma((double) a+x_length, 1/(sum(beta_abs, x_length)+b));
    PutRNGstate();

    return lambda;
}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, f, and d */

double calc_den_post(double *beta, double *freq, int **uniq_x_mat, int x_length, double D, int n_cases, int tot_hap) {
    int i, k;
    double den=0, term[tot_hap];
    for(i=0;i<tot_hap;++i) {
        term[i]=calc_a(freq, uniq_x_mat[i], x_length, D);
        for(k=0;k<x_length;++k) {
            term[i]=term[i]*exp(uniq_x_mat[i][k]*beta[k]);
        }
    }
    den=(n_cases)*log(sum(term, tot_hap));
    return den;
}

/* function to find sum of real numbers */
double sum(double *x, int n) {
    double sum=0.0;
    int i;

    for (i=0; i<n ; ++i)
        sum = sum + x[i];

    return sum;
}

/* function to calculate min. of an array of numbers of length n */

double find_min(double *arr, int n){
    int i;
    double min=arr[0];
    for(i=1;i<n; ++i) {
        if(min > arr[i])
            min = arr[i];
    }
    return min;
}

/* function to generate from double exponential distribution */

double gen_double_exp(double mean, double SD) {
    double x, b;

    GetRNGstate();
    //Rprintf("slaps\n");
    x = runif(0,1);
    b = SD / sqrt(2);
    PutRNGstate();
    //printf("in generating double exponetial, runif = %f \n", x);

    //Rprintf("gen_double, gen_exp=%f, x=%f, mean=%f, SD=%f\n",x,gen_exp,mean,SD);

    if(x > 0.5)
    return mean - b * log(2*(1-x));
    else
    return mean + b * log(2*x);

}

/* function to generate from Dirichet distribution */

void dirichlet(double *param, int dim, double *gen_sample, int n) {
    int i;
    double gen_gamma[dim], sum_gamma;

    GetRNGstate();
    //Rprintf("slaps\n");
    for (i=0; i<dim; ++i) {
        if(param[i]<=0) {
            Rprintf("Warning 3: iter %d,%d-th parameter, param=%lf\n", n, i+1,param[i]);
        }

        gen_gamma[i] = rgamma(param[i], 1);

        if (gen_gamma[i]<=0.000001) {
            Rprintf("Warning 4: gen_gamma=%lf, param=%lf, gen_gamma is reset as 0.000001\n", gen_gamma[i], param[i]);
            gen_gamma[i]=0.000001;
        }
    }

    sum_gamma = sum(gen_gamma, dim);
    for (i=0; i<dim; ++i) {
        gen_sample[i] = gen_gamma[i]/sum_gamma;
    }
    PutRNGstate();
}

/* Calculate a(F) that is in the denominator of the likelihood*/

double calc_a(double *freq, int *per_xz, int x_length, double D) {
    int i, j, num_haplo=0;
    double a;

    for (i=0; i < x_length; ++i) {
        if (per_xz[i]==1) {

            a = 2*(1-D)*freq[i];
            ++num_haplo;

            for (j=i+1; j<x_length; ++j) {
                if (per_xz[j]==1) {
                    a = a*freq[j];
                    ++num_haplo;
                    break;
                }
            }

            if (num_haplo==2) {
                break;
            }
            else {
                a = a*freq[x_length];
                ++num_haplo;
                break;
            }
        }
        else if (per_xz[i]==2) {
            a = D*freq[i]+(1-D)*pow(freq[i],2);
            num_haplo = num_haplo+2;
            break;
        }
    }

    if (num_haplo != 2) {
        a = D*freq[x_length] + (1-D)*pow(freq[x_length], 2);
    }

    return a;
}

void getWts(int *size, int *ID, double *wts, double *num_prob) {
  int i, j;

  for (i = 0; i < size[0]; i++) {
    double sum_num_prob = 0;
    for (j = 0; j < size[0]; j++) {
      if (ID[j] == ID[i])
        sum_num_prob += num_prob[j];
    }
    wts[i] = num_prob[i]/sum_num_prob;
  }

}
