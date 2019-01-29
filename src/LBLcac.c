#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*#include <lib.h>*/
#include <R.h>
#include <Rmath.h>
#include "util.h"

/*Meng 03/16/18: increased readability, added more error/warning messages for debugging */

void update_z_cc(double *beta, double *freq, double D, int per_y_new, int *per_xz, int **per_xhaplo, int num_haplo_per, int x_length);


void update_beta_cc(int **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int *y_new, int n_cases, int tot_uniq_x, int N);

void update_freq_cc(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz, int n_iter);

double update_D_cc(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz, int n_iter);

int **uniq_x_mat, *y_new;
//total_freq_update=0;

void LBLmcmc(int *x, int *tot_hap, int *y, int *N, int *num_haplo_id, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x,
  int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out, int monitor) {

  int i, j, k, l, m, n, ***xhaplo, **xz, x_mat[*tot_hap][*x_length], n_cases=0, which_beta, it=0, it1=0, it2=0, heat_it;
  y_new = calloc(*N, sizeof(int));
  xhaplo = calloc(*N, sizeof(int*));
    
  Rprintf("running LBL...\n");
  for (i =0; i<*N; ++i)
    xhaplo[i] = calloc(num_haplo_id[i], sizeof(int*));

  for (i =0; i<*N; ++i) {
    for (j = 0; j < num_haplo_id[i]; ++j) {
      xhaplo[i][j] = calloc(*x_length, sizeof(int));
      }
  }

  //tracking frequency proposal for debugging 
  // double * prop_freq_with_last;
  // prop_freq_with_last=calloc(*x_length+1,sizeof(double));

  /* separating x vector from R as x matrix */
  l=0;
  for (j=0; j<*x_length; ++j) {
    for (i=0; i<*tot_hap; ++i) {
      x_mat[i][j] = x[l];
      ++l;
    }
  }

  xz = calloc(*N, sizeof(int*));
  for (i=0; i<*N; ++i) {
    xz[i] = calloc(*x_length, sizeof(int));
  }
  
  /* separating haplotypes per person from the x matrix and assigning z (missing haplotype) for persons with only one compatible haplotype */
  l = 0;
  for (i =0; i<*N; ++i) {
    y_new[i] = y[l];
    n_cases += y_new[i];
    for (j = 0; j < num_haplo_id[i]; ++j) {
      m = 0;
      for (k = 0; k < *x_length; ++k) {
        xhaplo[i][j][k] = x_mat[l][m];
        if (num_haplo_id[i]==1)
		      xz[i][k] = x_mat[l][m];
	      else
		      xz[i][k] = 0;
	      ++m;
      }
      ++l;
    }
  }

   /* separating unique x vector (haplotypes) from R as unique x matrix (to be used in denominator calculation) */

  uniq_x_mat = calloc(*tot_uniq_x, sizeof(int*));
  for (i=0; i<*tot_uniq_x; ++i)
    uniq_x_mat[i] = calloc(*x_length, sizeof(int));

  l=0;
  for (i=0; i<*tot_uniq_x; ++i) {
    for (j=0; j<*x_length; ++j) {
      uniq_x_mat[i][j] = uniq_x[l];
	     ++l;
     }
   }

  /********************** start MCMC here ***************************/
  heat_it = 0; /* tracks # of heating iterations */
  for (n=0; n<*NUM_IT; ++n) {
    heat_it = heat_it+1;
    if (heat_it == 101) heat_it = 1;
      
      //Rprintf("iteration %d:\n", n);
      
     
    //allows user to interupt every 1024 iterations
    if (n % 1024 == 0) R_CheckUserInterrupt();
    /* assigning or updating z (missing haplotype) for persons with more than one compatible haplotype */

    for (i =0; i<*N; ++i) {
      if (num_haplo_id[i]>1) {
		    update_z_cc(beta, freq, *D, y_new[i], xz[i], xhaplo[i], num_haplo_id[i], *x_length);
	    }
	   }

    /* update beta parameters */
    
    for (i=0; i<*x_length; ++i) {
      which_beta=i;
      update_beta_cc(xz, beta, which_beta, *lambda, freq, *D, *x_length, y_new, n_cases, *tot_uniq_x, *N);
	   }
      
      
    /* update lambda */

    *lambda = update_lambda(beta, *a, *b, *x_length);
      
      //Rprintf("updated lambda=%f\n", lambda);

    /* update frequencies and D */

    //update_freq(freq, beta, *D, *x_length, *N, n_cases, *tot_uniq_x, xz, n, prop_freq_with_last);
      update_freq_cc(freq, beta, *D, *x_length, *N, n_cases, *tot_uniq_x, xz, n);

//      Rprintf("updated beta=");
//      for (i=0;i<*x_length;i++) {
//          Rprintf(" %f,", freq[i]);
//      }
//      Rprintf("\n");

      
    /* update D parameter */

    *D = update_D_cc(freq, beta, *D, *x_length, *N, n_cases, *tot_uniq_x, xz, n);

        
    /* add some updates every 5000 iterations */
    if (monitor == 1) {
      if ( (n+1) % 5000 == 0) Rprintf("Running iteration %d...\n", n+1);
    }
      
    if (n >= *BURN_IN) {
      for (i=0; i<*x_length; ++i) {
        beta_out[it] = beta[i];
        ++it;
      }

      lambda_out[it2] = *lambda;

      for (i=0; i<*x_length+1; ++i) {
        freq_out[it1] = freq[i];
	      ++it1;
	     }
	     D_out[it2] = *D;
	     ++it2;
     }
   }
 }


void update_z_cc(double *beta, double *freq, double D, int per_y_new, int *per_xz, int **per_xhaplo, int num_haplo_per, int x_length) {

  int i, k;
  double prob[num_haplo_per], cum_prob[num_haplo_per], sum_prob, x, a[num_haplo_per];

  for (i=0; i<num_haplo_per; ++i) {
    a[i] = calc_a(freq, per_xhaplo[i], x_length, D);
    prob[i] = a[i];

    if (per_y_new == 1) {
      for (k=0; k < x_length; ++k) {
        prob[i] = prob[i]*exp(per_xhaplo[i][k]*beta[k]);
      }
    }
  }

  sum_prob = sum(prob, num_haplo_per);

  for (i=0; i<num_haplo_per; ++i) {
    prob[i] = prob[i]/sum_prob;
  }

  
  cum_prob[0]= prob[0];
  for (i=1; i<num_haplo_per; ++i) {
    cum_prob[i] = cum_prob[i-1]+prob[i];
  }

  GetRNGstate();
  //Rprintf("slaps\n");
  x=runif(0,1);
  PutRNGstate();

  if (x < cum_prob[0]) {
    for (k=0; k<x_length; ++k)
      per_xz[k] = per_xhaplo[0][k];
  }

  for (i=1; i<num_haplo_per; ++i) {
    if (x > cum_prob[i-1] && x < cum_prob[i]) {
      for (k=0; k<x_length; ++k)
        per_xz[k] = per_xhaplo[i][k];
    }
  }
}

void update_beta_cc(int **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int *y_new, int n_cases, int tot_uniq_x, int N) {
  
  double beta_new, g_old, g_new, f_old, f_new, beta_new_vec[x_length], SD, accept_prob;
  int i, x;

  beta_new = gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
  g_old = -lambda*fabs(beta[which_beta]);
  g_new = -lambda*fabs(beta_new);
    
  for (i=0; i<N; ++i) {
    if (y_new[i]==1) {
      g_old += xz[i][which_beta]*beta[which_beta];
      g_new += xz[i][which_beta]*beta_new;
    }
  }

  g_old = g_old-calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases, tot_uniq_x);

  //Rprintf("i= %d, prosed beta= %f,g_old = %f, g_new= %f\n",which_beta, beta_new,g_old, g_new);
  for (i=0; i<x_length; ++i) {
    beta_new_vec[i] = beta[i];
  }
  
  beta_new_vec[which_beta] = beta_new;

  g_new = g_new-calc_den_post(beta_new_vec, freq, uniq_x_mat, x_length, D, n_cases, tot_uniq_x);

  SD = sqrt(fabs(beta_new));
  f_old = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);

  SD = sqrt(fabs(beta[which_beta]));
  f_new = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
  accept_prob = exp(g_new-g_old)*f_old/f_new;

//Rprintf("g_new2 = %f, f_old = %f, f_new = %f, accept_prob= %f, ", g_new, f_old, f_new, accept_prob);
  if (accept_prob > 1) {
    beta[which_beta] = beta_new;
  }
  else {
    GetRNGstate();
       ////Rprintf("slaps\n");
    x = rbinom(1,accept_prob);
    PutRNGstate();
    if (x == 1) beta[which_beta] = beta_new;
  }
}


void update_freq_cc(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz, int n_iter) {
  
  int i, C=10000,update=0;
  //double prop_freq_with_last[x_length+1], gen_gamma[x_length+1], sum_gamma, g_old, g_new, accept_prob=0, freq_last, alpha=0.1, f_old, f_new, b_old[x_length+1], b_new[x_length+1], min_f_old, min_f_new;
  double prop_freq_with_last[x_length+1], g_old, g_new, g_new_temp, accept_prob=0, f_old, f_new, b_old[x_length+1], b_new[x_length+1], min_f_old, min_f_new;
 
  for (i=0; i<x_length+1; ++i) {
    b_old[i] = freq[i]*C;
  }

  dirichlet(b_old, x_length+1, prop_freq_with_last,n_iter);
  min_f_old = find_min(freq, x_length+1);

  /* check if the constraint -min fk/(1-min fk) < d is satisfied */

  min_f_new = find_min(prop_freq_with_last, x_length+1);

  if (-min_f_new/(1-min_f_new)  < D) {
  /* needed in acceptance prob. computation */

    for (i=0; i<x_length+1; ++i) {
      b_new[i] = prop_freq_with_last[i]*C;
      if(b_new[i]<=0) {
        Rprintf("Warning 1: iteration %d, pos = %d, prop_freq_with_last=%e, b_new = %lf\n",n_iter,i+1, prop_freq_with_last[i], b_new[i]);
        Rprintf("proposed freq=(");
        for(i=0;i<x_length;++i) {
          Rprintf("%e, ", prop_freq_with_last[i]);
        }
        Rprintf("%e)\n",prop_freq_with_last[x_length+1]);
      
        Rprintf("old_freq*C=(");
        for(i=0;i<x_length;++i) {
          Rprintf("%e, ", b_old[i]);
        }
        Rprintf("%e)\n",b_old[x_length+1]);
      
        Rprintf("old_freq=(");
        for(i=0;i<x_length;++i) {
          Rprintf("%e, ", freq[i]);
        }
        Rprintf("%e)\n",freq[x_length+1]);
      }
    }

      /* calculate g(f^(t)) and g(f*) */

    g_old = -calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases, tot_uniq_x)+log(1-min_f_old);

    g_new_temp=calc_den_post(beta, prop_freq_with_last, uniq_x_mat, x_length, D, n_cases, tot_uniq_x);

    g_new = -g_new_temp+log(1-min_f_new);

    if (g_new_temp < -pow(10,10))
      Rprintf("Warning 2: iternation %d, calc_den_new = %e, g_old = %lf, g_new = %lf",n_iter, g_new_temp, g_old, g_new);
        
    for (i=0; i<N; ++i) {
      g_old += log(calc_a(freq, xz[i], x_length, D));
      g_new += log(calc_a(prop_freq_with_last, xz[i], x_length, D));
    }

    /* calculate f(f*|f^(t)) = f_new and f(f^(t)|f*) = f_old */

    f_old = lgammafn(C);
    f_new = f_old;

    for (i=0; i<x_length+1; ++i) {
      f_old += (b_new[i]-1)*log(freq[i]) - lgammafn(b_new[i]);
      f_new += (b_old[i]-1)*log(prop_freq_with_last[i]) - lgammafn(b_old[i]);
    }

    accept_prob = exp(g_new-g_old+f_old-f_new);

    if (-min_f_new/(1-min_f_new)  > D) accept_prob = 0;

     if (accept_prob > 1) {
       update = 1;
     }
     else {
         GetRNGstate();
         ////Rprintf("slaps\n");
       update = rbinom(1, accept_prob);
          PutRNGstate();
     }

    if (update ==1) {
      for (i=0; i<x_length+1; ++i) {
        freq[i] = prop_freq_with_last[i];
      }
      //++total_freq_update;
    }
  }

 /* error messages */
 
  if (-min_f_old/(1-min_f_old)  > D) {
    Rprintf("Iteration %d error: in updating f: Min f_old and D constraint violated min_f_old=%f", n_iter, min_f_old);
    Rprintf("current frequncy=(");
    for (i=0; i<x_length; ++i) {
      Rprintf("%f, ", freq[i]);
    }
    Rprintf("%f). D=%f.\n", freq[x_length],D);  
    error("LBL ran into an error and has stopped\n");
  }

  if (-min_f_new/(1-min_f_new)  > D & accept_prob>0) {
    Rprintf("Iteration %d error: in updating f: Min f_new and D constraint violated min_f_new=%f", n_iter, min_f_new);
    Rprintf("current frequncy=(");
    for (i=0; i<x_length; ++i)
      Rprintf("%f, ", prop_freq_with_last[i]);
    Rprintf("%f). D=%f.\n", prop_freq_with_last[x_length],D);  
    error("LBL ran into an error and has stopped\n");
  }
 
}

double update_D_cc(double *freq, double *beta, double D, int x_length, int N, int n_cases, int tot_uniq_x, int **xz, int n_iter) {
  
  int i, update = 0;
  double prop_D, accept_prob, g_old, g_new, min_f, delta=0.05, lower, upper, f_old, f_new;

  

  min_f = find_min(freq, x_length+1);
  lower = D-delta;
  upper = D+delta;

  if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
  if (upper > 1) upper = 1;
    
    GetRNGstate();
    //Rprintf("slaps\n");

  prop_D = runif(lower, upper);
    PutRNGstate();
    

  g_old = -calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases, tot_uniq_x);
  g_new = -calc_den_post(beta, freq, uniq_x_mat, x_length, prop_D, n_cases, tot_uniq_x);

  for (i=0; i<N; ++i) {
    g_old += log(calc_a(freq, xz[i], x_length, D));
    g_new += log(calc_a(freq, xz[i], x_length, prop_D));
  }

  f_new = 1/(upper-lower);

  lower = prop_D-delta;
  upper = prop_D+delta;
  if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);
  if (upper > 1) upper = 1;

  f_old = 1/(upper-lower);

  accept_prob = exp(g_new-g_old)*f_old/f_new;

  if (-min_f/(1-min_f) > D || -min_f/(1-min_f) > prop_D) {
    Rprintf("Iteration %d error: in updating f: Min f_new and D constraint violated\n", n_iter);
    Rprintf("current frequncy=(");
    for (i=0; i<x_length; ++i)
      Rprintf("%f, ", freq[i]);
    Rprintf("%f). current_d=%f, proposed_d=%f.\n", freq[x_length],D, prop_D);  
    error("LBL has run into an error and has stopped\n");
  }

    
  if (accept_prob > 1) update = 1;
    else {
        GetRNGstate();
        //Rprintf("slaps\n");
        update = rbinom(1, accept_prob);
        PutRNGstate();
    }

  if (update == 1) {
    return prop_D;
  }
  else {
    return D;
  }

}
