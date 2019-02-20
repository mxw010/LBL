#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "util.h"

void update_z1(double *beta, double *freq, double D, int per_y_new, int *per_xzc, int *per_xzu, int **per_xchaplo, int **per_xuhaplo, int num_haplo_per, int x_length);

void update_z2(double *beta, double *freq, double D, int per_y_new, int *per_xz, int **per_xhaplo, int num_haplo_per, int x_length);

//void update_freq(double *freq, double *beta, double D, int x_length, int N1, int N2, int n_cases1, int n_cases2,
//    int tot_uniq_x, int **xz1, int **xz2, int n, double *prop_freq_with_last);

void update_freq_combined(double *freq, double *beta, double D, int x_length, int N1, int N2, int n_cases1, int n_cases2, int tot_uniq_x, int **xz1, int **xz2, int n);

void update_beta_combined(int **xz1, int **xz2, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int *y_new1, int *y_new2, int n_cases1, int n_cases2, int tot_uniq_x, int N1, int N2);

double update_D_combined(double *freq, double *beta, double D, int x_length, int N1, int N2, int n_cases1, int n_cases2,
    int tot_uniq_x, int **xz1, int **xz2);


int **uniq_x_mat, *y_new1, *y_new2, total_freq_update=0;
void cLBLmcmc(int *x1, int *x2, int *tot_hap1, int *tot_hap2, int *y1, int *y2, int *N1, int *N2, int *num_haplo_id1, int *num_haplo_id2, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x, int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out) {

  int i, j, k, l, m, n, ***xhaplo1, ***xhaplo2, **xz1, **xz2, x1_mat[*tot_hap1][*x_length],
      x2_mat[*tot_hap2][*x_length], n_cases1=0, n_cases2=0, which_beta, it=0, it1=0, it2=0, heat_it;
  //storing the proposed values
  //double * prop_freq_with_last;

  //prop_freq_with_last=calloc(*x_length+1,sizeof(double));

  //FILE* f_freq;
  //    FILE* f_prop;
  Rprintf("running cLBL...\n");
  y_new1=calloc(*N1, sizeof(int));
  y_new2=calloc(*N2, sizeof(int));
  xhaplo1=calloc(*N1, sizeof(int*));
  xhaplo2=calloc(*N2, sizeof(int*));
  for(i=0;i<*N1;++i)
    xhaplo1[i]=calloc(num_haplo_id1[i], sizeof(int*));
  for(i=0;i<*N2;++i)
    xhaplo2[i]=calloc(num_haplo_id2[i], sizeof(int*));
  for(i=0;i<*N1;++i) {
    for(j=0;j<num_haplo_id1[i];++j) {
      xhaplo1[i][j]=calloc(*x_length, sizeof(int));
    }
  }
  for(i=0;i<*N2;++i) {
    for(j=0;j<num_haplo_id2[i];++j) {
      xhaplo2[i][j]=calloc(*x_length, sizeof(int));
    }
  }
  /* separating x vector from R as x matrix */
  l=0;
  for(j=0;j<*x_length;++j) {
    for(i=0;i<*tot_hap1;++i) {
      x1_mat[i][j]=x1[l];
      ++l;
    }
  }
  l=0;
  for(j=0;j<*x_length;++j) {
    for(i=0;i<*tot_hap2;++i) {
      x2_mat[i][j]=x2[l];
      ++l;
    }
  }
  xz1=calloc(*N1, sizeof(int*));
  xz2=calloc(*N2, sizeof(int*));
  for(i=0;i<*N1;++i)
    xz1[i]=calloc(*x_length, sizeof(int));
  for(i=0;i<*N2;++i)
    xz2[i]=calloc(*x_length, sizeof(int));
  /* separating haplotypes per person from the x matrix and assigning z (missing haplotype) for persons with only one compatible haplotype */
  /* family */
  l=0;
  for(i=0;i<*N1;++i) {
    y_new1[i]=y1[l];
    n_cases1+=y_new1[i];
    for(j=0;j<num_haplo_id1[i];++j) {
      m=0;
      for(k=0;k<*x_length;++k) {
        xhaplo1[i][j][k]=x1_mat[l][m];
        if(num_haplo_id1[i]==1) xz1[i][k]=x1_mat[l][m];
        else xz1[i][k]=0;
        ++m;
      }
      ++l;
    }
  }
  /* case-control  */
  l=0;
  for(i=0;i<*N2;++i) {
    y_new2[i]=y2[l];
    n_cases2+=y_new2[i];
    for(j=0;j<num_haplo_id2[i];++j) {
      m=0;
      for(k=0;k<*x_length;++k) {
        xhaplo2[i][j][k]=x2_mat[l][m];
        if(num_haplo_id2[i]==1) xz2[i][k]=x2_mat[l][m];
        else xz2[i][k]=0;
        ++m;
      }
      ++l;
    }
  }
  /* separating unique x vector (haplotypes) from R as unique x matrix (to be used in denominator calculation) */
  uniq_x_mat=calloc(*tot_uniq_x, sizeof(int*));
  for(i=0;i<*tot_uniq_x;++i)
    uniq_x_mat[i]=calloc(*x_length, sizeof(int));
  l=0;
  for(i=0;i<*tot_uniq_x;++i) {
    for(j=0;j<*x_length;++j) {
      uniq_x_mat[i][j]=uniq_x[l];
      ++l;
    }
  }
  /********************** start MCMC here ***************************/


  /*exactly what is going on with the frequency??*/
  //f_freq=fopen("freq.txt","r");
  //f_prop=fopen("freq_prop","r");

  heat_it=0; /* tracks # of heating iterations */
  for(n=0;n<*NUM_IT;++n) {
    heat_it=heat_it+1;
    if(heat_it==101) heat_it=1;
    /* assigning or updating z (missing haplotype) for persons with more than one compatible haplotype */
    /* for family */
    for(i=0;i<n_cases1;++i) {
      if(num_haplo_id1[i]>1) {
        update_z1(beta, freq, *D, y_new1[i], xz1[i], xz1[i+n_cases1], xhaplo1[i], xhaplo1[i+n_cases1],
            num_haplo_id1[i], *x_length);
      }
    }
    /* for case-control */
    for(i=0;i<*N2;++i) {
      if(num_haplo_id2[i]>1) {
        update_z2(beta, freq, *D, y_new2[i], xz2[i], xhaplo2[i], num_haplo_id2[i], *x_length);
      }
    }
    /* update beta parameters */
    for(i=0;i<*x_length;++i) {
      which_beta=i;
      update_beta_combined(xz1, xz2, beta, which_beta, *lambda, freq, *D, *x_length, y_new1, y_new2, n_cases1, n_cases2,
          *tot_uniq_x, *N1, *N2);
    }
    /* update lambda */
    *lambda=update_lambda(beta, *a, *b, *x_length);


          for (i=0;i<*x_length;i++) {
    //  fprintf(f_freq,"%e ", freq[i]);
    }
    //fprintf(f_freq, "%e\n", freq[*x_length+1]);


    /* update frequencies and D */
    update_freq_combined(freq, beta, *D, *x_length, *N1, *N2, n_cases1, n_cases2, *tot_uniq_x, xz1, xz2,n);

        //  for (i=0;i<*x_length;i++) {
    //  fprintf(f_prop,"%e ", prop_freq_with_last[i]);
    //}
    //fprintf(f_prop, "%e\n", prop_freq_with_last[*x_length+1]);


    /* update D parameter */
    *D=update_D_combined(freq, beta, *D, *x_length, *N1, *N2, n_cases1, n_cases2, *tot_uniq_x, xz1, xz2);

    /* add some updates every 5000 iterations */
   // if (monitor == 1) {
   //   if ( (n+1) % 5000 == 0) Rprintf("Running iteration %d...\n", n+1);
   // }

    if(n>=*BURN_IN) {
      for(i=0;i<*x_length;++i) {
        beta_out[it]=beta[i];
        ++it;
      }
      lambda_out[it2]=*lambda;
      for(i=0;i<*x_length+1;++i) {
        freq_out[it1]=freq[i];
        ++it1;
      }
      D_out[it2]=*D;
      ++it2;
    }

  //  fclose(f_freq);
  //  fclose(f_prop);
  }
}
void update_z1(double *beta, double *freq, double D, int per_y_new, int *per_xzc, int *per_xzu, int **per_xchaplo,
    int **per_xuhaplo, int num_haplo_per, int x_length) {
  int i, k;
  double prob[num_haplo_per], cum_prob[num_haplo_per], sum_prob, x, ac[num_haplo_per], au[num_haplo_per];
  for(i=0;i<num_haplo_per;++i) {
    ac[i]=calc_a(freq, per_xchaplo[i], x_length, D);
    au[i]=calc_a(freq, per_xuhaplo[i], x_length, D);
    prob[i]=ac[i]*au[i];
    for(k=0;k<x_length;++k) {
      prob[i]=prob[i]*exp(per_xchaplo[i][k]*beta[k]);
    }
  }
  sum_prob=sum(prob, num_haplo_per);
  for(i=0;i<num_haplo_per;++i) {
    prob[i]=prob[i]/sum_prob;
  }
  cum_prob[0]=prob[0];
  for(i=1;i<num_haplo_per;++i) {
    cum_prob[i]=cum_prob[i-1]+prob[i];
  }
  GetRNGstate();
  x=runif(0, 1);
  PutRNGstate();
  if(x<cum_prob[0]) {
    for(k=0;k<x_length;++k) {
      per_xzc[k]=per_xchaplo[0][k];
      per_xzu[k]=per_xuhaplo[0][k];
    }
  }
  for(i=1;i<num_haplo_per;++i) {
    if(x>cum_prob[i-1]&&x<cum_prob[i]) for(k=0;k<x_length;++k) {
      per_xzc[k]=per_xchaplo[i][k];
      per_xzu[k]=per_xuhaplo[i][k];
    }
  }
}
void update_z2(double *beta, double *freq, double D, int per_y_new, int *per_xz, int **per_xhaplo, int num_haplo_per,
    int x_length) {
  int i, k;
  double prob[num_haplo_per], cum_prob[num_haplo_per], sum_prob, x, a[num_haplo_per];
  for(i=0;i<num_haplo_per;++i) {
    a[i]=calc_a(freq, per_xhaplo[i], x_length, D);
    prob[i]=a[i];
    if(per_y_new==1) {
      for(k=0;k<x_length;++k) {
        prob[i]=prob[i]*exp(per_xhaplo[i][k]*beta[k]);
      }
    }
  }
  sum_prob=sum(prob, num_haplo_per);
  for(i=0;i<num_haplo_per;++i) {
    prob[i]=prob[i]/sum_prob;
  }
  cum_prob[0]=prob[0];
  for(i=1;i<num_haplo_per;++i) {
    cum_prob[i]=cum_prob[i-1]+prob[i];
  }
  GetRNGstate();
  x=runif(0, 1);
  PutRNGstate();
  if(x<cum_prob[0]) {
    for(k=0;k<x_length;++k)
      per_xz[k]=per_xhaplo[0][k];
  }
  for(i=1;i<num_haplo_per;++i) {
    if(x>cum_prob[i-1]&&x<cum_prob[i]) for(k=0;k<x_length;++k)
      per_xz[k]=per_xhaplo[i][k];
  }
}
void update_beta_combined(int **xz1, int **xz2, double *beta, int which_beta, double lambda, double *freq, double D,
    int x_length, int *y_new1, int *y_new2, int n_cases1, int n_cases2, int tot_uniq_x, int N1, int N2) {
  double beta_new, g_old, g_new, f_old, f_new, beta_new_vec[x_length], SD, accept_prob;
  int i, x;
  beta_new=gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
  g_old=-lambda*fabs(beta[which_beta]);
  g_new=-lambda*fabs(beta_new);
  for(i=0;i<N1;++i) {
    if(y_new1[i]==1) {
      g_old+=xz1[i][which_beta]*beta[which_beta];
      g_new+=xz1[i][which_beta]*beta_new;
    }
  }
  for(i=0;i<N2;++i) {
    if(y_new2[i]==1) {
      g_old+=xz2[i][which_beta]*beta[which_beta];
      g_new+=xz2[i][which_beta]*beta_new;
    }
  }
  g_old=g_old-calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases1+n_cases2, tot_uniq_x);
  for(i=0;i<x_length;++i)
    beta_new_vec[i]=beta[i];
  beta_new_vec[which_beta]=beta_new;
  g_new=g_new-calc_den_post(beta_new_vec, freq, uniq_x_mat, x_length, D, n_cases1+n_cases2, tot_uniq_x);
  SD=sqrt(fabs(beta_new));
  f_old=exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
  SD=sqrt(fabs(beta[which_beta]));
  f_new=exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
  accept_prob=exp(g_new-g_old)*f_old/f_new;
  if(accept_prob>1) beta[which_beta]=beta_new;
  else{
    GetRNGstate();
    x=rbinom(1, accept_prob);
    PutRNGstate();
    if(x==1) beta[which_beta]=beta_new;
  }
}

void update_freq_combined(double *freq, double *beta, double D, int x_length, int N1, int N2, int n_cases1, int n_cases2, int tot_uniq_x, int **xz1, int **xz2, int n) {
  int i, C=10000,update=0;
  //double prop_freq_with_last[x_length+1], gen_gamma[x_length+1], sum_gamma, g_old, g_new, accept_prob=0, freq_last,
  double prop_freq_with_last[x_length+1], g_old, g_new, g_new_temp, accept_prob=0,
      f_old, f_new, b_old[x_length+1], b_new[x_length+1], min_f_old, min_f_new;
  GetRNGstate();
  for(i=0;i<x_length+1;++i) {
    b_old[i]=freq[i]*C;
  }

  /*prosing new frequencies, stored in prop_freq_with_last)*/
  dirichlet(b_old, x_length+1, prop_freq_with_last,n);


  min_f_old=find_min(freq, x_length+1);
  /* check if the constraint -min fk/(1-min fk) < d is satisfied */
  min_f_new=find_min(prop_freq_with_last, x_length+1);
    /* print out a more detailed warning */
  if(-min_f_new/(1-min_f_new)<D) {
    /* needed in acceptance prob. computation */
    for(i=0;i<x_length+1;++i) {
      b_new[i]=prop_freq_with_last[i]*C;
      if(b_new[i]<=0) {
        Rprintf("Warning 1: iteration %d, pos = %d, prop_freq_with_last=%e, b_new = %lf\n",n,i+1, prop_freq_with_last[i], b_new[i]);
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
    g_old=-calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases1+n_cases2, tot_uniq_x)+log(1-min_f_old);

    g_new_temp=calc_den_post(beta, prop_freq_with_last, uniq_x_mat, x_length, D, n_cases1+n_cases2, tot_uniq_x);

    g_new=-g_new_temp+log(1-min_f_new);
    if(g_new_temp<-pow(10, 10)) Rprintf(
        "Warning 2: iternation %d, calc_den_new = %e, g_old = %lf, g_new = %lf\n",n, g_new_temp, g_old, g_new);
    for(i=0;i<n_cases1;++i) {
      g_old+=log(calc_a(freq, xz1[i], x_length, D))+log(calc_a(freq, xz1[n_cases1+i], x_length, D));
      g_new+=log(calc_a(prop_freq_with_last, xz1[i], x_length, D))
          +log(calc_a(prop_freq_with_last, xz1[n_cases1+i], x_length, D));
    }
    for(i=0;i<N2;++i) {
      g_old+=log(calc_a(freq, xz2[i], x_length, D));
      g_new+=log(calc_a(prop_freq_with_last, xz2[i], x_length, D));
    }
    /* calculate f(f*|f^(t)) = f_new and f(f^(t)|f*) = f_old */
    f_old=lgammafn(C);
    f_new=f_old;
    for(i=0;i<x_length+1;++i) {
      f_old+=(b_new[i]-1)*log(freq[i])-lgammafn(b_new[i]);
      f_new+=(b_old[i]-1)*log(prop_freq_with_last[i])-lgammafn(b_old[i]);
    }
    accept_prob=exp(g_new-g_old+f_old-f_new);
    if(-min_f_new/(1-min_f_new)>D) accept_prob=0;

     if(accept_prob>1) update=1;
     else update=rbinom(1, accept_prob);
    if(update==1) {
      for(i=0;i<x_length+1;++i)
        freq[i]=prop_freq_with_last[i];
      ++total_freq_update;
    }
  }
  if(-min_f_old/(1-min_f_old)>D) {
    Rprintf("error in updating f: Min f_old and D constraint violated min_f_old=%f\n", min_f_old);
    for(i=0;i<x_length+1;++i)
      Rprintf("%f\n", freq[i]);
    Rprintf("%f\n", D);
    error("LBL has run into an error and has stopped\n");
  }
  if((-min_f_new/(1-min_f_new)>D)&(accept_prob>0)) {
    Rprintf("error in updating f: Min f_new and D constraint violated ");
    error("LBL has run into an error and has stopped\n");
  }
  PutRNGstate();
}
double update_D_combined(double *freq, double *beta, double D, int x_length, int N1, int N2, int n_cases1, int n_cases2,
    int tot_uniq_x, int **xz1, int **xz2) {
  int i, update=0;
  double prop_D, accept_prob, g_old, g_new, min_f, delta=0.05, lower, upper, f_old, f_new;
  GetRNGstate();
  min_f=find_min(freq, x_length+1);
  lower=D-delta;
  upper=D+delta;
  if(lower<-min_f/(1-min_f)) lower=-min_f/(1-min_f);
  if(upper>1) upper=1;
  prop_D=runif(lower, upper);
  g_old=-calc_den_post(beta, freq, uniq_x_mat, x_length, D, n_cases1+n_cases2, tot_uniq_x);
  g_new=-calc_den_post(beta, freq, uniq_x_mat, x_length, prop_D, n_cases1+n_cases2, tot_uniq_x);
  for(i=0;i<n_cases1;++i) {
    g_old+=log(calc_a(freq, xz1[i], x_length, D))+log(calc_a(freq, xz1[i+n_cases1], x_length, D));
    g_new+=log(calc_a(freq, xz1[i], x_length, prop_D))+log(calc_a(freq, xz1[i+n_cases1], x_length, prop_D));
  }
  for(i=0;i<N2;++i) {
    g_old+=log(calc_a(freq, xz2[i], x_length, D));
    g_new+=log(calc_a(freq, xz2[i], x_length, prop_D));
  }
  f_new=1/(upper-lower);
  lower=prop_D-delta;
  upper=prop_D+delta;
  if(lower<-min_f/(1-min_f)) lower=-min_f/(1-min_f);
  if(upper>1) upper=1;
  f_old=1/(upper-lower);
  accept_prob=exp(g_new-g_old)*f_old/f_new;
  if(-min_f/(1-min_f)>D||-min_f/(1-min_f)>prop_D) {
    Rprintf("error in updating D:  Min f and D constraint violated ");
    error("LBL has run into an error and has stopped\n");
  }
  if(accept_prob>1) update=1;
  else update=rbinom(1, accept_prob);
  if(update==1) return prop_D;
  else return D;
  PutRNGstate();
}

