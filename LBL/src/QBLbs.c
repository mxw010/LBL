#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libnew.h"
#include "R.h"
#include <Rmath.h>

double GXE_calc_a(double *freq, int *per_freq, double D);
void GXE_update_z(int i, int **per_map, int *per_hz, double *beta, double *freq, double D, double sigmas, double per_y_new,
                  double *per_xz, double **per_xhaplo, int num_haplo_per, int x_length, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *per_match_freq_map_to_uniq, int *match_hz_to_uniq, double *ele, double *xbeta);
double GXE_update_lambda(double *beta, double a, double b, int x_length);
double GXE_update_sigma(double *beta, double asig, double bsig, double **xz, int N, int x_length, double *ele);
void GXE_update_beta(double sigma, double **xz, double *beta, int which_beta, double lambda, double *freq, double D,
                     int x_length, int h_length, double *GXE_y_new, int N, int n_cov, double *ele, double *xbeta);
double GXE_gen_double_exp(double mean, double SD);
void GXE_update_freq(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov,
                     int Con, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *match_hz_to_uniq);
double GXE_update_D(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *match_hz_to_uniq);
void GXE_dirichlet(double *param, int dim, double *gen_sample);
double GXE_sum(double *x, int n);
double GXE_find_min(double *arr, int n);
int accept_count=0;
int **GXE_uniq_map, GXE_total_freq_update=0;/*define unique matrix x,
 response y, */
double *GXE_y_new;
void QBLmcmc(double *x, int *tot_hap, double *y, int *N, int *num_haplo_id, int *x_length, int *h_length,
              int *freq_num, int *uniq_freq_num, int *num_uniq_hap_pair, double *freq, double *D, double *sigmas, double *beta, double *a, double *b, double *asig,
              double *bsig, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out,
              double *sigmas_out, double *freq_out, double *D_out, int *n_cov, double *up_xz, int *CC, double *BDIC_out)
  /*
   *x:        design matrix as a vector(col1',col2'...)
   *tot_hap:  # of rows in design matrix
   *y:        case/control(affected)
   *N:        sum of weight(total # of subject)
   *num_haplo_id: the # of rows for each original ID
   *x_length: # of hyplotypes (m)
   *freq:     freq vector, the freq of the most common hyplotype is the last
   *D:        D=0
   *sigmas:   sigmas=1
   *beta:     all beta=start #(0.01)
   *a:        gamma prior of lamda a=20
   *b:        gamma prior of lamda b=20
   *asig:     inv-gamma prior of sigma squared asig=
   *bsig:     inv-gamma prior of sigma squared bsig=
   *uniq_x:   get all unique rows of design matrix, then do the transpose, then                   read it a vector(row1, row2,...)
   *tot_uniq_x: # of unique rows
   *lambda:   initial of lamda(prior of beta) =1
   *NUM_IT:   # of iteration = 50000
   *BURN_IN:  start burning iteration
   *beta_out: # of sampled beta's
   *lambda_out: # num of sampled lambda's(prior of beta)
   *freq_out: #num of sampled f's
   *D_out:    # number of 'd'
   
   *x_mat: design matrix, # of rows=total number of psedo subjects, # of cols=total number of haps except baseline
   *xhaplo[i][j][k]: transform x_mat to 3-dimension
   *xz: similar to x_mat, but with orignal # of subjects
   *h_mat: corresponds to freq_num, # of rows=total number of psedo subjects, 2 cols
#freq_map[i][j][k]: transform h_mat to 3-dimension
   *hz: similar to h_mat, but with orignal # of subjects
   
   uniq_h_mat or uniq_freq_num: nrow=number of unique hap pairs, ncol=2
   num_uniq_hap_pair: number of unique hap pairs
   */
{
  printf("Start MCMC sampling...\n");
  
  int i, j, k, l, m, n, p, ***freq_map, which_beta, **hz, h_mat[*tot_hap][2], **uniq_h_mat, **match_freq_map_to_uniq, *match_hz_to_uniq, it=0, it1=0, it2=0, heat_it;
  double **xz, x_mat[*tot_hap][*x_length], ***xhaplo, *A, *ele, *xbeta;
  GXE_y_new=calloc(*N, sizeof(double));
  xhaplo=calloc(*N, sizeof(double));
  freq_map=calloc(*N, sizeof(int*));
  uniq_h_mat=calloc(*num_uniq_hap_pair, sizeof(int*));
  match_freq_map_to_uniq=calloc(*N, sizeof(int*));
  A=calloc(*num_uniq_hap_pair, sizeof(double));
  ele=calloc(*N, sizeof(double));
  xbeta=calloc(*N, sizeof(double));
  match_hz_to_uniq=calloc(*N, sizeof(int));
  for(i=0;i<*N;++i){
    freq_map[i]=calloc(num_haplo_id[i], sizeof(int*));
    xhaplo[i]=calloc(num_haplo_id[i], sizeof(double));
    match_freq_map_to_uniq[i]=calloc(num_haplo_id[i], sizeof(int));
  }
  for(i=0;i<*N;++i){
    for(j=0;j<num_haplo_id[i];++j){
      xhaplo[i][j]=calloc(*x_length, sizeof(double));
      freq_map[i][j]=calloc(2, sizeof(int));
    }
  }
  for(i=0;i<*num_uniq_hap_pair;++i){
    uniq_h_mat[i]=calloc(2, sizeof(int));
  }
  /* separating x vector from R as x matrix */
  l=0;
  for(j=0;j<*x_length;++j){
    for(i=0;i<*tot_hap;++i){
      x_mat[i][j]=x[l];
      ++l;
    }
  }
  /* separating haplo map vector from R as h matrix */
  l=0;
  for(j=0;j<2;++j){
    for(i=0;i<*tot_hap;++i){
      h_mat[i][j]=freq_num[l];
      ++l;
    }
  }
  
  l=0;
  for(j=0;j<2;++j){
    for(i=0;i<*num_uniq_hap_pair;++i){
      uniq_h_mat[i][j]=uniq_freq_num[l];
      ++l;
    }
  }
  
  xz=calloc(*N, sizeof(double));
  hz=calloc(*N, sizeof(int*));
  for(i=0;i<*N;++i){
    xz[i]=calloc(*x_length, sizeof(double));
    hz[i]=calloc(2, sizeof(int));
  }
  // x_mat is the design matrix, xz is the matrix with orignal # of subject.
  /* separating haplotypes per person from the x matrix and assigning z (missing haplotype) for persons with only one compatible haplotype */
  l=0;
  for(i=0;i<*N;++i){
    GXE_y_new[i]=y[l];
    // n_cases += GXE_y_new[i];
    for(j=0;j<num_haplo_id[i];++j){
      m=0;
      for(k=0;k<*x_length;++k){
        xhaplo[i][j][k]=x_mat[l][m];
        if(num_haplo_id[i]==1) xz[i][k]=x_mat[l][m];
        else xz[i][k]=0;
        ++m;
      }
      ++l;
    }
  }
  l=0;
  for(i=0;i<*N;++i){
    for(j=0;j<num_haplo_id[i];++j){
      m=0;
      for(k=0;k<2;++k){
        freq_map[i][j][k]=h_mat[l][m];
        if(num_haplo_id[i]==1) hz[i][k]=h_mat[l][m];
        else hz[i][k]=0;
        ++m;
      }
      ++l;
      
      //initializly match all psedo persons' hap pairs to the storage
      for(p=0;p<*num_uniq_hap_pair;++p){
        if(uniq_h_mat[p][0]==freq_map[i][j][0] && uniq_h_mat[p][1]==freq_map[i][j][1]){
          //uniq_h_mat[p]==freq_map[i][j] does not work. not able to compare vectors
          match_freq_map_to_uniq[i][j]=p;
          //printf("Found p! p is %d\n", p);
          break;
        }
      }
      
    }
    
    if(num_haplo_id[i]==1){
      match_hz_to_uniq[i]=p;
    }else{
      match_hz_to_uniq[i]=0;
    }
    
  }
  //printf("Initial match_hz_to_uniq[0] is %d\n",match_hz_to_uniq[0]);
  //printf("Initial match_hz_to_uniq[1] is %d\n",match_hz_to_uniq[1]);
  //printf("Initial match_hz_to_uniq[2] is %d\n",match_hz_to_uniq[2]);
  //initialize A
  for(p=0;p<*num_uniq_hap_pair;++p){
    A[p]=GXE_calc_a(freq, uniq_h_mat[p], *D);
  }
  //initialize ele
  for(i=0;i<*N;++i){
    if(num_haplo_id[i]==1){
      xbeta[i]=beta[0];
      for(k=0;k<*x_length;++k){
        xbeta[i]=xbeta[i]+xz[i][k]*beta[k+1];
      }
      ele[i]=(GXE_y_new[i]-xbeta[i])*(GXE_y_new[i]-xbeta[i]);
    }else{ //xz=0
      xbeta[i]=0;
      ele[i]=0;
    }
  }
  //printf("Initial A[0] is %lf\n",A[0]); //printf *A, printf A[0]
  //printf("Initial A[1] is %lf\n",A[1]);
  //printf("Initial A[2] is %lf\n",A[2]);
  //printf("Initial A[3] is %lf\n",A[3]);
  //printf("Initial A[4] is %lf\n",A[4]);
  //printf("Initial ele[0] is %lf\n",ele[0]);
  //printf("Initial ele[1] is %lf\n",ele[1]);
  //printf("Initial ele[2] is %lf\n",ele[2]);
  //printf("Initial ele[3] is %lf\n",ele[3]);
  //printf("Initial ele[4] is %lf\n",ele[4]);
  
  /********************** start MCMC here ***************************/
  heat_it=0; /* tracks # of heating iterations */
  for(n=0;n<*NUM_IT;++n){
    heat_it=heat_it+1;
    if(heat_it==101) heat_it=1;
    int h=0;
    /* assigning or updating z (missing haplotype) for persons with more than one compatible haplotype */
    for(i=0;i<*N;++i){
      if(num_haplo_id[i]>1){
        GXE_update_z(i, freq_map[i], hz[i], beta, freq, *D, *sigmas, GXE_y_new[i], xz[i], xhaplo[i],
                     num_haplo_id[i], *x_length, *num_uniq_hap_pair, uniq_h_mat, A, match_freq_map_to_uniq[i], match_hz_to_uniq, ele, xbeta);
        
      }
      // find out updated xz to test the algorithm
      for(k=0;k<*x_length;++k){
        up_xz[h]=xz[i][k];
        // printf("%d th row, %d th column of xz is %lf\n",i+1,k+1,xz[i][k]);
        // printf("%d th up_xz is %lf\n",h+1,up_xz[h]);
        ++h;
      }
      //
    }
    /* update beta parameters */
    // printf("sigmas is  =%lf\n",*sigmas);
    for(i=0;i<*x_length+1;++i){
      which_beta=i;
      GXE_update_beta(*sigmas, xz, beta, which_beta, *lambda, freq, *D, *x_length, *h_length, GXE_y_new, *N,
                      *n_cov, ele, xbeta);
    }
    //printf("Updated beta[0] is %lf\n", beta[0]);
    //printf("After updating beta ele[0] is %lf\n",ele[0]);
    /* update sigma parameters */
    *sigmas=GXE_update_sigma(beta, *asig, *bsig, xz, *N, *x_length, ele);
    //printf("New sigmas is %lf\n",*sigmas);
    /* update lambda */
    *lambda=GXE_update_lambda(beta, *a, *b, *x_length);
    /* update frequencies and D */
    GXE_update_freq(freq, beta, *D, *h_length, xz, *N, hz, *n_cov, *CC, *num_uniq_hap_pair, uniq_h_mat, A, match_hz_to_uniq);
    
    /* update D parameter */
    //printf("Old D is: %lf\n", *D);
    //printf("Old A[0] is: %lf\n", A[0]);
    *D=GXE_update_D(freq, beta, *D, *h_length, xz, *N, hz, *n_cov, *num_uniq_hap_pair, uniq_h_mat, A, match_hz_to_uniq);
    //printf("After updating D A[0] is: %lf\n", A[0]);
    //printf("New D is: %lf\n", *D);
    if(n>=*BURN_IN){
      for(i=0;i<*x_length+1;++i){
        beta_out[it]=beta[i];
        ++it;
      }
      lambda_out[it2]=*lambda;
      sigmas_out[it2]=*sigmas;
      //printf("Updated sigma squared is %lf",*sigmas);
      for(i=0;i<*h_length;++i){
        freq_out[it1]=freq[i];
        ++it1;
      }
      D_out[it2]=*D;
      /* Calcuate BDIC*/
      double fit=0;
      for(i=0;i<*N;++i){
        fit=fit+ele[i]/(*sigmas*2);
      }
      BDIC_out[it2]=-*N*log(sqrt(2*M_PI*pow(*sigmas, 2)))-fit;
      /*printf("After xxbeta[0]=%lf\n",xxbeta[0]);
       printf("fit=%lf\n",fit);
       printf("N =%d\n",*N);
       printf("Pi=%lf\n",M_PI);
       printf("sigma=%lf\n",*sigmas);
       printf("Updated first BDIC_out=%lf\n",BDIC_out[it2]);
       */
      ++it2;
    }
    
  }
  //free memory
  for (i=0; i<*N; ++i)
  {
    free(hz[i]);
    free(xz[i]);
  }
  free(xz);
  free(hz);
  
  for (i =0; i<*N; ++i)
  {
    for (j = 0; j < num_haplo_id[i]; ++j)
    {
      free(freq_map[i][j]);
      free(xhaplo[i][j]);
    }
  }
  
  for (i =0; i<*N; ++i)
  {
    free(xhaplo[i]);
    free(freq_map[i] );
  }
  free(freq_map);
  free(xhaplo);
  free(GXE_y_new);
  
}
void GXE_update_z(int i, int **per_map, int *per_hz, double *beta, double *freq, double D, double sigmas, double per_y_new,
                  double *per_xz, double **per_xhaplo, int num_haplo_per, int x_length, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *per_match_freq_map_to_uniq, int *match_hz_to_uniq, double *ele, double *xbeta){
  int q, k;
  double prob[num_haplo_per], cum_prob[num_haplo_per], sum_prob, x, xbeta_compatible[num_haplo_per], ele_compatible[num_haplo_per];
  for(q=0;q<num_haplo_per;++q){
    prob[q]=A[per_match_freq_map_to_uniq[q]];
    xbeta_compatible[q]=beta[0];
    for(k=0;k<x_length;++k){
      xbeta_compatible[q]=xbeta_compatible[q]+per_xhaplo[q][k]*beta[k+1];
    }
    ele_compatible[q]=(per_y_new-xbeta_compatible[q])*(per_y_new-xbeta_compatible[q]);
    prob[q]=prob[q]*exp(-ele_compatible[q]/(2*sigmas));
  }
  sum_prob=GXE_sum(prob, num_haplo_per);
  for(q=0;q<num_haplo_per;++q){
    prob[q]=prob[q]/sum_prob;
  }
  cum_prob[0]=prob[0];
  for(q=1;q<num_haplo_per;++q){
    cum_prob[q]=cum_prob[q-1]+prob[q];
  }
  GetRNGstate();
  x=runif(0, 1);
  PutRNGstate();
  if(x<cum_prob[0]){ //accept the first hap pair in storage per_xhaplo
    //printf("Old match_hz_to_uniq is %d\n", per_match_hz_to_uniq);
    match_hz_to_uniq[i]=per_match_freq_map_to_uniq[0];
    xbeta[i]=xbeta_compatible[0];
    ele[i]=ele_compatible[0];
    for(k=0;k<x_length;++k){
      per_xz[k]=per_xhaplo[0][k];
    }
    for(k=0;k<2;++k){
      per_hz[k]=per_map[0][k];
    }
  }
  for(q=1;q<num_haplo_per;++q){
    if(x>cum_prob[q-1]&&x<cum_prob[q]){ //accept the (i+1)th hap pair in storage per_xhaplo
      match_hz_to_uniq[i]=per_match_freq_map_to_uniq[q];
      xbeta[i]=xbeta_compatible[q];
      ele[i]=ele_compatible[q];
      for(k=0;k<x_length;++k){
        per_xz[k]=per_xhaplo[q][k];
      }
      for(k=0;k<2;++k){
        per_hz[k]=per_map[q][k];
      }
    }
  }
}
void GXE_update_beta(double sigma, double **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int h_length, double *GXE_y_new, int N, int n_cov, double *ele, double *xbeta){
  double beta_new, g_1, g_2, f_old, f_new, beta_new_vec[x_length+1], SD, accept_prob, xbetan[N], ele_new[N];
  
  if(which_beta==0){
    double mu_0=0, std_0;
    int i;
    for(i=0;i<N;++i){
      mu_0=mu_0+GXE_y_new[i]-xbeta[i]+beta[0];
    }
    
    mu_0=(mu_0/sigma)/(1+N/sigma);  //prior for beta0 is N(0,1)
    std_0=sqrt(1/(1+N/sigma));
    
    GetRNGstate();
    beta_new=rnorm(mu_0, std_0);
    PutRNGstate();
    
    for(i=0;i<N;++i){
      xbeta[i]=xbeta[i]-beta[0]+beta_new; //y-sum(bnew,x0b1,x1b2..)
      ele[i]=(GXE_y_new[i]-xbeta[i])*(GXE_y_new[i]-xbeta[i]);
    }
    beta[0]=beta_new;
    
  }else{
    int i, x;
    beta_new=GXE_gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
    g_1=-lambda*fabs(beta[which_beta]);
    g_2=-lambda*fabs(beta_new);
    for(i=0;i<x_length+1;++i){
      beta_new_vec[i]=beta[i];
    }
    beta_new_vec[which_beta]=beta_new;
    for(i=0;i<N;++i){ xbetan[i]=xbeta[i]-xz[i][which_beta-1]*beta[which_beta]+xz[i][which_beta-1]*beta_new_vec[which_beta];
      ele_new[i]=(GXE_y_new[i]-xbetan[i])*(GXE_y_new[i]-xbetan[i]);
      g_1=g_1-ele[i]/(2*sigma);
      g_2=g_2-ele_new[i]/(2*sigma);
    }
    SD=sqrt(fabs(beta_new));
    f_old=exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
    SD=sqrt(fabs(beta[which_beta]));
    f_new=exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
    accept_prob=exp(g_2-g_1)*f_old/f_new;
    
    if(accept_prob>1){
      beta[which_beta]=beta_new;
      for(i=0;i<N;++i){
        xbeta[i]=xbetan[i];
        ele[i]=ele_new[i];
      }
    }else{
      GetRNGstate();
      x=rbinom(1, accept_prob);
      PutRNGstate();
      if(x==1){ //accept beta_new
        beta[which_beta]=beta_new;
        for(i=0;i<N;++i){
          xbeta[i]=xbetan[i];
          ele[i]=ele_new[i];
        }
        accept_count++;
      }
    }
  }
  
}

double GXE_update_lambda(double *beta, double a, double b, int x_length){
  double lambda, beta_abs[x_length];
  int i;
  for(i=1;i<x_length+1;++i)
    beta_abs[i-1]=fabs(beta[i]);
  GetRNGstate();
  lambda=rgamma((double)a+x_length, 1/(GXE_sum(beta_abs, x_length)+b));
  PutRNGstate();
  return lambda;
}
double GXE_update_sigma(double *beta, double asig, double bsig, double **xz, int N, int x_length, double *ele){
  double sigma;
  GetRNGstate();
  sigma=rgamma((double)asig+N/2, 1/(bsig+GXE_sum(ele, N)/2));
  PutRNGstate();
  sigma=1/sigma;
  return sigma;
}
void GXE_update_freq(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov, int Con, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *match_hz_to_uniq){
  int i, p, C=Con, update=0;
  double prop_freq_with_last[h_length], g_old=0, g_new=0, accept_prob=0, f_old, f_new, b_old[h_length],
                                                                                            b_new[h_length], min_f_old, min_f_new, A_new[num_uniq_hap_pair];
  GetRNGstate();
  for(i=0;i<h_length;++i){
    b_old[i]=freq[i]*C;
  }
  
  GXE_dirichlet(b_old, h_length, prop_freq_with_last);
  
  min_f_old=GXE_find_min(freq, h_length);
  /* check if the constraint -min fk/(1-min fk) < d is satisfied */
  min_f_new=GXE_find_min(prop_freq_with_last, h_length);
  if(-min_f_new/(1-min_f_new)<D){
    /* needed in acceptance prob. computation */
    for(i=0;i<h_length;++i){
      b_new[i]=prop_freq_with_last[i]*C;
      if(b_new[i]<=0){
        printf("Warning: prop_freq_with_last=%e, b_new = %lf\n", prop_freq_with_last[i], b_new[i]);
        for(i=0;i<h_length;++i){
          printf("%e ", prop_freq_with_last[i]);
        }
      }
    }
    
    for(i=0;i<num_uniq_hap_pair;++i){
      A_new[i]=GXE_calc_a(prop_freq_with_last, uniq_h_mat[i], D);
    }
    
    g_old=log(1-min_f_old);
    g_new=log(1-min_f_new);
    
    for(i=0;i<N;++i){
      g_old=g_old+log(A[match_hz_to_uniq[i]]);
      g_new=g_new+log(A_new[match_hz_to_uniq[i]]);
    }
    /* calculate f(f*|f^(t)) = f_new and f(f^(t)|f*) = f_old */
    f_old=lgammafn(C);
    f_new=f_old;
    for(i=0;i<h_length;++i){
      f_old+=(b_new[i]-1)*log(freq[i])-lgammafn(b_new[i]);
      f_new+=(b_old[i]-1)*log(prop_freq_with_last[i])-lgammafn(b_old[i]);
    }
    
    accept_prob=exp(g_new-g_old+f_old-f_new);
    //  printf("accept probability is %lf\n",accept_prob);
    if(-min_f_new/(1-min_f_new)>=D) accept_prob=0;
    if(accept_prob>1) update=1;
    else update=rbinom(1, accept_prob);
    if(update==1){
      for(i=0;i<h_length;++i){
        freq[i]=prop_freq_with_last[i];
      }
      //update A based on accepted freq
      for(p=0;p<num_uniq_hap_pair;++p){
        A[p]=A_new[p];
      }
    }
  }
  if(-min_f_old/(1-min_f_old)>D){
    printf("error in updating f: Min f_old and D constraint violated min_f_old=%f\n", min_f_old);
    for(i=0;i<h_length;++i)
      printf("%f\n", freq[i]);
    printf("%f\n", D);
    exit(1);
  }
  if(-min_f_new/(1-min_f_new)>D&&accept_prob>0){
    printf("error in updating f: Min f_new and D constraint violated ");
    exit(1);
  }
  PutRNGstate();
}
double GXE_update_D(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov, int num_uniq_hap_pair, int **uniq_h_mat, double *A, int *match_hz_to_uniq){
  int i, p, update=0;
  double prop_D, accept_prob, g_old=0, g_new=0, min_f, delta=0.05, lower, upper, f_old, f_new, A_new[num_uniq_hap_pair];
  GetRNGstate();
  min_f=GXE_find_min(freq, h_length);
  lower=D-delta;
  upper=D+delta;
  if(lower<-min_f/(1-min_f)) lower=-min_f/(1-min_f);
  if(upper>1) upper=1;
  prop_D=runif(lower, upper);
  for(i=0;i<num_uniq_hap_pair;++i){
    A_new[i]=GXE_calc_a(freq, uniq_h_mat[i], prop_D);
  }
  g_old=0;
  g_new=0;
  for(i=0;i<N;++i){
    g_old=g_old+log(A[match_hz_to_uniq[i]]);
    g_new=g_new+log(A_new[match_hz_to_uniq[i]]);
  }
  f_new=1/(upper-lower);
  lower=prop_D-delta;
  upper=prop_D+delta;
  if(lower<-min_f/(1-min_f)) lower=-min_f/(1-min_f);
  if(upper>1) upper=1;
  f_old=1/(upper-lower);
  accept_prob=exp(g_new-g_old)*f_old/f_new;
  if(-min_f/(1-min_f)>D||-min_f/(1-min_f)>prop_D){
    printf("error in updating D:  Min f and D constraint violated ");
    exit(1);
  }
  if(accept_prob>1) update=1;
  else update=rbinom(1, accept_prob);
  if(update==1){
    //update A based on accepted D
    for(p=0;p<num_uniq_hap_pair;++p){
      A[p]=A_new[p];
    }
    return prop_D;
  }else {
    return D;
  }
  PutRNGstate();
}

/* Calculate a(F) that is in the denominator of the likelihood*/
double GXE_calc_a(double *freq, int *per_freq, double D){
  int i, j, k;
  double a;
  i=per_freq[0];
  j=per_freq[1];
  if(i==j){
    k=1;
  }else{
    k=0;
  }
  a=k*D*freq[i-1]+(2-k)*(1-D)*freq[i-1]*freq[j-1];
  return a;
}
/* function to find sum of real numbers */
double GXE_sum(double *x, int n){
  double sum=0.0;
  int i;
  for(i=0;i<n;++i)
    sum=sum+x[i];
  return sum;
}
/* function to calculate min. of an array of numbers of length n */
double GXE_find_min(double *arr, int n){
  int i;
  double min=arr[0];
  for(i=1;i<n;++i){
    if(min>arr[i]) min=arr[i];
  }
  return min;
}
/* function to generate from double exponential distribution */
double GXE_gen_double_exp(double mean, double SD){
  double x, gen_exp;
  GetRNGstate();
  x=runif(0, 1);
  //    printf("uniform is %lf\n",x);
  gen_exp=rexp(SD/sqrt(2));
  //    printf("Exp is %lf\n",gen_exp);
  PutRNGstate();
  if(x>0.5) return gen_exp+mean;
  else return -gen_exp+mean;
}
/* function to generate from Dirichet distribution */
void GXE_dirichlet(double *param, int dim, double *gen_sample){
  int i, j=0;
  double gen_gamma[dim], sum_gamma;
  GetRNGstate();
  for(i=0;i<dim;++i){
    //if (param[i]<=0) //printf("Warning: param=%lf\n", param[i]);
    gen_gamma[i]=rgamma(param[i], 1);
    if(gen_gamma[i]<0.000001){
      //printf("WarningOld: gen_gamma=%lf\n", gen_gamma[i]);
      gen_gamma[i]=0.000001;
      //printf("WarningNew: gen_gamma=%lf\n", gen_gamma[i]);
    }
    j=j+1;
    //  if (gen_gamma[i]<=0) //printf("Warning: gen_gamma=%lf, param=%lf\n", gen_gamma[i], param[i]);
  }
  sum_gamma=GXE_sum(gen_gamma, dim);
  for(i=0;i<dim;++i){
    gen_sample[i]=gen_gamma[i]/sum_gamma;
  }
  PutRNGstate();
}
