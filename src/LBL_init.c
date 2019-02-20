#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void cLBLmcmc(int *x1, int *x2, int *tot_hap1, int *tot_hap2, int *y1, int *y2, int *N1, int *N2, int *num_haplo_id1, int *num_haplo_id2, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x, int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out, int monitor);
extern void famLBLmcmc(int *x, int *tot_hap, int *y, int *N, int *num_haplo_id, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x,
                       int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out, int monitor);
extern void getWts(int *size, int *ID, double *wts, double *num_prob);
extern void LBLmcmc(int *x, int *tot_hap, int *y, int *N, int *num_haplo_id, int *x_length, double *freq, double *D, double *beta, double *a, double *b, int *uniq_x,
                    int *tot_uniq_x, double *lambda, int *NUM_IT, int *BURN_IN, double *beta_out, double *lambda_out, double *freq_out, double *D_out, int monitor);

static const R_CMethodDef CEntries[] = {
    {"cLBLmcmc",   (DL_FUNC) &cLBLmcmc,   25},
    {"famLBLmcmc", (DL_FUNC) &famLBLmcmc, 21},
    {"getWts",     (DL_FUNC) &getWts,      5},
    {"LBLmcmc",    (DL_FUNC) &LBLmcmc,    21},
    {NULL, NULL, 0}
};

void R_init_LBL(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

