double find_max(double *arr, int n);

double GXE_find_min(double *arr, int n);

double mean(double *x, int n);

int distance_ind(double distance, double *distance_array, int DIM_DISTANCE);

int mode(int number_dist, double *dist, double *freq_table, double *mode_x);

double median(int n, int number_dist, double *dist, double *freq_table);

double std_dev(double *x, int n);

double corr(double *x, double *y, int n);

double quantile(int n, int number_dist, double *dist, double *freq_table, double which_quant);

double log_factorial(int n);

double comb(int n, int a);

double min(double a, double b);

double sum(double *x, int n);

double qgaus(double (*func)(double), double a, double b);

double quad2d(double (*func)(double, double), double x1, double x2);

double quad3d(double (*func)(double, double, double), double x1, double x2);

double f1(double x);

double f2(double y);

double f3(double z);

double beta(double z, double w);

double gammln(double xx);

double BetaInequality(double x_a, double x_b, double y_a, double y_b);

double betai(double a, double b, double x);

double betacf(double a, double b, double x);

void nrerror(char error_text[]);

float golden(float ax, float bx, float cx, float (*f)(float), float tol,
             float *xmin);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
            float (*func)(float));

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

double trapzd(double (*func)(double), double a, double b, int n);

double qromb(double (*func)(double), double a, double b);

double mode_discrete(int number_dist, double *dist, double *freq_table, int *index);
