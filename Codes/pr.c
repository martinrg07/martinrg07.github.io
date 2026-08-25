
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

double int_simp(double *f, int ng, double h);


/* Newton's predictive recursion (PR) algorithm.  */

void PR(double *lik, int *nobs, int *ngrid, double *gridsize, double *w, double *f, double *logmix, int *nperm, int *perm) {

  int i, j, k, n=nobs[0], ng=ngrid[0], np=nperm[0];
  double *post, o, int_post;
  post = (double *)calloc(ng * np, sizeof(double));
  for(k = 0; k < np; k++) {

    logmix[k] = 0.0;

    for(i = 0; i < n; i++) {

      for(j = 0; j < ng; j++) {

        o = f[j + k*ng] * lik[perm[i + k*n] + j*n];
        if(o < 1e-100) post[j + k*ng] = 1e-100;
        else post[j + k*ng] = o;

      }
      int_post = int_simp(post + k * ng, ng, gridsize[0]);
      if(int_post < 1e-100) int_post = 1e-100;
      logmix[k] += log(int_post);
      for(j = 0; j < ngrid[0]; j++) {

        post[j + k*ng] /= int_post;
        f[j + k*ng] += w[i] * (post[j + k*ng] - f[j + k*ng]);

      }

    }

  }
  free(post);

}



/* Simpson's rule for one-dimensional integration */

double int_simp(double *f, int ng, double h) {

  double *u, out0=0, out;
  int i;
  u = (double *)calloc(ng, sizeof(double));
  u[0] = 1;
  u[ng-1] = 1;
  for(i = 1; i < ng-1; i++) {

    u[i] = 2 * pow(2, i % 2);
    out0 += u[i] * f[i];

  }
  out = h * (out0 + f[0] + f[ng-1]) / 3;
  free(u);
  return out;

}

