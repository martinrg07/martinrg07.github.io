
/*=======================================================================================*/
// Author: Ryan Martin
// C subroutine for computing PR estimates and marginal likelihood for finite mixtures.
// See corresponding paper for details of the method.
// These subroutines are to be called by "sasa.R" -- an R front-end function.
/*=======================================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>


double vec_sum(double *x, int len) {

  int i;
  double s = 0.0;
  for(i = 0; i < len; i++) s += x[i];
  return s;

}


void PR_finite(double *lik, int *nobs, int *ngrid, double *w, double *f, double *logmix, int *nperm, int *perm) {

  int i, j, k;
  double sum_post, L = 0.0, LL;
  double *post = (double *)R_alloc(ngrid[0], sizeof(double));
  double *F = (double *)R_alloc(ngrid[0], sizeof(double));
  double *ff = (double *)R_alloc(ngrid[0], sizeof(double));

  for(k = 0; k < nperm[0]; k++) {

    LL = logmix[0];
    for(j = 0; j < ngrid[0]; j++) ff[j] = f[j];
    for(i = 0; i < nobs[0]; i++) {

      for(j = 0; j < ngrid[0]; j++) post[j] = ff[j] * lik[perm[i + nobs[0] * k] + nobs[0] * j];

      sum_post = vec_sum(post, ngrid[0]);
      LL += log(sum_post);

      for(j = 0; j < ngrid[0]; j++) ff[j] += w[i] * (post[j] / sum_post  - ff[j]);

    }

    L += (LL - L) / (k + 1);
    for(j = 0; j < ngrid[0]; j++) F[j] += (ff[j] - F[j]) / (k + 1);

  }

  logmix[0] = -L;
  for(j = 0; j < ngrid[0]; j++) f[j] = F[j];

}



