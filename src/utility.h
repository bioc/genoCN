/**********************************************************************
 * 
 * utility.h
 *
 * copyright (c) 2008, Wei Sun, UNC-CH
 *
 * last modified Jun 13, 2008
 * first written Jun 02, 2008
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * utility C functions for the xCNV package
 *
 * Hidden Markov Model
 *
 **********************************************************************/

void zRestrict(double *x, double lower, double upper);

void max(double *v, int nl, int nh, double *val, int *idx);

int is_infinite (const double x);

void logsumexp(double* v, int *RL, double *lse);

void reorg(double *v, double ***m, int nrow, int ncol);

void readfile(double** mat, char *str, int nrl, int nrh, int ncl, int nch);

void Rprint_v(double* v, int nrl, int nrh);

void Rprint_ve(double* v, int nrl, int nrh);

void Rprint_vi(int* v, int nrl, int nrh);

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch);

void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch);

void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch);

double dtnorm(double x, double mu, double sigma, double truncate, int left);

double logL(double *r, int n, double mean_r, double sd_r);

