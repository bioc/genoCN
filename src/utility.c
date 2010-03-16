/**********************************************************************
 *
 * utility.c
 *
 * copyright (c) 2008, Wei Sun, UNC-CH
 *
 * last modified Jun 13, 2008
 * first written Jun 02, 2008
 *
 *
 * utilities C functions for the xCNV package
 *
 * Hidden Markov Model
 *
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/**********************************************************************
 *
 * zRestrict(x, lower, upper)
 * zRestrict x into the interval [lower, upper]
 *
 **********************************************************************/
 
void zRestrict(double *x, double lower, double upper){
 	  if(*x < lower){
 	  	  *x = lower;
 	  }else if(*x > upper){
 	  	  *x = upper;
 	  }
}

/**********************************************************************
 *
 * max(v, nl, nh, val, idx)
 * in vector x[nl:nh], find the maximum value and the corresponding index
 *
 **********************************************************************/

void max(double *v, int nl, int nh, double *val, int *idx)
{
    int i;
    *idx = nl;
    *val = v[nl];

    for(i=nl+1; i<nh; i++){
        if(v[i] > *val){
            *val = v[i];
            *idx = i;
        }
    }
}

/**********************************************************************
 *
 * is_infinite(x)
 *
 * return +1 if x is positive infinity, -1 if x is negative infinity
 * and 0 otherwise
 *
 **********************************************************************/
int is_infinite (const double x){
    double y = x - x;
    int s = (y!=y);

    if(s && x >0)
        return +1;
    else if(s && x < 0)
        return -1;
    else
        return 0;
}

/**********************************************************************
 *
 * logsumexp
 *
 * log(sum(exp(v)))
 *
 **********************************************************************/
void logsumexp(double* v, int *RN, double *lse)
{
		int i, idx=0, N=RN[0];
		double res, val=0.0;

		if(N==0){
				*lse = -1.0/0.0;
		}else if(N==1){
				*lse = v[0];
		}else{
				max(v, 0, N, &val, &idx);
				if(val==1.0/0.0){
						error("positive infinite value in v\n");
				}
				res = 0;
				for(i=0; i<N; i++){
				  if(i==idx || v[i]==-1.0/0.0){
				    continue;
				  }
          res = res + exp(v[i] - val);
		    }

				*lse = val + log(1+res);
    }
}

/**********************************************************************
 *
 * reorg
 *
 * Reorganize a vector to a matrix of given size.
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg(double *v, double ***m, int nrow, int ncol)
{
    int i;

    *m = (double **)R_alloc(nrow, sizeof(double*));

    (*m)[0] = v;
    if(nrow>1){
        for(i=1; i<nrow; i++){
            (*m)[i] = (*m)[i-1] + ncol;
        }
    }
}

/**********************************************************************
 *
 * readfile
 *
 * read data into a matrix with given rows and columns
 *
 **********************************************************************/

void readfile(double** mat, char *str, int nrl, int nrh, int ncl, int nch) {
    FILE *file;
    int i,j;
    char temp[255];
    file = fopen(str,"r+t");
    for (i = nrl; i <= nrh; i++) {
        for (j = ncl; j <= nch; j++) {
            fscanf(file,"%s",temp);
            mat[i][j] = (double) atof(temp);
        }
    }
    fclose(file);
}

/**********************************************************************
 *
 * print_v
 *
 * print out a vector with different forms
 *
 **********************************************************************/
void Rprint_v(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%f\t", v[i]);
	}
	Rprintf ("%f\n", v[i]);
}

void Rprint_ve(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%.2e\t", v[i]);
	}
	Rprintf ("%.2e\n", v[i]);
}

void Rprint_vi(int* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%d\t", v[i]);
	}
	Rprintf ("%d\n", v[i]);
}

/**********************************************************************
 *
 * print_m
 *
 * print out a matrix with different forms
 *
 **********************************************************************/

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i < nrh; i++){
        for(j = ncl; j < nch; j++){
            Rprintf ("%.2e\t", m[i][j]);
        }
        Rprintf("\n");
	}
}

void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i < nrh; i++){
        for(j = ncl; j < nch; j++){
            Rprintf ("%i\t", m[i][j]);
        }
        Rprintf("\n");
	}
}

void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i < nrh; i++){
        for(j = ncl; j < nch; j++){
            Rprintf ("%.4f\t", m[i][j]);
        }
        Rprintf("\n");
	}
}


/**********************************************************************
 *
 * dtnorm
 *
 * density of truncated normal distribution
 *
 * some distribution functions from R API:
 *
 * dnorm(double x, double mu, double sigma, int given_log)
 * pnorm(double x, double mu, double sigma, int lower_tail, int give_log) 
 *
 **********************************************************************/
double dtnorm(double x, double mu, double sigma, double truncate, int left){
	double y = dnorm(x, mu, sigma, 0);
	if(left){
		y = y/pnorm(x, mu, sigma, 0, 0);
	}else{
		y = y/pnorm(x, mu, sigma, 1, 0);
	}
	return(y);
}

/**********************************************************************
 *
 * logL
 *
 * log Likelihood, assume normal density
 *
 **********************************************************************/

double logL(double *r, int n, double mean_r, double sd_r){
    int i;
    double l = 0.0;

    for(i=0; i<n; i++){
        l += dnorm(r[i], mean_r, sd_r, 1);
    }

    return(l);
}
