/**********************************************************************
 *
 * tnorm_mle.h
 *
 * copyright (c) 2008, Wei Sun, UNC-CH
 *
 * last modified Aug 20, 2008
 * first written Aug 20, 2008
 *
 * C functions for the xCNV package
 *
 * MLE of truncated normal distribution
 *
 **********************************************************************/

void tnorm_mleR(double* T, double* M1, double* M2, int* dims, 
  double* mu, double* sigma, int* left, double* roundoff);

double g(double h, double p, double Vpn2);

void tnorm_mle(double T, double M1, double M2, double p, 
  double* mu, double* sigma, int left, double roundoff, int maxIt);
  
double zeroin(double h0, double p, double Vpn2, double step, 
  double roundoff, int maxIt);
