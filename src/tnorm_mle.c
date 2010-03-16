#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "tnorm_mle.h"

/**********************************************************************
 *
 * tnorm_mleR
 *
 * MLE of truncated normal, with input from an R function,
 * accept input and then call the c function tnorm_mle
 *
 **********************************************************************/

void tnorm_mleR(double* T, double* M1, double* M2, int* dims, 
  double* mu, double* sigma, int* left, double* roundoff){
  int s = dims[0];
  int n = dims[1];
  int maxIt = dims[2];
  double p = ((double)s)/((double)n);
  tnorm_mle(*T, *M1, *M2, p, mu, sigma, *left, *roundoff, maxIt);
}

/**********************************************************************
 *
 * g
 *
 * the function of h = (T-mu)/sigma, MLE of h, denoted by hh
 * can be solved by g(h) = 0
 * g should be a monotone decreaing function
 *
 * some distribution functions from R API:
 *
 * dnorm(double x, double mu, double sigma, int given_log)
 * pnorm(double x, double mu, double sigma, int lower_tail, 
 *   int given_log)
 * qnorm(double p, double mu, double sigma, int lower_tail, int log_p)
 **********************************************************************/

double g(double h, double p, double Vpn2){
	double y;
	y  = (p/((1-p)*Vpn2))*((2-Vpn2)*h + 2*sqrt(h*h + Vpn2));
	y -= dnorm(h, 0.0, 1.0, 0)/pnorm(h, 0.0, 1.0, 0, 0);
	return(y);
}

/**********************************************************************
 *
 * tnorm_mle
 *
 * T is the truncation value
 * M1 is the first moment of observed data
 * M2 is the second moment of observed data
 * p is proportion of observed data among all the data
 * mu is mean
 * sigma is standard deviation
 * left = 1 if trcate the left tail
 * 
 * MLE of trucated normal
 **********************************************************************/

void tnorm_mle(double T, double M1, double M2, double p, 
  double* mu, double* sigma, int left, double roundoff, int maxIt){
  
  double Vpn2, h0, hh, step=0.01;

  Vpn2 = 4*(M2 - 2*T*M1 + T*T)/(T-M1)/(T-M1);

  /*
   * definitino of qnorm:
   * double qnorm(double p, double mu, double sigma, int lower_tail, 
   *              int log_p)
   *
   * h0 is the initial estimate of hh.
   * here p is proportion of observed data, so if we truncate the 
   * left tail, p is the upper (right) tail probability, therefore  
   * the parameter lower_tail should be 1-left
   */
  h0 = qnorm(p, 0.0, 1.0, 1-left, 0); 
  hh = zeroin(h0, p, Vpn2, step, roundoff, maxIt);
    
  if(left){
  	hh = -hh;
  	*sigma = 0.5*(T - M1)*(-hh - sqrt(hh*hh + Vpn2));
  }else{
  	*sigma = 0.5*(T - M1)*(-hh + sqrt(hh*hh + Vpn2));  	
  }
  
  *mu = T - *sigma*hh;
}

/**********************************************************************
 *
 * zeroin
 *
 * identify hh such that g(hh) = 0
 **********************************************************************/

double zeroin(double h0, double p, double Vpn2, double step, 
  double roundoff, int maxIt){
  	
  int i;
  double h1, h2, y0, y1, y2;
  	
	h1 = h0 - step;
  h2 = h0 + step;
  
  y0 = g(h0, p, Vpn2);
  y1 = g(h1, p, Vpn2);
  y2 = g(h2, p, Vpn2);
  
  /*
   * g should be a mononte decreasing function 
   * therefore y1 > y0 > y2
   */
  if(y0 > y1 || y2 > y0){
  	error("function g is not montone decreasing\n");
  }
  
  /* ------------------------------------------------ 
   * Suppose the solution is if hh, where g(hh) = 0, 
   * because of the monotone of function g, there 
   * are four possibilities:
   * (1) hh < h1 < h0 < h2, if(y1 < 0)
   * (2) h1 < hh < h0 < h2, if(y1 > 0 && y0 < 0)
   * (3) h1 < h0 < hh < h2, if(y0 > 0 && y2 < 0)
   * (4) h1 < h0 < h2 < hh, if(y2 > 0)
   *
   * adjust the h0, h2 so that h0 < hh < h2
   * ------------------------------------------------
   */
   	
  if(fabs(y0) < roundoff){
  	return(h0);
  }else if(fabs(y1) < roundoff){
  	return(h1);
  }else if(fabs(y2) < roundoff){
  	return(h2);
  }else if(y1 < 0){
		while(y1 < 0){
			h1 = h1 - step;
			y1 = g(h1, p, Vpn2);
		}
  	h0 = h1;
  	h2 = h0 + step;
  }else if(y1 > 0 && y0 < 0){
  	h2 = h0;
  	h0 = h1;
  }else if(y0 > 0 && y2 < 0){
  	// do nothing
  }else if(y2 > 0){
		while(y2 > 0){
			h2 = h2 + step;
			y2 = g(h2, p, Vpn2);
		}
  	h0 = h2 - step;
  }else{
  	error("hm, I do not think there is anything else");
  }
 
  y0 = g(h0, p, Vpn2);
  y2 = g(h2, p, Vpn2);
 
  for(i=0; i < maxIt; i++){
  	h1 = 0.5*(h0 + h2);
    y1 = g(h1, p, Vpn2);	
  	if(fabs(y1) < roundoff){
  		return(h1);
  	}else{
  		if(y1 > 0){
  			h0 = h1;
  			y0 = y1;
  		}else{
  			h2 = h1;
  			y2 = y1;
  		}
  	}
  }
  error("zeroin fail to converge, h0=%f, h1=%f, roundoff=%f, maxIt=%d", 
   h0, h1, roundoff, maxIt);
  return(h0);
}
