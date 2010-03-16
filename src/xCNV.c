/**********************************************************************
 *
 * xCNV.c
 *
 * copyright (c) 2008,2009, Wei Sun, UNC-CH
 *
 * last modified Oct 26, 2009
 * first written Jun 16, 2008
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/genoCN package
 *
 * Hidden Markov Model for CNV/CNA detection using SNP array data
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
#include "tnorm_mle.h"
#include "xCNV.h"

/**********************************************************************
 *   HMM: 9 states:
 *
 # -----------------------------------------------------------------
 #  HMM state  Genotype state
 # -----------------------------------------------------------------
 #  1          1 -> AA      2 -> AB        3 -> BB
 #  2          4 -> AA      5 -> (AA,AB)   6 -> (BB,AB)   7 -> BB
 #  3          8 -> NULL
 #  4          9 -> A      10 -> (A,AB)   11 -> (B,AB)   12 -> B
 #  5          13 -> AAA   14 -> AAB      15 -> ABB      16 -> BBB
 #  6/CNV      17 -> AAAA  18 -> AAAB     19 -> AABB     20 -> ABBB 21 -> BBBB
 #  6/CNA      22 -> AAA   23 -> (AAA,AB) 24 -> (BBB,AB) 25 -> BBB
 #  7          26 -> AAAA  27 -> AABB     28 -> BBBB
 #  8          29 -> AAAA  30 -> (AAAA,AB)31 -> (BBBB,AB)32 -> BBBB
 #  9          33 -> AAAA  34 -> AAAB     35 -> ABBB     36 -> BBBB
 # -----------------------------------------------------------------
 #
 **********************************************************************/


/**********************************************************************
 * function weights
 * 
 * update on Oct 26, 2009
 *
 * calculate weights for each genotype class within a HMM state
 * 
 * In previous versions of genoCN (v<=1.05), weights are pre-calculated 
 * using an R function, get.weights. However it requires too much memory. 
 * So I use this function to calculate the weights on the fly. 
 *
 * redundant computation will be done for each iteration of the 
 * Baulm-Welch (EM) algorithm
 *
 * NOTE, HERE z indicates stats, which START FROM 0, NOT FROM 1
 **********************************************************************/


int weights(double* ws, double pbf, int z, int CNA, int contamination, 
            int normalGtp, double geno_error){
  
  double paf = 1 - pbf;
  
  if (z==2) {
    ws[0]=1.0;
    return(1);
  }
  
  /****************************
   * for CNV
   ****************************/
  if(!CNA){
    if (z==0) {
      // state 1 (AA, AB, BB)
      ws[0] = paf*paf;
      ws[1] = 2*paf*pbf;
      ws[2] = pbf*pbf;
    }else if (z==1 || z==3) {
      // state 2 (AA, (AA,AB), (BB,AB), BB)
      // state 4 (A, (A,AB), (B,AB), B)        
      ws[0] = paf;
      ws[3] = pbf;
    }else if (z==4) {
      // state 5 (AAA, AAB, ABB, BBB)
      ws[0] = paf*paf*paf;
      ws[1] = 3.0*paf*paf*pbf;
      ws[2] = 3.0*paf*pbf*pbf;
      ws[3] = pbf*pbf*pbf;
    }else if (z==5) {
      // state 6 (AAAA, AAAB, AABB, ABBB, BBBB)
      ws[0] = paf*paf*paf*paf;
      ws[1] = 4.0*paf*paf*paf*pbf;
      ws[2] = 6.0*paf*paf*pbf*pbf;
      ws[3] = 4.0*paf*pbf*pbf*pbf;
      ws[4] = pbf*pbf*pbf*pbf;
    }
  }else {
  /****************************
   * for CNA
   ****************************/
    if (normalGtp == -1) {
    /*******************************
     * if normal genotype is missing 
     *******************************/
      
      if (z%2 == 1) {
        // state 2 (AA, (AA,AB), (BB,AB), BB)
        // state 4 (A, (A,AB), (B,AB), B)
        // state 6 (AAA, (AAA, AB), (AB, BBB), BBB)
        // state 8 (AAAA, (AAAA,AB), (BBBB,AB), BBBB)
        if (contamination) {
          ws[0] = paf*paf;
          ws[1] = ws[2] = paf*pbf;
          ws[3] = pbf*pbf;
        }else {
          ws[0] = paf;
          ws[3] = pbf;
        }
      }else if (z==0 || z==6) {
        // state 1 (AA, AB, BB)
        // state 7 (AAAA, AABB, BBBB)
        ws[0] = paf*paf;
        ws[1] = 2*paf*pbf;
        ws[2] = pbf*pbf;
      }else if (z==4 || z==8) {
        // state 5 (AAA, AAB, ABB, BBB)
        // state 9 (AAAA, AAAB, ABBB, BBBB)
        ws[0] = paf*paf;
        ws[1] = ws[2] = paf*pbf;
        ws[3] = pbf*pbf;
      }    
    }else if (normalGtp == 0) {
      /*******************************
       * if normal genotype is AA 
       *******************************/

      ws[0] = 1 - geno_error;
      
      if (z==0 || z==6) {
        // state 1 (AA, AB, BB)
        // state 7 (AAAA, AABB, BBBB)
        ws[1] = ws[2] = geno_error/2;
      }else if (z%2 == 1) {
        // state 2 (AA, (AA,AB), (BB,AB), BB)
        // state 4 (A, (A,AB), (B,AB), B)
        // state 6 (AAA, (AAA, AB), (AB, BBB), BBB)
        // state 8 (AAAA, (AAAA,AB), (BBBB,AB), BBBB)
        if (contamination) {
          ws[1] = ws[2] = ws[3] = geno_error/3;
        }else {
          ws[3] = geno_error;
        }
      }else if (z==4 || z==8) {
        // state 5 (AAA, AAB, ABB, BBB)
        // state 9 (AAAA, AAAB, ABBB, BBBB)
        ws[1] = ws[2] = ws[3] = geno_error/3;
      }    
      
    }else if (normalGtp == 1) {
      /*******************************
       * if normal genotype is AB 
       *******************************/
            
      if (z==0 || z==6) {
        // state 1 (AA, AB, BB)
        // state 7 (AAAA, AABB, BBBB)
        ws[1] = 1 - geno_error;
        ws[0] = ws[2] = geno_error/2;
      }else if (z%2 == 1) {
        // state 2 (AA, (AA,AB), (BB,AB), BB)
        // state 4 (A, (A,AB), (B,AB), B)
        // state 6 (AAA, (AAA, AB), (BBB, AB), BBB)
        // state 8 (AAAA, (AAAA,AB), (BBBB,AB), BBBB)
        if (contamination) {
          ws[1] = ws[2] = (1 - geno_error)/2;
          ws[0] = ws[3] = geno_error/2;
        }else {
          ws[0] = ws[3] = 0.5;
        }
      }else if (z==4 || z==8) {
        // state 5 (AAA, AAB, ABB, BBB)
        // state 9 (AAAA, AAAB, ABBB, BBBB)
        ws[1] = ws[2] = (1 - geno_error)/2;
        ws[0] = ws[3] = geno_error/2;
      }    
      
    }else if (normalGtp == 2) {
      /*******************************
       * if normal genotype is BB 
       *******************************/
      
      if (z==0 || z==6) {
        // state 1 (AA, AB, BB)
        // state 7 (AAAA, AABB, BBBB)
        ws[2] = 1 - geno_error;
        ws[0] = ws[1] = geno_error/2;
      }else if (z%2 == 1) {
        // state 2 (AA, (AA,AB), (BB,AB), BB)
        // state 4 (A, (A,AB), (B,AB), B)
        // state 6 (AAA, (AAA, AB), (BBB, AB), BBB)
        // state 8 (AAAA, (AAAA,AB), (BBBB,AB), BBBB)
        ws[3] = 1 - geno_error;
        if (contamination) {
          ws[0] = ws[1] = ws[2] = geno_error/3;
        }else {
          ws[0] = geno_error;
        }
      }else if (z==4 || z==8) {
        // state 5 (AAA, AAB, ABB, BBB)
        // state 9 (AAAA, AAAB, ABBB, BBBB)
        ws[3] = 1 - geno_error;
        ws[0] = ws[1] = ws[2] = geno_error/3;
      }    
      
    }
  }
  return(1);
}

/**********************************************************************
 * function emiss4
 *
 * emission probabiity for one observation of BAF
 * for those states with homozygous genotypes and tissue contamination
 *
 **********************************************************************/

double emiss4(double b, double EPSILON, double* ws, double* mu_b, 
 double* sd_b, double pi_b, int contam){
 
  double eB1, eBaf=0.0, pn;
  
  if(b>EPSILON && b<1.0-EPSILON){
    eB1 = 0.0;
    eB1 += ws[0]*dnorm(b,mu_b[0],sd_b[0],0);
    eB1 += ws[3]*dnorm(b,mu_b[3],sd_b[3],0);
    if(contam){
      eB1 += ws[1]*dnorm(b,mu_b[1],sd_b[1],0);
      eB1 += ws[2]*dnorm(b,mu_b[2],sd_b[2],0);
    }        
    eBaf = pi_b + (1-pi_b)*eB1;
  }else if(b<EPSILON){
    /* pnorm(double x, double mu, double sigma,  
     *       int lower_tail, int give_log)
     */
    if(contam){   
      pn = pnorm(0.0, mu_b[1], sd_b[1], 1, 0);
      eBaf = (1-pi_b)*(0.5*ws[0] + ws[1]*pn);
    }else{
      eBaf = (1-pi_b)*0.5*ws[0];
    }
  }else if(b>1.0-EPSILON){
    if(contam){
      pn = pnorm(1.0, mu_b[2], sd_b[2], 0, 0);
      eBaf = (1-pi_b)*(0.5*ws[3] + ws[2]*pn);
    }else{
      eBaf = (1-pi_b)*0.5*ws[3];
    }
  }
  return(eBaf);
}

/**********************************************************************
 * function emissK
 *
 * emission probabiity for one observation of BAF
 * for those states with at least one heterzygous genotype class
 *
 **********************************************************************/

double emissK(double b, double EPSILON, double* ws, double* mu_b, 
 double* sd_b, double pi_b, int Hz){
 
  double eB1, eBaf=0.0;
  int h;
  
  if(b>EPSILON && b<1.0-EPSILON){
    eB1 = 0.0;
    for(h=0; h<Hz; h++){
      eB1 += ws[h]*dnorm(b,mu_b[h],sd_b[h],0);
    }
    eBaf = pi_b + (1-pi_b)*eB1;
  }else if(b<EPSILON){
    eBaf = 0.5*(1-pi_b)*ws[0];
  }else if(b>1.0-EPSILON){
    eBaf = 0.5*(1-pi_b)*ws[Hz-1];
  }
  
  return(eBaf);
}

/**********************************************************************
 * function emiss
 *
 * emission probabiity for one chromosome
 *
 * function call:  "reorg@utility.c"
 *
 **********************************************************************/
void emiss(double* LRR, double* BAF, double* pBs, int* normalGtp,
           double* geno_error, double* pi_r, double* mu_r, double *sd_r, 
           double *pi_rcn, double *mu_rcn, double *sd_rcn, double *R_mR, 
           double *pi_b, double *mu_bR, double *sd_bR, int *dims, 
           int* is_snp_probe, int* CNAR, int* contamR, double *eLRRR, 
           double *eBAFR, double *emR)
{
  int L, N, M, K, i, z, k, nGtp, contam=*contamR, CNA=*CNAR;
  double b, r, pbf;
  double **eLRR, **eBAF, **em, R_m=*R_mR, **mu_b, **sd_b;
  // number of genotype classes for each state
  int H[9] = {3, 4, 1, 4, 4, 4, 3, 4, 4};
  double ws[10];
  double EPSILON = 1e-13; // rounding error

  if(!CNA){ 
    contam = 0; H[5] = 5; H[6] = H[7] = H[8] = 0; 
  }

  L = dims[0]; // number of observations
  N = dims[1]; // number of states
  M = dims[2]; // maximum number of genotypes per state
  
  reorg(eLRRR, &eLRR, L, N);
  reorg(eBAFR, &eBAF, L, N);
  reorg(emR,   &em,   L, N);
  reorg(mu_bR, &mu_b, N, M);
  reorg(sd_bR, &sd_b, N, M);
  
/*
  Rprintf("\nmu_b\n");
  Rprint_mf(mu_b, 0, N, 0, M);
  Rprintf("\nsd_b\n");
  Rprint_mf(sd_b, 0, N, 0, M);
*/
  for(i=0; i<L; i++){
    r = LRR[i];
    b = BAF[i]; 
    pbf  = pBs[i];
    nGtp = normalGtp[i];
    
    if(is_snp_probe[i]){
      /* emission prob for LRR as well as overall emission prob */
      for(z=0; z<N; z++){
        eLRR[i][z] = pi_r[z]/R_m+(1-pi_r[z])*dnorm(r,mu_r[z],sd_r[z],0);        
        em[i][z]   = log(eLRR[i][z]);
        
        weights(ws, pbf, z, CNA, contam, nGtp, *geno_error);
        
        if (z==5) {
          if (CNA) {
            eBAF[i][z] = emiss4(b, EPSILON, ws, mu_b[z], sd_b[z], pi_b[z], contam);
          }else {
            eBAF[i][z] = emissK(b, EPSILON, ws, mu_b[z], sd_b[z], pi_b[z], H[z]);
          }
        }else {
          if (z%2==1) {
            eBAF[i][z] = emiss4(b, EPSILON, ws, mu_b[z], sd_b[z], pi_b[z], contam);
          }else {
            eBAF[i][z] = emissK(b, EPSILON, ws, mu_b[z], sd_b[z], pi_b[z], H[z]);
          }
        }
        
        em[i][z] = em[i][z] + log(eBAF[i][z]);
      }
      
    }else{
      for(z=0; z<N; z++){
        eLRR[i][z] = pi_rcn[z]/R_m+(1-pi_rcn[z])*dnorm(r,mu_rcn[z],sd_rcn[z],0);        
        em[i][z]   = log(eLRR[i][z]);
      }
    }
  }    
}

/**********************************************************************
 *
 * transition
 *
 * evaluate transition matrix given position specific information
 *
 * function call: reorg@utility.c
 *
 **********************************************************************/

void transition_c(double** trans_m, double di, double* Ds, int N,
                double** trans1, double* transB, double distThreshold)
{
  double tmp;
  int j, k;
  
  if (di > distThreshold) {
    for(j=0; j<N; j++){
      for(k=0; k<N; k++){
        trans1[j][k] = transB[k];
      }
    }
  }else{
    for(j=0; j<N; j++){
      tmp = 1-exp(-di/Ds[j]);
      for(k=0; k<N; k++){
          trans1[j][k] = trans_m[j][k]*tmp;
      }
      trans1[j][j] = 1 - tmp;
    }
  }
}

/**********************************************************************
 *
 * viterbi
 *
 * find the best path
 *
 * function call: reorg@utility.c, "transition", max@utility.c
 *
 **********************************************************************/
void viterbi(double* pos, double* emR, int * em_inf, 
             double* trans_mR, double* trans_begin, 
             double* Ds, double* distThreshold, double* transB,
             int* dims, int* path, double *logP)
{
    int L, N, em_noinf, i, j, z, **path_m, idx=0, DEBUG=0;
    double **em, *v_old, *v_new, *tmp, di,
           **trans_m, **trans1, val=0.0;

    L = dims[0]; // number of observations
    N = dims[1]; // number of states
    em_noinf = dims[2]; // number of -Inf in em
    
    // Rprintf("L=%d, N=%d\n", L, N);

    for(i=0; i<em_noinf; i++){
        emR[em_inf[i]-1] = -1.0/0.0;
    }

    reorg(emR, &em, L, N);
    reorg(trans_mR, &trans_m, N, N);

    trans1 = (double**) Calloc(N, double*);
    trans1[0] = (double*) Calloc(N*N, double);
    for(i=1; i<N; i++){
        trans1[i] = trans1[0] + i*N;
    }

    v_old = (double *)Calloc(N, double);
    v_new = (double *)Calloc(N, double);

    path_m = (int**) Calloc(L, int*);
    path_m[0] = (int*) Calloc(L*N, int);
    for(i=1; i<L; i++){
        path_m[i] = path_m[0] + i*N;
    }

    tmp = (double *)Calloc(N, double);

    // initialization
    for(j=0; j<N; j++){
        v_old[j] = log(trans_begin[j]) + em[0][j];
    }
    
    if(DEBUG){
      Rprintf("trans_begin\n");
      for(j=0; j<N; j++){//transition from * to z
        Rprintf("%.2e ", trans_begin[j]);
      }
      Rprintf("\n");

      Rprintf("trans_m\n");
      for(j=0; j<N; j++){//transition from * to z
        for(z=0; z<N; z++){
          Rprintf("%.2e ", trans_m[j][z]);
        }
        Rprintf("\n");
      }
    }
  
    // recursion
    for(i=1; i<L; i++){
        // Rprintf("i=%d\n",i);
        di = pos[i] - pos[i-1];
        transition_c(trans_m, di, Ds, N, trans1, transB, *distThreshold);

        for(z=0; z<N; z++){//transition from * to z
            // debug, for transition to the LOH state
            /*
            if(z==1){
              for(j=0; j<N; j++){
                if(trans1[j][z] > 1e-16){
                  Rprintf("j=%d, z=%d, trans1[j][z]=%.2e\n", j, z, trans1[j][z]);
                }
              }
            }
            */
            for(j=0; j<N; j++){
                tmp[j] = v_old[j] + log(trans1[j][z]);
            }
            max(tmp, 0, N, &val, &idx);
            v_new[z] = val + em[i][z];
            path_m[i-1][z] = idx+1;
        }

        for(j=0; j<N; j++){
            v_old[j] = v_new[j];
        }
    }

    // termination
    max(v_old, 0, N, &val, &idx);
    *logP = val;
    path_m[L-1][0] = idx+1;


    // extract the final states
    path[L-1] = path_m[L-1][0];
    for(i=L-2; i>=0; i--){
        path[i] = path_m[i][path[i+1]-1];
    }

    // Rprint_vi(path, 0, L);
    // Rprint_mi(path_m, 0, L, 0, N);

    Free(v_old);
    Free(v_new);
    Free(tmp);

    Free(trans1[0]);
    Free(trans1);

    Free(path_m[0]);
    Free(path_m);
}


/**********************************************************************
 *
 * forward
 *
 * forward probability
 *
 * L    total number of observations
 * N    total number of states
 *
 * f    forward probability
 * f[i][z] = Prob(x_(1),..., x_(i), q_i=z | Theta)
 * i: # of observation, z: # of state
 *
 * function call: "reorg@utility.c", "transition", "logsumexp"
 *
 **********************************************************************/

void forward(double* pos, double* emR, int* em_inf, int* Rem_noinf,
             double* trans_mR, double* trans_begin, double* Ds,
             double* fR, int* dims, double* distThreshold, 
             double* transB)
{
    int i, j, z, L, N;
    double di, **f, **em, *v, **trans_m, **trans1, lse=0.0;

    for(i=0; i<*Rem_noinf; i++){
        emR[em_inf[i]-1] = -1.0/0.0;
    }

    L = dims[0]; // number of observations
    N = dims[1]; // number of states
    reorg(fR, &f, L, N);
    reorg(emR, &em, L, N);
    reorg(trans_mR, &trans_m, N, N);

    v = (double *)Calloc(N, double);

    trans1 = (double**) Calloc(N, double*);
    trans1[0] = (double*) Calloc(N*N, double);
    for(i=1; i<N; i++){
        trans1[i] = trans1[0] + i*N;
    }

    // Rprintf("L=%d, N=%d\n", L, N);

    //initialization
    // log(f[1,]) = log(trans.begin) + log(em[1,])
    for(j=0; j<N; j++){
      f[0][j] = log(trans_begin[j]) + em[0][j];
    }
    /*
    Rprintf("trans_begin[6]=%f\n", trans_begin[5]);
    Rprintf("em[1,6]=%f\n", em[0][5]);
    Rprintf("f[1,6]=%f\n", f[0][5]);
    */

    //recursion
    for(i=1; i<L; i++){
        di = pos[i] - pos[i-1];
        transition_c(trans_m, di, Ds, N, trans1, transB, *distThreshold);
        for(z=0; z<N; z++){// transition from j to z
            for(j=0; j<N; j++){
                v[j] = f[i-1][j] + log(trans1[j][z]);
            }
            logsumexp(v, &N, &lse);
            f[i][z] = em[i][z] + lse;
        }
    }

    Free(v);
    Free(trans1[0]);
    Free(trans1);
}

/**********************************************************************
 *
 * backward
 *
 * backward probability
 *
 * L    total number of observations
 * N    total number of states
 *
 * b    backward probability
 * b[i][z] = Prob(x_(i+1),..., x_(L)| q_i=z, Theta)
 * i: # of observation, z: # of state
 *
 * function call: "reorg@utility.c", "transition", "logsumexp"
 *
 **********************************************************************/

void backward(double* pos, double* emR, int* em_inf, int* Rem_noinf,
              double* trans_mR, double* Ds, double* bR, int* dims, 
              double* distThreshold, double* transB)
{
    int i, j, z, L, N;
    double di, **b, **em, **trans_m, *v, **trans1, lse=0.0;

    for(i=0; i<*Rem_noinf; i++){
        emR[em_inf[i]-1] = -1.0/0.0;
    }

    L = dims[0]; // number of observations
    N = dims[1]; // number of states
    reorg(bR, &b, L, N);
    reorg(emR, &em, L, N);
    reorg(trans_mR, &trans_m, N, N);

    v = (double *)Calloc(N, double);

    trans1 = (double**) Calloc(N, double*);
    trans1[0] = (double*) Calloc(N*N, double);
    for(i=1; i<N; i++){
        trans1[i] = trans1[0] + i*N;
    }

    /*
    Rprintf("L=%d, N=%d\n", L, N);
    */

    //initialization
    // b[L,z]=1, log(b[L,z])=0
    for(j=0; j<N; j++){
      b[L-1][j] = 0;
    }

    //recursion
    for(i=L-1; i>0; i--){
        di = pos[i] - pos[i-1];
        transition_c(trans_m, di, Ds, N, trans1, transB, *distThreshold);
        for(z=0; z<N; z++){// transition from z to j
            for(j=0; j<N; j++){
                v[j] = log(trans1[z][j]) + em[i][j] + b[i][j];
            }
            logsumexp(v, &N, &lse);
            b[i-1][z] = lse;
      }
    }

  Free(v);
  Free(trans1[0]);
  Free(trans1);
}

/**********************************************************************
 * posterior probability
 * pP = P(q_{i}=z|X,Theta)
 *
 * L is number of observations
 * N is number of states
 *
 * Function call: "logsumexp"
 *
 **********************************************************************/

void postP(double* fR, double* bR, int* f_inf, int* b_inf, int* no_inf,
           int* dims, double* pPR)
{
    double **f, **b, logP, **pP, lse=0.0;
    int L, N, i, j, f_noinf = no_inf[0], b_noinf=no_inf[1];

    for(i=0; i<f_noinf; i++){
        fR[f_inf[i]-1] = -1.0/0.0;
    }
    for(i=0; i<b_noinf; i++){
        bR[b_inf[i]-1] = -1.0/0.0;
    }

    L = dims[0];
    N = dims[1];
    reorg(fR, &f, L, N);
    reorg(bR, &b, L, N);
    reorg(pPR, &pP, L, N);

    // Rprintf("L=%d, N=%d\n", L, N);

    logsumexp(f[L-1], &N, &lse);
    logP = lse;

    for(i=0; i<L; i++){
        for(j=0; j<N; j++){
            pP[i][j] = exp(f[i][j] + b[i][j] - logP);
        }
    }
}

/**********************************************************************
 *
 * EM algorithm: update parameters
 *
 * Function call: "reorg@utility.c", "logsumexp@utility.c"
 *
 **********************************************************************/

void baum_welch(double* pos, double* LRR, double* BAF,
    double* pBs, int* normalGtp, double* geno_error, 
    int* is_snp_probe, int* RCNA, double* Ds, 
    double* RR_m, double* pi_r, double* pi_rcn, 
    double* mu_r, double* mu_r_upper, double* mu_r_lower,
    double* sd_r, double* sd_r_upper, double* sd_r_lower,
    double* mu_rcn, double* mu_rcn_upper, double* mu_rcn_lower,
    double* sd_rcn, double* sd_rcn_upper, double* sd_rcn_lower,
    double* pi_b, double* trans_mR, double* min_tp, 
    double* max_diffR, double* mu_bR, double* mu_bR_upper, 
    double* mu_bR_lower, double* sd_bR, 
    double* sd_bR_upper, double* sd_bR_lower,
    double* emR, int* em_inf, double* eLRRR,
    double* eBAFR, double* fR, double* bR, 
    int* f_inf, int* b_inf, int* no_inf, 
    double* pPR, int* contamR, 
    int* dims, int* len, double* distThreshold, 
    double* pi_r_new, double* mu_r_new, 
    double* sd_r_new, double* pi_rcn_new, 
    double* mu_rcn_new, double* sd_rcn_new, 
    double* pi_b_new, 
    double* mu_b_newR, double* sd_b_newR,
    double* trans_m_newR, double* trans_begin_newR)
{
  int L, LL, N, NN, M, MM, S, LS, CNA = *RCNA,
      i, j, k, k1, k2, k3, x, z, hh, contam=*contamR, 
      flag, zStart, zEnd, *Omiga, Oi, em_noinf = no_inf[0],
      f_noinf = no_inf[1], b_noinf = no_inf[2], DEBUG, 
      left, maxIt = 200, ri, ci, estimate_pi_r, estimate_pi_b, 
      estimate_trans_m;
  int k1s[4], rid[4], cid[4], nGtp;
  double ws[10], pbf, baf;
  double sum0, sum1, pnz0, pnz1;
  
  double **mu_b, **mu_b_upper, **mu_b_lower, 
         **sd_b, **sd_b_upper, **sd_b_lower, 
         **mu_b_new, **sd_b_new,           
         **em, **eLRR, **eBAF, **f, **b, **pP, transum0,
         transum1, **trans_m,  **trans_m_new, **trans_begin_new,
         *dis, *ps, *vjk, *vjk0, *lvjk, *lvjk_j,
         **rU, **rN, *rUS, *pPS, *rUS_cn, *pPS_cn, *rNcb, 
         sumz, muz, sdz, sumz_cn, muz_cn, sdz_cn, tmp,
         **bU, **bN, *bUS, *BAF01, d0, d1, dajust, 
         *bNcb, max_diff = *max_diffR, pi_b_upp=0.05,
         R_m = *RR_m, lse=0.0, lsum=0.0, pi_r_upp; 

  // varaibles for tnorm_mle: truncated normal MLE
  double M1, M2, p=0.0, T, roundoff=1e-6;

  double mu_b_tmp, sd_b_tmp;
  double EPSILON = 1e-13; // rounding error
  // minimum probability of one state, if the cumulative 
  // posterior probability of one state is smaller than 
  // minStateP, stop updating it
  double minStateP = 3.0; 

  // number of genotype classes for each state
  int H[9] = {3, 4, 1, 4, 4, 4, 3, 4, 4};
  if(!CNA){
    H[5] = 5; H[6] = H[7] = H[8] = 0; 
  }

  // starting index of genotype class for each state
  int zIdxStart[9] = {0, 3, 7, 8, 12, 21, 25, 28, 32};
  if(!CNA){
    zIdxStart[5] = 16; zIdxStart[6] = zIdxStart[7] = zIdxStart[8] = 0; 
  }

  for(i=0; i<em_noinf; i++){
    emR[em_inf[i]-1] = -1.0/0.0;
  }

  for(i=0; i<f_noinf; i++){
    fR[f_inf[i]-1] = -1.0/0.0;
  }

  for(i=0; i<b_noinf; i++){
    bR[b_inf[i]-1] = -1.0/0.0;
  }

  if(CNA){ 
    pi_r_upp=0.05; 
  }else{
    pi_r_upp=0.10; 
  }
  
  L  = dims[0]; /* (total) sequence length  */
  LL = L-1;
  N  = dims[1]; /* number of states         */
  NN = N-1;
  /* number of columns of mu.b, i.e., number of normal components
     in the mixture distribution.           */
  M  = dims[2]; 
  S  = dims[3]; /* number of sequence segments, e.g., chromosomes */
  LS = L-S;
  DEBUG = dims[4];
  estimate_pi_r = dims[5];
  estimate_pi_b = dims[6];
  estimate_trans_m = dims[7];
  
  // the total number of genotype classes
  if (CNA) { MM = 36; } else { MM = 21; }
  
   // Rprintf("L=%d, N=%d, M=%d, S=%d \n", L, N, M, S);
  reorg(emR, &em, L, N);
  reorg(eLRRR, &eLRR, L, N);
  reorg(eBAFR, &eBAF, L, N);
  reorg(fR, &f, L, N);
  reorg(bR, &b, L, N);
  reorg(pPR, &pP, L, N);
  reorg(mu_bR, &mu_b, N, M);
  reorg(mu_bR_upper, &mu_b_upper, N, M);
  reorg(mu_bR_lower, &mu_b_lower, N, M);
  reorg(sd_bR, &sd_b, N, M);
  reorg(sd_bR_upper, &sd_b_upper, N, M);
  reorg(sd_bR_lower, &sd_b_lower, N, M);
  reorg(trans_mR, &trans_m, N, N);
  reorg(mu_b_newR, &mu_b_new, N, M);
  reorg(sd_b_newR, &sd_b_new, N, M);
  reorg(trans_m_newR, &trans_m_new, N, N);
  reorg(trans_begin_newR, &trans_begin_new, S, N);
  
  if(DEBUG){
    Rprintf("mu_b:\n");
    Rprint_mf(mu_b,0,N,0,M);
    
    Rprintf("mu_b_upper:\n");
    Rprint_mf(mu_b_upper,0,N,0,M);

    Rprintf("mu_b_lower:\n");
    Rprint_mf(mu_b_lower,0,N,0,M);

    Rprintf("sd_b\n");
    Rprint_mf(sd_b,0,N,0,M);
    
    Rprintf("em\n");
    Rprint_mf(em,0,10,0,N);
  }

  /* variables for transition probability */
  dis    = (double *)Calloc(L-1, double);
  vjk    = (double *)Calloc(L-1, double);
  vjk0   = (double *)Calloc(LS, double);
  lvjk   = (double *)Calloc(N, double);
  lvjk_j = (double *)Calloc(N-1, double);
  ps     = (double *)Calloc(L, double);

  /* variables for emission probability of LRR */
  rUS  = (double *)Calloc(N, double);
  pPS  = (double *)Calloc(N, double);
  rNcb = (double *)Calloc(L, double);

  /* variables for emission probability of LRR of CN-only probes */
  rUS_cn  = (double *)Calloc(N, double);
  pPS_cn  = (double *)Calloc(N, double);

  /* variables for emission probability of BAF */
  bUS   = (double *)Calloc(N, double);
  Omiga = (int *)Calloc(L, int);
  BAF01 = (double *)Calloc(L, double);
  bNcb  = (double *)Calloc(L, double);

  /* posterior probability of r is from uniform distribution
   * each row corresponds to one probe
   * each column corresponds to one state
   */
  rU  = (double**) Calloc(L, double*);
  rU[0] = (double*) Calloc(L*N, double);
  for(i=1; i<L; i++){
    rU[i] = rU[0] + i*N;
  }

  /* posterior probability of r is from normal distribution
   * each row corresponds to one probe
   * ecah column corresponds to one state 
   */
  rN  = (double**) Calloc(L, double*);
  rN[0] = (double*) Calloc(L*N, double);
  for(i=1; i<L; i++){
    rN[i] = rN[0] + i*N;
  }

  /* posterior probability of b is from uniform distribution
   */
  bU  = (double**) Calloc(L, double*);
  bU[0] = (double*) Calloc(L*N, double);
  for(i=1; i<L; i++){
    bU[i] = bU[0] + i*N;
  }

  /* posterior probability of b is from normal distribution
   */
  bN  = (double**) Calloc(L, double*);
  bN[0] = (double*) Calloc(L*MM, double);
  for(i=1; i<L; i++){
    bN[i] = bN[0] + i*MM;
  }
  
  if(DEBUG){
    Rprintf("update transition matrix \n");
  }
  
  // update transition matrix
  // only trans_m_new[j][k], j!=k are estimated
  for(i=0; i<L-1; i++){
    dis[i] = pos[i+1] - pos[i];
  }
  
  if(estimate_trans_m){
    // Rprintf("estimate trans_m\n");
    for(j=0; j<N; j++){
      x = 0;
      for(i=0; i<S; i++){
        for(k=x; k<x+len[i]-1; k++){
          ps[k] = log(1-exp(-dis[k]/Ds[j]));
        }
        // note that ps[x + len[i] - 1] should never be used
        ps[x+len[i]-1] = 0.0;
        x = x + len[i];
      }

      for(k=0; k<N; k++){
        lvjk[k] = 0;
      }

      for(k=0; k<N; k++){
        if(k!=j){
          for(i=0; i<L-1; i++){
            if (dis[i] > *distThreshold) {
              vjk[i] = -1.0/0.0;
            }else {
              vjk[i] = f[i][j] + ps[i] + em[i+1][k] + b[i+1][k];
            }
          }

          x = 0;
          for(i=0; i<S; i++){
            for(z=x; z<x+len[i]-1; z++){
              vjk0[z-i] = vjk[z];
            }
            x = x + len[i];
          }

          logsumexp(vjk0, &LS, &lse);
          lvjk[k] = log(trans_m[j][k]) + lse;
          // Rprintf("k=%d,lse=%f\n",k,lse);
        }
      }

      for(k=0; k<N; k++){
        if(k<j){
          lvjk_j[k] = lvjk[k];
        }
        if(k>j){
          lvjk_j[k-1] = lvjk[k];
        }
      }

      logsumexp(lvjk_j, &NN, &lsum);
      
      // Rprintf("lsum=%f\n",lsum);
      // Rprintf("\n");

      if(lsum < log(EPSILON)){
        /* jth row of the transition matrix is close to 0, skip updating */
        for(k=0; k<N; k++){
          trans_m_new[j][k] = trans_m[j][k];
        }
      }else{
        for(k=0; k<N; k++){
          trans_m_new[j][k] = exp(lvjk[k] - lsum);
        }
      }
    }

    /* 
     * (1) calculate trans_m_new[j][j].
     * (2) update transition matrix so that the minimum 
     *     transition probability is min_tp. 
     */
    for(j=0; j<N; j++){
      transum0 = transum1 = 0.0;

      for(k=0; k<N; k++){
        if(k!=j){ 
          transum0 += trans_m_new[j][k];
          if(trans_m_new[j][k] < min_tp[k]){
            trans_m_new[j][k] = min_tp[k];
          }
          transum1 += trans_m_new[j][k];
        }
      }
      
      if(fabs(1 - transum0) > 1e-5){
        Rprintf("trans_m_new:\n");
        Rprint_me(trans_m_new,0,N,0,N);
        error("transition matrix row sum is not 1\n");
      }
      trans_m_new[j][j] = 0.0;
      
      for(k=0; k<N; k++){
        if(k!=j){ 
          trans_m_new[j][k] = trans_m_new[j][k]/transum1;
        }
      }
    }
    
    
    // update the parameter for initial probability
    x = 0;
    for(i=0; i<S; i++){
      for(j=0; j<N; j++){
        trans_begin_new[i][j] = pP[x][j];
      }
      x = x + len[i];
    }

  }

  if(DEBUG){
    Rprintf("update initial probability\n");
  }

  if(DEBUG){
    Rprintf("update parameters for LRR\n");
  }
 
  // update the parameters for mixture distribution of LRR
  for(i=0; i<L; i++){
    for(z=0; z<N; z++){
      if(eLRR[i][z] < EPSILON){
        rU[i][z] = rN[i][z] = 0.0;
      }else{
        if(is_snp_probe[i]){
          rU[i][z] = pP[i][z]*(pi_r[z]/R_m)/eLRR[i][z];
        }else{
          rU[i][z] = pP[i][z]*(pi_rcn[z]/R_m)/eLRR[i][z];
        }
        rN[i][z] = pP[i][z] - rU[i][z];
        if(rN[i][z] < 0){
          if(rN[i][z] < -EPSILON){ 
            Rprintf("rN[%d][%d] = %e\n", i, z, rN[i][z]);
            error("negative posterior probability of LRR\n"); 
          }
          rN[i][z] = 0.0;
        }
      }
    }
  }

  if(DEBUG){
    Rprintf("  estimate pi_r\n");
  }

  if(estimate_pi_r){
    // estimate pi_r, the proportion of uniform component
    for(z=0; z<N; z++){
      rUS[z] = 0.0;
      pPS[z] = 0.0;
      rUS_cn[z] = 0.0;
      pPS_cn[z] = 0.0;
      for(i=0; i<L; i++){
        if(is_snp_probe[i]){
          rUS[z] = rUS[z] + rU[i][z];
          pPS[z] = pPS[z] + pP[i][z];
        }else{
          rUS_cn[z] = rUS_cn[z] + rU[i][z];
          pPS_cn[z] = pPS_cn[z] + pP[i][z];
        }
      }
      if(pPS[z] < minStateP){
        pi_r_new[z] = pi_r[z];
      }else{
        pi_r_new[z] = rUS[z]/pPS[z];
        if(pi_r_new[z] > pi_r_upp){
          pi_r_new[z] = pi_r_upp;
        }
      }
      if(pPS_cn[z] < minStateP){
        pi_rcn_new[z] = pi_rcn[z];
      }else{
        pi_rcn_new[z] = rUS_cn[z]/pPS_cn[z];
        if(pi_rcn_new[z] > pi_r_upp){
          pi_rcn_new[z] = pi_r_upp;
        }
      }
    }
  }
  
  if(DEBUG){
    Rprintf("  estimate mu and sigma for LRR\n");
  }

  // estimate mu and sigma for LRR (r)
  zStart = 3;
  
  // combine state 1 and 2 here to estimate LRR
  sumz = 0.0;
  muz  = 0.0;
  sdz  = 0.0;
  sumz_cn = 0.0;
  muz_cn  = 0.0;
  sdz_cn  = 0.0;
  
  for(i=0; i<L; i++){
    rNcb[i] = rN[i][0] + rN[i][1];
    if(is_snp_probe[i]){
      sumz   += rNcb[i];
      muz    += rNcb[i]*LRR[i];
    }else{
      sumz_cn += rNcb[i];
      muz_cn  += rNcb[i]*LRR[i];
    }
  }
  
  if(sumz < minStateP){
    // warning("no update for state 1 and 2\n");
    mu_r_new[0] = mu_r_new[1] = mu_r[0];
    sd_r_new[0] = sd_r_new[1] = sd_r[0];
  }else{
    muz = muz/sumz;
    zRestrict(&muz, mu_r_lower[0], mu_r_upper[0]);
    for(i=0; i<L; i++){
      if(is_snp_probe[i]){
        sdz += rNcb[i]*(LRR[i]-muz)*(LRR[i]-muz);
      }
    }
    sdz = sqrt(sdz/sumz);
    zRestrict(&sdz, sd_r_lower[0], sd_r_upper[0]);
    mu_r_new[0] = mu_r_new[1] = muz;
    sd_r_new[0] = sd_r_new[1] = sdz;
  }

  if(sumz_cn < minStateP){
    // warning("no update for state 1 and 2\n");
    mu_rcn_new[0] = mu_rcn_new[1] = mu_rcn[0];
    sd_rcn_new[0] = sd_rcn_new[1] = sd_rcn[0];
  }else{
    muz_cn = muz_cn/sumz_cn;
    zRestrict(&muz_cn, mu_rcn_lower[0], mu_rcn_upper[0]);
    for(i=0; i<L; i++){
      if(!is_snp_probe[i]){
        sdz_cn += rNcb[i]*(LRR[i]-muz_cn)*(LRR[i]-muz_cn);
      }
    }
    sdz_cn = sqrt(sdz_cn/sumz_cn);
    zRestrict(&sdz_cn, sd_rcn_lower[0], sd_rcn_upper[0]);
    mu_rcn_new[0] = mu_rcn_new[1] = muz_cn;
    sd_rcn_new[0] = sd_rcn_new[1] = sdz_cn;
  }

  if(!CNA){ 
    zEnd = 6;
  }else{
    zEnd = 4;

    // combine state 5 and 6 to estimate LRR when CN=3
    sumz = 0.0;
    muz  = 0.0;
    sdz  = 0.0;
    sumz_cn = 0.0;
    muz_cn  = 0.0;
    sdz_cn  = 0.0;
    
    for(i=0; i<L; i++){
      rNcb[i] = rN[i][4] + rN[i][5];
      if(is_snp_probe[i]){
        sumz   += rNcb[i];
        muz    += rNcb[i]*LRR[i];
      }else{
        sumz_cn += rNcb[i];
        muz_cn  += rNcb[i]*LRR[i];
      }
    }

    if(sumz < minStateP){
      // warning("no update for state 5 and 6\n");
      mu_r_new[4] = mu_r_new[5] = mu_r[4];
      sd_r_new[4] = sd_r_new[5] = sd_r[4];
    }else{
      muz = muz/sumz;
      zRestrict(&muz, mu_r_lower[4], mu_r_upper[4]);
      for(i=0; i<L; i++){
        if(is_snp_probe[i]){
          sdz += rNcb[i]*(LRR[i]-muz)*(LRR[i]-muz);
        }
      }
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_r_lower[4], sd_r_upper[4]);
      mu_r_new[4] = mu_r_new[5] = muz;
      sd_r_new[4] = sd_r_new[5] = sdz;
    }
    
    if(sumz_cn < minStateP){
      // warning("no update for state 5 and 6\n");
      mu_rcn_new[4] = mu_rcn_new[5] = mu_rcn[4];
      sd_rcn_new[4] = sd_rcn_new[5] = sd_rcn[4];
    }else{
      muz_cn = muz_cn/sumz_cn;
      zRestrict(&muz_cn, mu_rcn_lower[4], mu_rcn_upper[4]);
      for(i=0; i<L; i++){
        if(!is_snp_probe[i]){
          sdz_cn += rNcb[i]*(LRR[i]-muz_cn)*(LRR[i]-muz_cn);
        }
      }
      sdz_cn = sqrt(sdz_cn/sumz_cn);
      zRestrict(&sdz_cn, sd_rcn_lower[4], sd_rcn_upper[4]);
      mu_rcn_new[4] = mu_rcn_new[5] = muz_cn;
      sd_rcn_new[4] = sd_rcn_new[5] = sdz_cn;
    }

    // combine state 7, 8, and 9 to estimate LRR when CN=4
    sumz = 0.0;
    muz  = 0.0;
    sdz  = 0.0;
    sumz_cn = 0.0;
    muz_cn  = 0.0;
    sdz_cn  = 0.0;

    for(i=0; i<L; i++){
      rNcb[i] = rN[i][6] + rN[i][7] + rN[i][8];
      if(is_snp_probe[i]){
        sumz   += rNcb[i];
        muz    += rNcb[i]*LRR[i];
      }else{
        sumz_cn += rNcb[i];
        muz_cn  += rNcb[i]*LRR[i];
      }
    }

    if(sumz < minStateP){
      // warning("no update for state 5 and 6\n");
      mu_r_new[6] = mu_r_new[7] = mu_r_new[8] = mu_r[6];
      sd_r_new[6] = sd_r_new[7] = sd_r_new[8] = sd_r[6];
    }else{
      muz = muz/sumz;
      zRestrict(&muz, mu_r_lower[6], mu_r_upper[6]);
      for(i=0; i<L; i++){
        if(is_snp_probe[i]){
          sdz += rNcb[i]*(LRR[i]-muz)*(LRR[i]-muz);
        }
      }
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_r_lower[6], sd_r_upper[6]);
      mu_r_new[6] = mu_r_new[7] = mu_r_new[8] = muz;
      sd_r_new[6] = sd_r_new[7] = sd_r_new[8] = sdz;
    }

    if(sumz_cn < minStateP){
      // warning("no update for state 5 and 6\n");
      mu_rcn_new[6] = mu_rcn_new[7] = mu_rcn_new[8] = mu_rcn[6];
      sd_rcn_new[6] = sd_rcn_new[7] = sd_rcn_new[8] = sd_rcn[6];
    }else{
      muz_cn = muz_cn/sumz_cn;
      zRestrict(&muz_cn, mu_rcn_lower[6], mu_rcn_upper[6]);
      for(i=0; i<L; i++){
        if(!is_snp_probe[i]){
          sdz_cn += rNcb[i]*(LRR[i]-muz_cn)*(LRR[i]-muz_cn);
        }
      }
      sdz_cn = sqrt(sdz_cn/sumz_cn);
      zRestrict(&sdz_cn, sd_rcn_lower[6], sd_rcn_upper[6]);
      mu_rcn_new[6] = mu_rcn_new[7] = mu_rcn_new[8] = muz_cn;
      sd_rcn_new[6] = sd_rcn_new[7] = sd_rcn_new[8] = sdz_cn;
    }

  }
  
  for(z=zStart-1; z<zEnd; z++){
    sumz = 0.0;
    muz  = 0.0;
    sdz  = 0.0;
    sumz_cn = 0.0;
    muz_cn  = 0.0;
    sdz_cn  = 0.0;

    for(i=0; i<L; i++){
      if(is_snp_probe[i]){
        sumz += rN[i][z];
        muz  += rN[i][z]*LRR[i];
      }else{
        sumz_cn += rN[i][z];
        muz_cn  += rN[i][z]*LRR[i];
      }
    }

    if(sumz < minStateP){
      // warning("no update for state %d\n",z+1);
      mu_r_new[z] = mu_r[z];
      sd_r_new[z] = sd_r[z];
    }else{
      muz = muz/sumz;
      zRestrict(&muz, mu_r_lower[z], mu_r_upper[z]);
      for(i=0; i<L; i++){
        if(is_snp_probe[i]){
          sdz += rN[i][z]*(LRR[i]-muz)*(LRR[i]-muz);
        }
      }
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_r_lower[z], sd_r_upper[z]);
      sd_r_new[z] = sdz;
      mu_r_new[z] = muz;
    }

    if(sumz_cn < minStateP){
      // warning("no update for state %d\n",z+1);
      mu_rcn_new[z] = mu_rcn[z];
      sd_rcn_new[z] = sd_rcn[z];
    }else{
      muz_cn = muz_cn/sumz_cn;
      zRestrict(&muz_cn, mu_rcn_lower[z], mu_rcn_upper[z]);
      for(i=0; i<L; i++){
        if(!is_snp_probe[i]){
          sdz_cn += rN[i][z]*(LRR[i]-muz_cn)*(LRR[i]-muz_cn);
        }
      }
      sdz_cn = sqrt(sdz_cn/sumz_cn);
      zRestrict(&sdz_cn, sd_rcn_lower[z], sd_rcn_upper[z]);
      sd_rcn_new[z] = sdz_cn;
      mu_rcn_new[z] = muz_cn;
    }
  }
      
  if(DEBUG){
      Rprintf("update parameters for BAF\n");
  }
  
  // update the parameters for mixture distribuiton of BAF
  // where 0<BAF<1

  flag = 0;
  for(i=0; i<L; i++){
    if(BAF[i]>EPSILON && BAF[i]<1.0-EPSILON && is_snp_probe[i]){
      Omiga[flag] = i;
      BAF01[flag]= BAF[i];
      flag++;
    }
  }

  /* bU[i][z] = postP(probe i from uniform component of state z) 
   * where i is in Omiga
   */

  for(i=0; i<flag; i++){
    Oi = Omiga[i];
    for(z=0; z<N; z++){
      if(eBAF[Oi][z] < EPSILON){
        bU[Oi][z] = 0.0;
      }else{
        bU[Oi][z] = pi_b[z]*pP[Oi][z]/eBAF[Oi][z];
      }
    }
  }

  if(DEBUG){
    Rprintf("  estimate pi_b\n");
  }
  
  if(estimate_pi_b){
    // estimate pi for BAF (b)
    for(z=0; z<N; z++){
      bUS[z] = 0.0;
      pPS[z] = 0.0;
      
      for(i=0; i<L; i++){
        if(is_snp_probe[i]){
          bUS[z] += bU[i][z];
          pPS[z] += pP[i][z];
        }
      }
    
      if(pPS[z] < minStateP){
        pi_b_new[z] = pi_b[z];
      }else{
        pi_b_new[z] = bUS[z]/pPS[z];
        if(z!=2 && pi_b_new[z] > pi_b_upp){
          pi_b_new[z] = pi_b_upp;
        }
      }
    }
  }
  
  if(DEBUG){
    Rprintf("  estimate mu and sigma for BAF \n");
  }

  // estimate mu and sigma for BAF (b)
  /*
   * note bN is a matrix of size L*MM, where MM=36 (CNA) or 21 (CNV), 
   * one state corresponds to a group of consecutive columns in bN
   */
  k = 0;

  for(i=0; i<flag; i++){
    Oi  = Omiga[i];
    baf = BAF[Oi];
    pbf = pBs[Oi];
    nGtp = normalGtp[Oi];
    
    for(z=0; z<N; z++){
      weights(ws, pbf, z, CNA, contam, nGtp, *geno_error);
      
      for(hh=0; hh<H[z]; hh++){
        k  = zIdxStart[z] + hh;
        if(eBAF[Oi][z] < EPSILON){
          bN[i][k] = 0.0;
        }else{
          tmp = ws[hh]*dnorm(baf,mu_b[z][hh],sd_b[z][hh],0);
          bN[i][k] = (1.0 - pi_b[z])*tmp*pP[Oi][z]/eBAF[Oi][z];
        }
      }
    }

  }

/**
*  the columns of matrix bN
*
*  HMM state  Genotype state
*  1          1 -> AA      2 -> AB        3 -> BB
*  2          4 -> AA      5 -> (AA,AB)   6 -> (BB,AB)   7 -> BB
*  3          8 -> NULL
*  4          9 -> A      10 -> (A,AB)   11 -> (B,AB)   12 -> B
*  5          13 -> AAA   14 -> AAB      15 -> ABB      16 -> BBB
*  6/CNV      17 -> AAAA  18 -> AAAB     19 -> AABB     20 -> ABBB 21 -> BBBB
*  6/CNA      22 -> AAA   23 -> (AAA,AB) 24 -> (BBB,AB) 25 -> BBB
*  7          26 -> AAAA  27 -> AABB     28 -> BBBB
*  8          29 -> AAAA  30 -> (AAAA,AB)31 -> (BBBB,AB)32 -> BBBB
*  9          33 -> AAAA  34 -> AAAB     35 -> ABBB     36 -> BBBB
**/

  if(DEBUG){
    Rprintf("    BAF combine state12 \n");
  }
  
  /**
   * combine state 1 2 for BAF parameter estimation
   */
  zStart = 3;

  /* AB */
  sumz = 0.0;
  muz  = 0.0;
  sdz  = 0.0;
  k = 1;
  for(j=0; j<flag; j++){
    sumz += bN[j][k];
    muz  += bN[j][k]*BAF01[j];
  }
  if(sumz < minStateP){
    mu_b_new[0][1] = mu_b[0][1];
    sd_b_new[0][1] = sd_b[0][1];
  }else{
    muz = muz/sumz;
    zRestrict(&muz, mu_b_lower[0][1], mu_b_upper[0][1]);
    for(j=0; j<flag; j++){
      sdz += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
    }
    sdz = sqrt(sdz/sumz);
    zRestrict(&sdz, sd_b_lower[0][1], sd_b_upper[0][1]);
    mu_b_new[0][1] = muz;
    sd_b_new[0][1] = sdz;
  }

  /* AA */
  sdz  = 0.0;
  sumz = 0.0;
  k1 = 0;
  k2 = 3;
  for(j=0; j<flag; j++){
    bNcb[j] = bN[j][k1] + bN[j][k2];
    sdz  += bNcb[j]*BAF01[j]*BAF01[j];
    sumz += bNcb[j];
  }

  mu_b_new[0][0] = mu_b_new[1][0] = 0.0;
  if(sumz < minStateP){
    sd_b_new[0][0] = sd_b_new[1][0] = sd_b[0][0];
  }else{
    sdz = sqrt(sdz/sumz);
    zRestrict(&sdz, sd_b_lower[0][0], sd_b_upper[0][0]);
    sd_b_new[0][0] = sd_b_new[1][0] = sdz;            
  }

  /* BB */
  sdz  = 0.0;
  sumz = 0.0;
  k1 = 2;
  k2 = 6;
  for(j=0; j<flag; j++){
    bNcb[j] = bN[j][k1] + bN[j][k2];
    sdz    += bNcb[j]*(BAF01[j]-1.0)*(BAF01[j]-1.0);
    sumz   += bNcb[j];
  }

  mu_b_new[0][2] = mu_b_new[1][3] = 1.0;
  if(sumz < minStateP){
    sd_b_new[0][2] = sd_b_new[1][3] = sd_b[0][2];
  }else{
    sdz = sqrt(sdz/sumz);
    zRestrict(&sdz, sd_b_lower[0][2], sd_b_upper[0][2]);
    sd_b_new[0][2] = sd_b_new[1][3] = sdz;
  }

  if(DEBUG){
    Rprintf("    BAF combine states for CNA \n");
  }
  
  /**
   * combine some states for BAF parameter estimation
   */
  
  if (!CNA) {
    zEnd = 6;
  }else {
    zEnd = 4;
    
    /* AAA state 5 and 6 */
    sdz  = 0.0;
    sumz = 0.0;
    k1 = 12;
    k2 = 21;
    for(j=0; j<flag; j++){
      bNcb[j] = bN[j][k1] + bN[j][k2];
      sdz    += bNcb[j]*BAF01[j]*BAF01[j];
      sumz   += bNcb[j];
    }

    mu_b_new[4][0] = mu_b_new[5][0] = 0.0;
    if(sumz < minStateP){         
      sd_b_new[4][0] = sd_b_new[5][0] = sd_b[4][0];
    }else{
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_b_lower[4][0], sd_b_upper[4][0]);
      sd_b_new[4][0] = sd_b_new[5][0] = sdz;
    }

    /* BBB state 5 and 6 */
    sdz  = 0.0;
    sumz = 0.0;
    k1 = 15;
    k2 = 24;
    for(j=0; j<flag; j++){
      bNcb[j] = bN[j][k1] + bN[j][k2];
      sdz    += bNcb[j]*(BAF01[j] - 1.0)*(BAF01[j] - 1.0);
      sumz   += bNcb[j];
    }

    mu_b_new[4][3] = mu_b_new[5][3] = 1.0;
    if(sumz < minStateP){         
      sd_b_new[4][3] = sd_b_new[5][3] = sd_b[4][3];
    }else{
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_b_lower[4][3], sd_b_upper[4][3]);
      sd_b_new[4][3] = sd_b_new[5][3] = sdz;
    }

    /* AAB, ABB, AAAB, ABBB */
    k1s[0]=13, k1s[1]=14, k1s[2]=33, k1s[3]=34;
    rid[0]=4,  rid[1]=4,  rid[2]=8,  rid[3]=8;
    cid[0]=1,  cid[1]=2,  cid[2]=1,  cid[3]=2; 
    
    for(i=0; i<4; i++){
      sumz = 0.0;
      muz  = 0.0;
      sdz  = 0.0;

      k1 = k1s[i];
      ri = rid[i];
      ci = cid[i];
      
      for(j=0; j<flag; j++){
        bNcb[j] = bN[j][k1];
        sumz += bNcb[j];
        muz  += bNcb[j]*BAF01[j];
      }
      
      if(sumz < minStateP){
        mu_b_new[ri][ci] = mu_b[ri][ci];
        sd_b_new[ri][ci] = sd_b[ri][ci];
      }else{
        muz = muz/sumz;
        zRestrict(&muz, mu_b_lower[ri][ci], mu_b_upper[ri][ci]);
        for(j=0; j<flag; j++){
          sdz += bNcb[j]*(BAF01[j]-muz)*(BAF01[j]-muz);
        }
        sdz = sqrt(sdz/sumz);
        zRestrict(&sdz, sd_b_lower[ri][ci], sd_b_upper[ri][ci]);
        mu_b_new[ri][ci] = muz;
        sd_b_new[ri][ci] = sdz;
      }
    }
    
    /* check the symetric of AAB vs. ABB, and AAAB vs. ABBB*/
    for(ri=4; ri<=8; ri+=4){
      d0 = mu_b_new[ri][1];
      d1 = 1.0 - mu_b_new[ri][2];
      if(d1 - d0 > max_diff){
        dajust = 0.5*(d1 - d0 - max_diff);
        mu_b_new[ri][1] = mu_b_new[ri][1] + dajust;
        mu_b_new[ri][2] = mu_b_new[ri][2] + dajust;            
      }else if(d0 - d1 > max_diff){
        dajust = 0.5*(d0 - d1 - max_diff);
        mu_b_new[ri][1] = mu_b_new[ri][1] - dajust;
        mu_b_new[ri][2] = mu_b_new[ri][2] - dajust;            
      }
    }
    
    /* AABB */
    sumz = 0.0;
    muz  = 0.0;
    sdz  = 0.0;
    k    = 26;
    for(j=0; j<flag; j++){
      sumz += bN[j][k];
      muz  += bN[j][k]*BAF01[j];
    }
    if(sumz < minStateP){
      sd_b_new[6][1] = sd_b[6][1];
      mu_b_new[6][1] = mu_b[6][1];
    }else{
      muz = muz/sumz;
      zRestrict(&muz, mu_b_lower[6][1], mu_b_upper[6][1]);
      for(j=0; j<flag; j++){
          sdz +=  bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
      }
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_b_lower[6][1], sd_b_upper[6][1]);
      sd_b_new[6][1] = sdz;
      mu_b_new[6][1] = muz;
    }

    /* AAAA */
    sdz  = 0.0;
    sumz = 0.0;
    k1 = 25;
    k2 = 28;
    k3 = 32;
    for(j=0; j<flag; j++){
      bNcb[j] = bN[j][k1] + bN[j][k2] + bN[j][k3];
      sdz    += bNcb[j]*BAF01[j]*BAF01[j];
      sumz   += sumz + bNcb[j];
    }

    mu_b_new[6][0] = mu_b_new[7][0] = mu_b_new[8][0] = 0.0;
    if(sumz < minStateP){         
      sd_b_new[6][0] = sd_b_new[7][0] = sd_b_new[8][0] = sd_b[6][0];
    }else{
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_b_lower[6][0], sd_b_upper[6][0]);
      sd_b_new[6][0] = sd_b_new[7][0] = sd_b_new[8][0] = sdz;
    }

    /* BBBB */
    sdz  = 0.0;
    sumz = 0.0;
    k1 = 27;
    k2 = 31;
    k3 = 35;
    for(j=0; j<flag; j++){
      bNcb[j] = bN[j][k1] + bN[j][k2] + bN[j][k3];
      sdz    += bNcb[j]*(BAF01[j] -1.0)*(BAF01[j] - 1.0);
      sumz   += bNcb[j];
    }

    mu_b_new[6][2] = mu_b_new[7][3] = mu_b_new[8][3] = 1.0;
    if(sumz < minStateP){         
      sd_b_new[6][2] = sd_b_new[7][3] = sd_b_new[8][3] = sd_b[6][2];
    }else{
      sdz = sqrt(sdz/sumz);
      zRestrict(&sdz, sd_b_lower[6][2], sd_b_upper[6][2]);
      sd_b_new[6][2] = sd_b_new[7][3] = sd_b_new[8][3] = sdz;
    }
  }
  
  if(CNA && contam){
    if(DEBUG){
      Rprintf("    If CNA and there is contamination");
      Rprintf(" handle the mixture components for state 2, 4, 6, and 8\n");
    }

    // for state 2, 4, 6, 8, handle the mixture components
    
    for(z=1; z<8; z+=2){
      // hh = 1, ie., the mixture of (AA, AB), (A, AB), (AAA, AB), or (AAAA, AB)
      hh = 1;
      k  = zIdxStart[z] + hh;

      //sumz is summation across those probes with 0<BAF<1
      sumz = 0.0; 
      for(j=0; j<flag; j++){
        sumz += bN[j][k];
      }
      
      if(sumz < minStateP){
        mu_b_new[z][hh] = mu_b[z][hh];
        sd_b_new[z][hh] = sd_b[z][hh];                	
      }else{
        M1 = M2 = 0.0;
        for(j=0; j<flag; j++){
          M1 += bN[j][k]*BAF01[j];
          M2 += bN[j][k]*BAF01[j]*BAF01[j];
        } 
        M1 = M1/sumz;
        M2 = M2/sumz;
        
        // sum0 is the sum of the probability of being state z
        // and normalGtp is AB and BAF = 0
        sum0 = 0.0;
        for (i=0; i<L; i++) {
          if (is_snp_probe[i]) {
            if (BAF[i] < EPSILON) {
              if(normalGtp[i] == 1){
                sum0 += pP[i][z];
              }else if (normalGtp[i] == -1){
                weights(ws, pbf, z, CNA, contam, nGtp, *geno_error);
                pnz0  = pnorm(0.0, mu_b[z][1], sd_b[z][1], 1, 0);
                sum0 += pP[i][z]*pnz0*ws[1]/(pnz0*ws[1] + 0.5*ws[0]);
              }
            }
          }
        }
        p = sumz/(sumz + sum0);  T = 0.0;  left = 1;
        if(p < 0.9999){
          tnorm_mle(T, M1, M2, p, &muz, &sdz, left, roundoff, maxIt);
        }else{
          muz = M1;
          sdz = sqrt(M2 - M1*M1);
        }
        
        /*
        if(DEBUG && z==5){
          Rprintf("z=%d, hh=%d, p=%e\n", z, hh, p);
          Rprintf("  sumz=%e, sum0=%e, sum1=%e\n", sumz, sum0, sum1);
          Rprintf("  M1=%e, M2=%e, muz=%e, sdz=%e\n", M1, M2, muz, sdz);
        }
        */
        
        // if muz need to be restricted, sdz need to be recalculated 
        if(muz < mu_b_lower[z][hh] || muz > mu_b_upper[z][hh]){
          zRestrict(&muz, mu_b_lower[z][hh], mu_b_upper[z][hh]);
          // update sdz given the new muz
          if(p < 0.9999){
            sdz  = 0.0;
            sumz = 0.0;
            for(j=0; j<flag; j++){
              if(BAF01[j] > muz){
                sumz += bN[j][k];
                sdz  += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
              }
            }
          }else{
            sdz = 0.0;
            for(j=0; j<flag; j++){
              sdz  += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
            }
          }
          sdz = sqrt(sdz/sumz);
          zRestrict(&sdz,sd_b_lower[z][hh],sd_b_upper[z][hh]);         
        }
        /*
        if (DEBUG && z==5) {
          Rprintf("  After restriction: muz=%e, sdz=%e\n", muz, sdz);
        }
        */
        
        mu_b_new[z][hh] = muz; 
        sd_b_new[z][hh] = sdz;
      }
      
      // hh = 2, ie., the mixture of (BB, AB), (B, AB), (BBB, AB), or (BBBB, AB)
      hh = 2;
      k  = zIdxStart[z] + hh;
      
      //sumz is summation across those probes with 0<BAF<1
      sumz = 0.0; 
      for(j=0; j<flag; j++){
        sumz += bN[j][k];
      }
      
      if(sumz < minStateP){
        mu_b_new[z][hh] = mu_b[z][hh];
        sd_b_new[z][hh] = sd_b[z][hh];                	
      }else{
        M1 = M2 = 0.0;
        for(j=0; j<flag; j++){
          M1 += bN[j][k]*BAF01[j];
          M2 += bN[j][k]*BAF01[j]*BAF01[j];
        } 
        M1 = M1/sumz;
        M2 = M2/sumz;
        
        // sum1 is the sum of the probability of being state z
        // and normalGtp is AB and BAF = 1
        sum1 = 0.0;
        for (i=0; i<L; i++) {
          if (is_snp_probe[i]) {
            if (BAF[i] > 1.0 - EPSILON) {
              if(normalGtp[i] == 1){
                sum1 += pP[i][z];
              }else if (normalGtp[i] == -1){
                weights(ws, pbf, z, CNA, contam, nGtp, *geno_error);
                pnz1  = pnorm(1.0, mu_b[z][2], sd_b[z][2], 0, 0);
                sum1 += pP[i][z]*pnz1*ws[2]/(pnz1*ws[2] + 0.5*ws[3]);
              }
            }
          }
        }
        
        p = sumz/(sumz + sum1);  T = 1.0;  left = 0;
        if(p < 0.9999){
          tnorm_mle(T, M1, M2, p, &muz, &sdz, left, roundoff, maxIt);
        }else{
          muz = M1;
          sdz = sqrt(M2 - M1*M1);
        }
        
        /*
         Rprintf("z=%d, hh=%d, p=%e\n", z, hh, p);
         Rprintf("  sumz=%e, sum0=%e, sum1=%e\n", sumz, sum0, sum1);
         Rprintf("  M1=%e, M2=%e, muz=%e, sdz=%e\n", M1, M2, muz, sdz);
         */
        
        // if muz need to be restricted, 
        // sdz need to be recalculated 
        if(muz < mu_b_lower[z][hh] || muz > mu_b_upper[z][hh]){
          zRestrict(&muz, mu_b_lower[z][hh], mu_b_upper[z][hh]);
          // update sdz given the new muz
          if(p < 0.9999){
            sdz  = 0.0;
            sumz = 0.0;
            for(j=0; j<flag; j++){
              if(BAF01[j] < muz){
                sumz += bN[j][k];
                sdz  += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
              }
            }
          }else{
            sdz = 0.0;
            for(j=0; j<flag; j++){
              sdz  += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
            }
          }
          sdz = sqrt(sdz/sumz);
          zRestrict(&sdz,sd_b_lower[z][hh],sd_b_upper[z][hh]);         
        }
        
        mu_b_new[z][hh] = muz; 
        sd_b_new[z][hh] = sdz;
        /*
        Rprintf("  After restriction: muz=%e, sdz=%e\n", muz, sdz);        
        */
      }
      
      d0 = mu_b_new[z][1];
      d1 = 1.0 - mu_b_new[z][2];
      if(d1 - d0 > max_diff){
        dajust = 0.5*(d1 - d0 - max_diff);
        mu_b_new[z][1] = mu_b_new[z][1] + dajust;
        mu_b_new[z][2] = mu_b_new[z][2] + dajust;            
      }else if(d0 - d1 > max_diff){
        dajust = 0.5*(d0 - d1 - max_diff);
        mu_b_new[z][1] = mu_b_new[z][1] - dajust;
        mu_b_new[z][2] = mu_b_new[z][2] - dajust;            
      }
    }
  }

  if(DEBUG){
    Rprintf("    BAF other states \n");
  }

  // k is the column index,  
  // if zStart > 1, we need to skip some columns
  k = -1;
  if(zStart > 1){
    for(z=0; z<zStart-1; z++){
      k += H[z];
    }
  }
  
  for(z=zStart-1; z<zEnd; z++){
    for(hh=0; hh<H[z]; hh++){
      k++;
      // state 2/4, mixture of (A/AA, AB) and (B/BB, AB)
      if((z==1 || z==3) && contam && (hh==1 || hh==2)){
        continue;
      }
      sumz = 0.0;
      for(j=0; j<flag; j++){
        sumz += bN[j][k];
      }
      
      if(sumz < minStateP){
        mu_b_new[z][hh] = mu_b[z][hh];
        sd_b_new[z][hh] = sd_b[z][hh];
      }else{
        if(z==2){
          muz = 0.5;
        }else{
          if(hh==0){
            muz = 0.0;
          }else if(hh==H[z]-1){
            muz = 1.0;
          }else{
            muz = 0.0;
            for(j=0; j<flag; j++){
              muz += bN[j][k]*BAF01[j];
            }
            muz = muz/sumz;
            zRestrict(&muz,mu_b_lower[z][hh],mu_b_upper[z][hh]);
          }
        }
        sdz = 0.0;
        for(j=0; j<flag; j++){
          sdz += bN[j][k]*(BAF01[j]-muz)*(BAF01[j]-muz);
        }
        sdz = sqrt(sdz/sumz);
        zRestrict(&sdz,sd_b_lower[z][hh],sd_b_upper[z][hh]);
        mu_b_new[z][hh] = muz;
        sd_b_new[z][hh] = sdz;
      }           
    }
  }
  
  if(!CNA){
    /* check the symetric of AAB vs. ABB*/
    ri = 4; //state 5
    d0 = mu_b_new[ri][1];
    d1 = 1.0 - mu_b_new[ri][2];
    if(d1 - d0 > max_diff){
      dajust = 0.5*(d1 - d0 - max_diff);
      mu_b_new[ri][1] = mu_b_new[ri][1] + dajust;
      mu_b_new[ri][2] = mu_b_new[ri][2] + dajust;            
    }else if(d0 - d1 > max_diff){
      dajust = 0.5*(d0 - d1 - max_diff);
      mu_b_new[ri][1] = mu_b_new[ri][1] - dajust;
      mu_b_new[ri][2] = mu_b_new[ri][2] - dajust;            
    }

    /* check the symetric of AAAB vs. ABBB*/
    ri = 5; //state 6
    d0 = mu_b_new[ri][1];
    d1 = 1.0 - mu_b_new[ri][3];
    if(d1 - d0 > max_diff){
      dajust = 0.5*(d1 - d0 - max_diff);
      mu_b_new[ri][1] = mu_b_new[ri][1] + dajust;
      mu_b_new[ri][3] = mu_b_new[ri][3] + dajust;            
    }else if(d0 - d1 > max_diff){
      dajust = 0.5*(d0 - d1 - max_diff);
      mu_b_new[ri][1] = mu_b_new[ri][1] - dajust;
      mu_b_new[ri][3] = mu_b_new[ri][3] - dajust;            
    }
  }
  
  if(CNA && k != 11){
    error("CNA=%d, k=%d, MM=%d\n", CNA, k, MM);
  }

  if(!CNA && k != 20){
    error("CNA=%d, k=%d, MM=%d\n", CNA, k, MM);
  }

  // Rprintf("Free memory\n");

  Free(dis);
  Free(vjk);
  Free(vjk0);
  Free(lvjk);
  Free(lvjk_j);
  Free(ps);

  Free(rUS);
  Free(pPS);
  Free(rUS_cn);
  Free(pPS_cn);
  Free(rNcb);

  Free(bUS);
  Free(Omiga);
  Free(BAF01);
  Free(bNcb);

  Free(rU[0]);
  Free(rU);
  Free(rN[0]);
  Free(rN);
  Free(bU[0]);
  Free(bU);
  Free(bN[0]);
  Free(bN);
}
