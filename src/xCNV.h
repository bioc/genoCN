/**********************************************************************
 *
 * xCNV.h
 *
 * copyright (c) 2008, Wei Sun, UNC-CH
 *
 * last modified Aug 29, 2008
 * first written Jun 02, 2008
 *
 * C functions for the xCNV package
 *
 * Hidden Markov Model
 *
 **********************************************************************/

int weights(double* ws, double pbf, int z, int CNA, int contamination, 
            int normalGtp, double geno_error);

double emiss4(double b, double EPSILON, double* ws, double* mu_b, 
              double* sd_b, double pi_b, int contam);

double emissK(double b, double EPSILON, double* ws, double* mu_b, 
              double* sd_b, double pi_b, int Hz);

void emiss(double* LRR, double* BAF, double* pBs, int* normalGtp,
           double* geno_error, double* pi_r, double* mu_r, double *sd_r, 
           double *pi_rcn, double *mu_rcn, double *sd_rcn, double *R_mR, 
           double *pi_b, double *mu_bR, double *sd_bR, int *dims, 
           int* is_snp_probe, int* CNAR, int* contamR, double *eLRRR, 
           double *eBAFR, double *emR);


void transition_c(double** trans_m, double di, double* Ds, int N,
                  double** trans1, double* transB, double distThreshold);

void viterbi(double* pos, double* emR, int * em_inf, 
             double* trans_mR, double* trans_begin, 
             double* Ds, double* distThreshold, double* transB, 
             int* dims, int* path, double *logP);

void backward(double* pos, double* emR, int* em_inf, int* Rem_noinf,
              double* trans_mR, double* Ds, double* bR, int* dims, 
              double* distThreshold, double* transB);

void forward(double* pos, double* emR, int* em_inf, int* Rem_noinf,
             double* trans_mR, double* trans_begin, double* Ds,
             double* fR, int* dims, double* distThreshold, 
             double* transB);

void postP(double* fR, double* bR, int* f_inf, int* b_inf, int* no_inf,
           int* dims, double* pPR);

void baum_welch(double* pos, double* LRR, double* BAF,
                double* pBs, int* normalGtp, double* geno_error, 
                int* is_snp_probe, int* RCNA, double* Ds, 
                double* RR_m, double* pi_r, double* pi_rcn, 
                double* mu_r, double* mu_r_upper, double* mu_r_lower,
                double* sd_r, double* sd_r_upper, double* sd_r_lower,
                double* mu_rcn, double* mu_rcn_upper, double* mu_rcn_lower,
                double* sd_rcn, double* sd_rcn_upper, double* sd_rcn_lower,
                double* pi_b, double* trans_mR, double* min_tpR, 
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
                double* trans_m_newR, double* trans_begin_newR);
