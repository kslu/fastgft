#ifndef TXFM_H_
#define TXFM_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define CONFIG_DEBUG 0
#define MAX_NUM_DATA 20000

#define SQRT2 1.414213564
#define INVSQRT2 0.7071067812

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n);
#endif
void mat_times_vec(const double *input, double *output, const double *gftmtx,
                   int n);

void gft_star10_mat(const double *input, double *output);
void gft_star10_btf(const double *input, double *output);
void gft_bd4x4_mat(const double *input, double *output);
void gft_bd4x4_btf(const double *input, double *output);
void gft_dct4x4_mat(const double *input, double *output);
void gft_dct4x4_btf(const double *input, double *output);
void gft_skeleton15_mat(const double *input, double *output);
void gft_skeleton15_btf(const double *input, double *output);

#endif
