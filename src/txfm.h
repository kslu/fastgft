#ifndef TXFM_H_
#define TXFM_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define CONFIG_DEBUG 0
#define CONFIG_USE_LIFTING_FOR_GIVENS 1
#define MAX_NUM_DATA 2000
#define BATCH_SIZE 1000
#define MAX_N 100

#define SQRT2 1.414213564
#define INVSQRT2 0.7071067812

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n);
#endif

void gft_star10_mat(const double *input, double *output);
void gft_star10_btf(const double *input, double *output);
void gft_star100_mat(const double *input, double *output);
void gft_star100_btf(const double *input, double *output);

void gft_cycle12_mat(const double *input, double *output);
void gft_cycle12_btf(const double *input, double *output);
void gft_cycle80_mat(const double *input, double *output);
void gft_cycle80_btf(const double *input, double *output);

void gft_bd4x4_mat(const double *input, double *output);
void gft_bd4x4_btf(const double *input, double *output);
void gft_bd8x8_mat(const double *input, double *output);
void gft_bd8x8_btf(const double *input, double *output);
void gft_bd8x8_tj(const double *input, double *output, int n_givens);
void gft_bd8x8_ptj(const double *input, double *output, int n_layers);
void gft_bd8x8_btf_tj(const double *input, double *output, int n_givens);
void gft_bd8x8_btf_ptj(const double *input, double *output, int n_layers);

void gft_dct4x4_mat(const double *input, double *output);
void gft_dct4x4_btf(const double *input, double *output);
void gft_dct4x4_sep(const double *input, double *output);
void gft_dct8x8_mat(const double *input, double *output);
void gft_dct8x8_btf(const double *input, double *output);
void gft_dct8x8_sep(const double *input, double *output);

void gft_skeleton15_mat(const double *input, double *output);
void gft_skeleton15_btf(const double *input, double *output);
void gft_skeleton25_mat(const double *input, double *output);
void gft_skeleton25_btf(const double *input, double *output);

void gft_z4x4_mat(const double *input, double *output);
void gft_z4x4_btf(const double *input, double *output);
void gft_z8x8_mat(const double *input, double *output);
void gft_z8x8_btf(const double *input, double *output);
void gft_z8x8_tj(const double *input, double *output, int n_givens);
void gft_z8x8_ptj(const double *input, double *output, int n_layers);
void gft_z8x8_btf_tj(const double *input, double *output, int n_givens);
void gft_z8x8_btf_ptj(const double *input, double *output, int n_layers);

#endif
