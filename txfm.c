#include "txfm.h"
#include "gftmatrices.h"

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n) {
  fprintf(stderr, "(");
  for (int i = 0; i < n - 1; i++)
    fprintf(stderr, "%.4lf,", input[i]);
  fprintf(stderr, "%.4lf)\n[", input[n - 1]);
  for (int i = 0; i < n - 1; i++)
    fprintf(stderr, "%.4lf,", output[i]);
  fprintf(stderr, "%.4lf]\n", output[n - 1]);
}
#endif

void mat_times_vec(const double *input, double *output, const double *gftmtx,
                   int n) {
  for (int i = 0; i < n; i++) {
    output[i] = 0;
    for (int j = 0; j < n; j++)
      output[i] += input[j] * gftmtx[i * n + j];
  }
}

void gft_star10_mat(const double *input, double *output) {
  mat_times_vec(input, output, star10, 10);
#if CONFIG_DEBUG
  show_coefficients(input, output, 10);
#endif
}

void gft_star10_btf(const double *input, double *output) {
  double temp[10] = {0.0};
  double temp1[6] = {0.0};
  double temp2[4] = {0.0};
  double temp_out[10] = {0.0};

  int output_idx[10] = {0, 1, 9, 2, 3, 4, 5, 6, 7, 8};

  // stage 1
  temp2[0] = 2 * SQRT2 * input[0];
  temp2[1] = 2 * SQRT2 * input[1];
  temp[2] = input[2] + input[9];
  temp[9] = input[2] - input[9];
  temp[3] = input[3] + input[8];
  temp[8] = input[3] - input[8];
  temp[4] = input[4] + input[7];
  temp[7] = input[4] - input[7];
  temp[5] = input[5] + input[6];
  temp[6] = input[5] - input[6];

  // stage 2
  temp1[2] = temp[2] + temp[5];
  temp1[5] = temp[2] - temp[5];
  temp1[3] = temp[3] + temp[4];
  temp1[4] = temp[3] - temp[4];
  temp_out[6] = INVSQRT2 * temp[6];
  temp_out[7] = INVSQRT2 * temp[7];
  temp_out[8] = INVSQRT2 * temp[8];
  temp_out[9] = INVSQRT2 * temp[9];

  // stage 3
  temp2[2] = temp1[2] + temp1[3];
  temp2[3] = temp1[2] - temp1[3];
  temp_out[4] = temp1[4] / 2;
  temp_out[5] = temp1[5] / 2;

  // stage 4
  mat_times_vec(temp2, temp_out, star10_ppp, 3);
  temp_out[3] = temp2[3] * INVSQRT2 / 2;

  // output
  for (int i = 0; i < 10; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 10);
#endif
}

void gft_bd4x4_mat(const double *input, double *output) {
  mat_times_vec(input, output, bd4x4, 16);
#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void gft_bd4x4_btf(const double *input, double *output) {
  double temp[16] = {0.0};
  double temp1[16] = {0.0};
  double temp_out[16] = {0.0};

  int stage3_idx[16] = {0, 1, 2, 3, 5, 6, 7, 10, 11, 15, 4, 8, 9, 12, 13, 14};
  int output_idx[16] = {0, 3, 6, 8, 12, 15, 2, 7, 10, 14, 1, 4, 9, 13, 5, 11};

  // stage 1
  temp[0] = SQRT2 * input[0];
  temp[5] = SQRT2 * input[5];
  temp[10] = SQRT2 * input[10];
  temp[15] = SQRT2 * input[15];
  temp[1] = input[1] + input[4];
  temp[4] = input[1] - input[4];
  temp[2] = input[2] + input[8];
  temp[8] = input[2] - input[8];
  temp[3] = input[3] + input[12];
  temp[12] = input[3] - input[12];
  temp[6] = input[6] + input[9];
  temp[9] = input[6] - input[9];
  temp[7] = input[7] + input[13];
  temp[13] = input[7] - input[13];
  temp[11] = input[11] + input[14];
  temp[14] = input[11] - input[14];

  // stage 2
  temp1[3] = SQRT2 * temp[3];
  temp1[6] = SQRT2 * temp[6];
  temp1[9] = SQRT2 * temp[9];
  temp1[12] = SQRT2 * temp[12];
  temp1[0] = temp[0] + temp[15];
  temp1[15] = temp[0] - temp[15];
  temp1[1] = temp[1] + temp[11];
  temp1[11] = temp[1] - temp[11];
  temp1[2] = temp[2] + temp[7];
  temp1[7] = temp[2] - temp[7];
  temp1[5] = temp[5] + temp[10];
  temp1[10] = temp[5] - temp[10];
  temp1[4] = temp[4] + temp[14];
  temp1[14] = temp[4] - temp[14];
  temp1[8] = temp[8] + temp[13];
  temp1[13] = temp[8] - temp[13];

  // permutation
  for (int i = 0; i < 16; i++)
    temp[i] = temp1[stage3_idx[i]];

  // stage 3
  mat_times_vec(temp, temp_out, bd4x4_pp, 6);
  mat_times_vec(&temp[6], &temp_out[6], bd4x4_pm, 4);
  mat_times_vec(&temp[10], &temp_out[10], bd4x4_mp, 4);
  mat_times_vec(&temp[14], &temp_out[14], bd4x4_mm, 2);

  // output
  for (int i = 0; i < 16; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void gft_dct4x4_mat(const double *input, double *output) {
  mat_times_vec(input, output, dct4x4, 16);
#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void gft_dct4x4_btf(const double *input, double *output) {
  double temp[16] = {0.0};
  double temp1[16] = {0.0};
  double temp_out[16] = {0.0};

  int stage4_idx[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 12, 9, 13, 10, 11, 15, 14};
  int output_idx[16] = {0, 4, 5, 10, 1, 8, 6, 13, 2, 9, 7, 14, 3, 11, 15, 12};

  // stage 1
  temp[0] = input[0] + input[12];
  temp[12] = input[0] - input[12];
  temp[1] = input[1] + input[13];
  temp[13] = input[1] - input[13];
  temp[2] = input[2] + input[14];
  temp[14] = input[2] - input[14];
  temp[3] = input[3] + input[15];
  temp[15] = input[3] - input[15];
  temp[4] = input[4] + input[8];
  temp[8] = input[4] - input[8];
  temp[5] = input[5] + input[9];
  temp[9] = input[5] - input[9];
  temp[6] = input[6] + input[10];
  temp[10] = input[6] - input[10];
  temp[7] = input[7] + input[11];
  temp[11] = input[7] - input[11];

  // stage 2
  temp1[0] = temp[0] + temp[3];
  temp1[3] = temp[0] - temp[3];
  temp1[1] = temp[1] + temp[2];
  temp1[2] = temp[1] - temp[2];
  temp1[4] = temp[4] + temp[7];
  temp1[7] = temp[4] - temp[7];
  temp1[5] = temp[5] + temp[6];
  temp1[6] = temp[5] - temp[6];
  temp1[8] = temp[8] + temp[11];
  temp1[11] = temp[8] - temp[11];
  temp1[9] = temp[9] + temp[10];
  temp1[10] = temp[9] - temp[10];
  temp1[12] = temp[12] + temp[15];
  temp1[15] = temp[12] - temp[15];
  temp1[13] = temp[13] + temp[14];
  temp1[14] = temp[13] - temp[14];

  // stage 3
  temp[10] = SQRT2 * temp1[10];
  temp[15] = SQRT2 * temp1[15];
  temp[0] = temp1[0] + temp1[4];
  temp[4] = temp1[0] - temp1[4];
  temp[1] = temp1[1] + temp1[5];
  temp[5] = temp1[1] - temp1[5];
  temp[2] = temp1[2] + temp1[6];
  temp[6] = temp1[2] - temp1[6];
  temp[3] = temp1[3] + temp1[7];
  temp[7] = temp1[3] - temp1[7];
  temp[8] = temp1[8] + temp1[9];
  temp[9] = temp1[8] - temp1[9];
  temp[12] = temp1[12] + temp1[13];
  temp[13] = temp1[12] - temp1[13];
  temp[11] = temp1[11] + temp1[14];
  temp[14] = temp1[11] - temp1[14];

  // stage 4
  for (int i = 0; i < 16; i++)
    temp1[i] = temp[stage4_idx[i]];
  temp_out[0] = (temp1[0] + temp1[1]) / 4;
  temp_out[1] = (temp1[0] - temp1[1]) / 4;
  temp_out[2] = (temp1[2] + temp1[3]) / 4;
  temp_out[3] = (temp1[2] - temp1[3]) / 4;
  mat_times_vec(&temp1[4], &temp_out[4], dct4x4_pmp, 2);
  mat_times_vec(&temp1[6], &temp_out[6], dct4x4_pmm, 2);
  mat_times_vec(&temp1[8], &temp_out[8], dct4x4_pmp, 2);
  mat_times_vec(&temp1[10], &temp_out[10], dct4x4_pmm, 2);
  mat_times_vec(&temp1[12], &temp_out[12], dct4x4_mmp, 3);
  temp_out[15] = temp1[15] * INVSQRT2 / 2;

  // output
  for (int i = 0; i < 16; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void dct4_btf(const double *input, double *output) {
  double temp[4] = {0.0};
  double cos_factors[2] = {0.6532814824, 0.2705980501};

  // stage 1
  temp[0] = input[0] + input[3];
  temp[3] = input[0] - input[3];
  temp[1] = input[1] + input[2];
  temp[2] = input[1] - input[2];

  // stage 2
  output[0] = (temp[0] + temp[1]) / 2;
  output[2] = (temp[0] - temp[1]) / 2;
  output[1] = cos_factors[0] * temp[3] + cos_factors[1] * temp[2];
  output[3] = cos_factors[1] * temp[3] - cos_factors[0] * temp[2];
}

void gft_dct4x4_sep(const double *input, double *output) {
  double temp[16] = {0.0};
  double temp1[16] = {0.0};

  int flip_idx[16] = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
  int output_idx[16] = {0, 2, 5, 9, 1, 3, 7, 12, 4, 6, 10, 14, 8, 11, 13, 15};

  // column transform
  dct4_btf(input, temp);
  dct4_btf(&input[4], &temp[4]);
  dct4_btf(&input[8], &temp[8]);
  dct4_btf(&input[12], &temp[12]);

  for (int i = 0; i < 16; i++)
    temp1[i] = temp[flip_idx[i]];

  // row transform
  dct4_btf(temp1, temp);
  dct4_btf(&temp1[4], &temp[4]);
  dct4_btf(&temp1[8], &temp[8]);
  dct4_btf(&temp1[12], &temp[12]);

  // output (ordered by frequencies)
  for (int i = 0; i < 16; i++)
    output[output_idx[i]] = temp[i];
}

void gft_skeleton15_mat(const double *input, double *output) {
  mat_times_vec(input, output, sk15, 15);
}

void gft_skeleton15_btf(const double *input, double *output) {
  double temp[15] = {0.0};
  double temp_out[15] = {0.0};

  int output_idx[15] = {0, 1, 4, 5, 6, 9, 10, 13, 14, 2, 7, 11, 3, 8, 12};

  // stage 1
  temp[0] = SQRT2 * input[0];
  temp[1] = SQRT2 * input[1];
  temp[2] = SQRT2 * input[2];
  temp[3] = input[3] + input[6];
  temp[4] = input[4] + input[7];
  temp[5] = input[5] + input[8];
  temp[6] = input[9] + input[12];
  temp[7] = input[10] + input[13];
  temp[8] = input[11] + input[14];
  temp[9] = input[3] - input[6];
  temp[10] = input[4] - input[7];
  temp[11] = input[5] - input[8];
  temp[12] = input[9] - input[12];
  temp[13] = input[10] - input[13];
  temp[14] = input[11] - input[14];

  // stage 2
  mat_times_vec(temp, temp_out, sk15_p, 9);
  mat_times_vec(&temp[9], &temp_out[9], sk15_m, 3);
  mat_times_vec(&temp[12], &temp_out[12], sk15_m, 3);

  // output
  for (int i = 0; i < 15; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 15);
#endif
}
