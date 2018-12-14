#include "txfm.h"
#include "gftmatrices.h"
#if CONFIG_USE_LIFTING_FOR_GIVENS
#include "lifting_ptj.h"
#include "lifting_tj.h"
#else
#include "givens_ptj.h"
#include "givens_tj.h"
#endif

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n) {
  /*
    fprintf(stderr, "(");
    for (int i = 0; i < n - 1; i++)
      fprintf(stderr, "%.4lf,", input[i]);
    fprintf(stderr, "%.4lf)\n[", input[n - 1]);
    for (int i = 0; i < n - 1; i++)
      fprintf(stderr, "%.4lf,", output[i]);
    fprintf(stderr, "%.4lf]\n", output[n - 1]);
  */
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

#if CONFIG_USE_LIFTING_FOR_GIVENS
void serial_givens_lifting_tx(const double *input, double *output, int n,
                              const int *coords, const double *lifts,
                              const int *idx, int n_givens) {
  int ii, jj;
  double xi, xj;
  double temp[MAX_N] = {0.0};

  // copy data for in-place operations
  assert(n <= MAX_N);
  for (int i = 0; i < n; i++)
    temp[i] = input[i];

  for (int i = 0; i < n_givens; i++) {
    // i * 2 + j <==> [i][j]
    ii = coords[i * 2];
    jj = coords[i * 2 + 1];
    if (ii == -1 && jj == -1)
      break;
    // in-place lifting operations
    temp[ii] += lifts[i * 2] * temp[jj];
    temp[jj] += lifts[i * 2 + 1] * temp[ii];
    temp[ii] += lifts[i * 2] * temp[jj];
  }

  // order the frequency
  for (int i = 0; i < n; i++)
    output[idx[n_givens * n + i]] = temp[i];
}

void layered_givens_lifting_tx(const double *input, double *output, int n,
                               const int *coords, const double *lifts,
                               const int *idx, int n_layers) {
  int ii, jj;
  double xi, xj;
  double temp[MAX_N] = {0.0};

  // copy data for in-place operations
  assert(n < MAX_N);
  for (int i = 0; i < n; i++)
    temp[i] = input[i];

  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n / 2; j++) {
      // i * (n/2 * 2) + j * 2 + k  <==> [i][j][k]
      ii = coords[i * n + j * 2];
      jj = coords[i * n + j * 2 + 1];
      if (ii == -1 && jj == -1)
        continue;
      // in-place lifting operations
      temp[ii] += lifts[i * n + j * 2] * temp[jj];
      temp[jj] += lifts[i * n + j * 2 + 1] * temp[ii];
      temp[ii] += lifts[i * n + j * 2] * temp[jj];
    }
  }
  // order the frequency
  for (int i = 0; i < n; i++)
    output[idx[n_layers * n + i]] = temp[i];
}
#else
void serial_givens_tx(const double *input, double *output, int n,
                      const int *coords, const double *angles, const int *idx,
                      int n_givens) {
  int ii, jj;
  double xi, xj, anglecos, anglesin;
  double temp[MAX_N] = {0.0};

  // copy data for in-place operations
  assert(n <= MAX_N);
  for (int i = 0; i < n; i++)
    temp[i] = input[i];

  for (int i = 0; i < n_givens; i++) {
    // i * 2 + j <==> [i][j]
    ii = coords[i * 2];
    jj = coords[i * 2 + 1];
    if (ii == -1 && jj == -1)
      break;
    xi = temp[ii];
    xj = temp[jj];
    temp[ii] = xi * angles[i * 2] - xj * angles[i * 2 + 1];
    temp[jj] = xi * angles[i * 2 + 1] + xj * angles[i * 2];
  }

  // order the frequency
  for (int i = 0; i < n; i++)
    output[idx[n_givens * n + i]] = temp[i];
}

void layered_givens_tx(const double *input, double *output, int n,
                       const int *coords, const double *angles, const int *idx,
                       int n_layers) {
  int ii, jj;
  double xi, xj, anglecos, anglesin;
  double temp[MAX_N] = {0.0};

  // copy data for in-place operations
  assert(n < MAX_N);
  for (int i = 0; i < n; i++)
    temp[i] = input[i];

  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n / 2; j++) {
      // i * (n/2 * 2) + j * 2 + k  <==> [i][j][k]
      ii = coords[i * n + j * 2];
      jj = coords[i * n + j * 2 + 1];
      if (ii == -1 && jj == -1)
        continue;
      xi = temp[ii];
      xj = temp[jj];
      temp[ii] = xi * angles[i * n + j * 2] - xj * angles[i * n + j * 2 + 1];
      temp[jj] = xi * angles[i * n + j * 2 + 1] + xj * angles[i * n + j * 2];
    }
  }
  // order the frequency
  for (int i = 0; i < n; i++)
    output[idx[n_layers * n + i]] = temp[i];
}
#endif // CONFIG_USE_LIFTING_FOR_GIVENS

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

void gft_star100_mat(const double *input, double *output) {
  mat_times_vec(input, output, star100, 100);
#if CONFIG_DEBUG
  show_coefficients(input, output, 100);
#endif
}

void gft_star100_btf(const double *input, double *output) {
  double temp[51] = {0.0};
  double temp1[27] = {0.0};
  double temp2[15] = {0.0};
  double temp3[9] = {0.0};
  double temp4[6] = {0.0};
  double temp5[5] = {0.0};
  double temp_out[100] = {0.0};

  // stage 1 (2->50 51->99)
  temp[0] = SQRT2 * input[0];
  temp[1] = SQRT2 * input[1];
  for (int i = 2; i < 51; i++) {
    temp[i] = input[i] + input[101 - i];
    temp_out[101 - i] = (input[i] - input[101 - i]) * INVSQRT2;
  }

  // stage 2 (3->26 27->50)
  temp1[0] = SQRT2 * temp[0];
  temp1[1] = SQRT2 * temp[1];
  temp1[2] = SQRT2 * temp[2];
  for (int i = 3; i < 27; i++) {
    temp1[i] = temp[i] + temp[53 - i];
    temp_out[53 - i] = (temp[i] - temp[53 - i]) / 2;
  }

  // stage 3 (3->14 15->26)
  temp2[0] = SQRT2 * temp1[0];
  temp2[1] = SQRT2 * temp1[1];
  temp2[2] = SQRT2 * temp1[2];
  for (int i = 3; i < 15; i++) {
    temp2[i] = temp1[i] + temp1[29 - i];
    temp_out[29 - i] = (temp1[i] - temp1[29 - i]) / 2 * INVSQRT2;
  }

  // stage 4 (3->8 9->14)
  temp3[0] = SQRT2 * temp2[0];
  temp3[1] = SQRT2 * temp2[1];
  temp3[2] = SQRT2 * temp2[2];
  for (int i = 3; i < 9; i++) {
    temp3[i] = temp2[i] + temp2[17 - i];
    temp_out[17 - i] = (temp2[i] - temp2[17 - i]) / 4;
  }

  // stage 5 (3->5 6->8)
  temp4[0] = SQRT2 * temp3[0];
  temp4[1] = SQRT2 * temp3[1];
  temp4[2] = SQRT2 * temp3[2];
  for (int i = 3; i < 6; i++) {
    temp4[i] = temp3[i] + temp3[11 - i];
    temp_out[11 - i] = (temp3[i] - temp3[11 - i]) / 4 * INVSQRT2;
  }

  // stage 6 (4, 5)
  temp5[0] = SQRT2 * temp4[0];
  temp5[1] = SQRT2 * temp4[1];
  temp5[2] = SQRT2 * temp4[2];
  temp5[3] = SQRT2 * temp4[3];
  temp5[4] = temp4[4] + temp4[5];
  temp_out[5] = (temp4[4] - temp4[5]) / 8;

  // stage 7
  mat_times_vec(temp5, temp_out, star100_pppppp, 5);

  // output
  output[0] = temp_out[0];
  output[1] = temp_out[1];
  output[2] = temp_out[2];
  output[3] = temp_out[3];
  output[99] = temp_out[4];
  for (int i = 4; i < 99; i++)
    output[i] = temp_out[i + 1];

#if CONFIG_DEBUG
  show_coefficients(input, output, 10);
#endif
}

void gft_cycle12_mat(const double *input, double *output) {
  mat_times_vec(input, output, cycle12, 12);
#if CONFIG_DEBUG
  show_coefficients(input, output, 12);
#endif
}

void gft_cycle12_btf(const double *input, double *output) {
  double temp[12] = {0.0};
  double temp1[12] = {0.0};
  double temp_out[12] = {0.0};

  int inv1[2][6] = {{0, 1, 2, 3, 4, 5}, {11, 10, 9, 8, 7, 6}};
  int inv2[2][6] = {{0, 1, 2, 11, 10, 9}, {5, 4, 3, 6, 7, 8}};

  int output_idx[12] = {0, 8, 3, 1, 5, 9, 7, 4, 11, 2, 6, 10};

  // stage 1
  for (int i = 0; i < 6; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  for (int i = 0; i < 6; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // stage 3
  temp[1] = SQRT2 * temp1[1];
  temp[7] = SQRT2 * temp1[7];
  temp[0] = temp1[0] + temp1[2];
  temp[2] = temp1[0] - temp1[2];
  temp[8] = temp1[8] + temp1[6];
  temp[6] = temp1[8] - temp1[6];
  mat_times_vec(&temp1[3], &temp_out[3], cycle12_pm, 3);
  mat_times_vec(&temp1[9], &temp_out[9], cycle12_mp, 3);

  // stage 4
  mat_times_vec(temp, temp_out, cycle12_ppp, 2);
  temp_out[2] = temp[2] / 2 * INVSQRT2; // cycle12_ppm
  mat_times_vec(&temp[7], &temp_out[7], cycle12_mmp, 2);
  temp_out[6] = temp[6] / 2 * INVSQRT2; // cycle12_mmm

  // output
  for (int i = 0; i < 12; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 12);
#endif
}

void gft_cycle80_mat(const double *input, double *output) {
  mat_times_vec(input, output, cycle80, 80);
#if CONFIG_DEBUG
  show_coefficients(input, output, 80);
#endif
}

void gft_cycle80_btf(const double *input, double *output) {
  double temp[80] = {0.0};
  double temp1[80] = {0.0};
  double temp_out[80] = {0.0};

  int inv1[2][40] = {{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                      14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                      28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
                     {79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66,
                      65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52,
                      51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40}};
  int inv2[2][40] = {{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                      14, 15, 16, 17, 18, 19, 79, 78, 77, 76, 75, 74, 73, 72,
                      71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60},
                     {39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26,
                      25, 24, 23, 22, 21, 20, 40, 41, 42, 43, 44, 45, 46, 47,
                      48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59}};
  int inv3[2][20] = {
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50},
      {19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
       40, 41, 42, 43, 44, 45, 46, 47, 48, 49}};
  int inv4[2][10] = {{0, 1, 2, 3, 4, 49, 48, 47, 46, 45},
                     {9, 8, 7, 6, 5, 40, 41, 42, 43, 44}};

  int output_idx[80] = {0,  31, 63, 16, 47, 7,  23, 39, 55, 71, 3,  11, 19, 27,
                        35, 43, 51, 59, 67, 75, 1,  5,  9,  13, 17, 21, 25, 29,
                        33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 32, 64,
                        15, 48, 79, 8,  24, 40, 56, 72, 4,  12, 20, 28, 36, 44,
                        52, 60, 68, 76, 2,  6,  10, 14, 18, 22, 26, 30, 34, 38,
                        42, 46, 50, 54, 58, 62, 66, 70, 74, 78};

  // stage 1
  for (int i = 0; i < 40; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  for (int i = 0; i < 40; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // stage 3
  for (int i = 0; i < 20; i++) {
    temp[inv3[0][i]] = temp1[inv3[0][i]] + temp1[inv3[1][i]];
    temp[inv3[1][i]] = temp1[inv3[0][i]] - temp1[inv3[1][i]];
  }
  mat_times_vec(&temp1[20], &temp_out[20], cycle80_pm, 20);
  mat_times_vec(&temp1[60], &temp_out[60], cycle80_mp, 20);

  // stage 4
  for (int i = 0; i < 10; i++) {
    temp1[inv4[0][i]] = temp[inv4[0][i]] + temp[inv4[1][i]];
    temp1[inv4[1][i]] = temp[inv4[0][i]] - temp[inv4[1][i]];
  }
  mat_times_vec(&temp[10], &temp_out[10], cycle80_ppm, 10);
  mat_times_vec(&temp[50], &temp_out[50], cycle80_mmp, 10);

  // stage 5
  temp[2] = SQRT2 * temp1[2];
  temp[42] = SQRT2 * temp1[42];
  temp[0] = temp1[0] + temp1[4];
  temp[4] = temp1[0] - temp1[4];
  temp[1] = temp1[1] + temp1[3];
  temp[3] = temp1[1] - temp1[3];
  temp[44] = temp1[44] + temp1[40];
  temp[40] = temp1[44] - temp1[40];
  temp[43] = temp1[43] + temp1[41];
  temp[41] = temp1[43] - temp1[41];
  mat_times_vec(&temp1[5], &temp_out[5], cycle80_pppm, 5);
  mat_times_vec(&temp1[45], &temp_out[45], cycle80_mmmp, 5);

  // stage 6
  mat_times_vec(temp, temp_out, cycle80_ppppp, 3);
  mat_times_vec(&temp[3], &temp_out[3], cycle80_ppppm, 2);
  mat_times_vec(&temp[42], &temp_out[42], cycle80_mmmmp, 3);
  mat_times_vec(&temp[40], &temp_out[40], cycle80_mmmmm, 2);

  // output
  for (int i = 0; i < 80; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 80);
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

void gft_bd8x8_mat(const double *input, double *output) {
  mat_times_vec(input, output, bd8x8, 64);
#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
#endif
}

void gft_bd8x8_btf(const double *input, double *output) {
  double temp[64] = {0.0};
  double temp1[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][28] = {{1,  2,  3,  4,  5,  6,  7,  10, 11, 12, 13, 14, 15, 19,
                      20, 21, 22, 23, 28, 29, 30, 31, 37, 38, 39, 46, 47, 55},
                     {8,  16, 24, 32, 40, 48, 56, 17, 25, 33, 41, 49, 57, 26,
                      34, 42, 50, 58, 35, 43, 51, 59, 44, 52, 60, 53, 61, 62}};
  int inv2[2][28] = {{0,  1,  2,  3,  4,  5,  6,  8,  9,  10, 11, 12, 13, 16,
                      17, 18, 19, 20, 24, 25, 26, 27, 32, 33, 34, 40, 41, 48},
                     {63, 55, 47, 39, 31, 23, 15, 62, 54, 46, 38, 30, 22, 61,
                      53, 45, 37, 29, 60, 52, 44, 36, 59, 51, 43, 58, 50, 57}};

  int stage3_idx[64] = {0,  1,  2,  3,  4,  5,  6,  7,  9,  10, 11, 12, 13,
                        14, 18, 19, 20, 21, 27, 28, 15, 22, 23, 29, 30, 31,
                        36, 37, 38, 39, 45, 46, 47, 54, 55, 63, 8,  16, 17,
                        24, 25, 26, 32, 33, 34, 35, 40, 41, 42, 48, 49, 56,
                        43, 44, 50, 51, 52, 53, 57, 58, 59, 60, 61, 62};
  int output_idx[64] = {0,  3,  6,  8,  12, 15, 18, 22, 26, 29, 33, 37, 40,
                        43, 46, 51, 55, 57, 60, 63, 2,  7,  11, 14, 19, 23,
                        27, 30, 34, 39, 44, 47, 50, 54, 58, 62, 1,  4,  9,
                        13, 16, 21, 25, 28, 32, 36, 41, 45, 48, 53, 56, 61,
                        5,  10, 17, 20, 24, 31, 35, 38, 42, 49, 52, 59};

  // stage 1
  temp[0] = SQRT2 * input[0];
  temp[9] = SQRT2 * input[9];
  temp[18] = SQRT2 * input[18];
  temp[27] = SQRT2 * input[27];
  temp[36] = SQRT2 * input[36];
  temp[45] = SQRT2 * input[45];
  temp[54] = SQRT2 * input[54];
  temp[63] = SQRT2 * input[63];
  for (int i = 0; i < 28; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  temp1[7] = SQRT2 * temp[7];
  temp1[14] = SQRT2 * temp[14];
  temp1[21] = SQRT2 * temp[21];
  temp1[28] = SQRT2 * temp[28];
  temp1[35] = SQRT2 * temp[35];
  temp1[42] = SQRT2 * temp[42];
  temp1[49] = SQRT2 * temp[49];
  temp1[56] = SQRT2 * temp[56];
  for (int i = 0; i < 28; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // permutation
  for (int i = 0; i < 64; i++)
    temp[i] = temp1[stage3_idx[i]];

  // stage 3
  mat_times_vec(temp, temp_out, bd8x8_pp, 20);
  mat_times_vec(&temp[20], &temp_out[20], bd8x8_pm, 16);
  mat_times_vec(&temp[36], &temp_out[36], bd8x8_mp, 16);
  mat_times_vec(&temp[52], &temp_out[52], bd8x8_mm, 12);

  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
#endif
}

void gft_bd8x8_tj(const double *input, double *output, int n_givens) {
#if CONFIG_USE_LIFTING_FOR_GIVENS
  serial_givens_lifting_tx(input, output, 64, bd8x8_tj_coords, bd8x8_tj_lifts,
                           bd8x8_tj_idx, n_givens);
#else
  serial_givens_tx(input, output, 64, bd8x8_tj_coords, bd8x8_tj_angles,
                   bd8x8_tj_idx, n_givens);
#endif
}

void gft_bd8x8_ptj(const double *input, double *output, int n_layers) {
#if CONFIG_USE_LIFTING_FOR_GIVENS
  layered_givens_lifting_tx(input, output, 64, bd8x8_ptj_coords,
                            bd8x8_ptj_lifts, bd8x8_ptj_idx, n_layers);
#else
  layered_givens_tx(input, output, 64, bd8x8_ptj_coords, bd8x8_ptj_angles,
                    bd8x8_ptj_idx, n_layers);
#endif
}

void gft_bd8x8_btf_tj(const double *input, double *output, int n_givens) {
  double temp[64] = {0.0};
  double temp1[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][28] = {{1,  2,  3,  4,  5,  6,  7,  10, 11, 12, 13, 14, 15, 19,
                      20, 21, 22, 23, 28, 29, 30, 31, 37, 38, 39, 46, 47, 55},
                     {8,  16, 24, 32, 40, 48, 56, 17, 25, 33, 41, 49, 57, 26,
                      34, 42, 50, 58, 35, 43, 51, 59, 44, 52, 60, 53, 61, 62}};
  int inv2[2][28] = {{0,  1,  2,  3,  4,  5,  6,  8,  9,  10, 11, 12, 13, 16,
                      17, 18, 19, 20, 24, 25, 26, 27, 32, 33, 34, 40, 41, 48},
                     {63, 55, 47, 39, 31, 23, 15, 62, 54, 46, 38, 30, 22, 61,
                      53, 45, 37, 29, 60, 52, 44, 36, 59, 51, 43, 58, 50, 57}};

  int stage3_idx[64] = {0,  1,  2,  3,  4,  5,  6,  7,  9,  10, 11, 12, 13,
                        14, 18, 19, 20, 21, 27, 28, 15, 22, 23, 29, 30, 31,
                        36, 37, 38, 39, 45, 46, 47, 54, 55, 63, 8,  16, 17,
                        24, 25, 26, 32, 33, 34, 35, 40, 41, 42, 48, 49, 56,
                        43, 44, 50, 51, 52, 53, 57, 58, 59, 60, 61, 62};
  int output_idx[64] = {0,  3,  6,  8,  12, 15, 18, 22, 26, 29, 33, 37, 40,
                        43, 46, 51, 55, 57, 60, 63, 2,  7,  11, 14, 19, 23,
                        27, 30, 34, 39, 44, 47, 50, 54, 58, 62, 1,  4,  9,
                        13, 16, 21, 25, 28, 32, 36, 41, 45, 48, 53, 56, 61,
                        5,  10, 17, 20, 24, 31, 35, 38, 42, 49, 52, 59};

  // stage 1
  temp[0] = SQRT2 * input[0];
  temp[9] = SQRT2 * input[9];
  temp[18] = SQRT2 * input[18];
  temp[27] = SQRT2 * input[27];
  temp[36] = SQRT2 * input[36];
  temp[45] = SQRT2 * input[45];
  temp[54] = SQRT2 * input[54];
  temp[63] = SQRT2 * input[63];
  for (int i = 0; i < 28; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  temp1[7] = SQRT2 * temp[7];
  temp1[14] = SQRT2 * temp[14];
  temp1[21] = SQRT2 * temp[21];
  temp1[28] = SQRT2 * temp[28];
  temp1[35] = SQRT2 * temp[35];
  temp1[42] = SQRT2 * temp[42];
  temp1[49] = SQRT2 * temp[49];
  temp1[56] = SQRT2 * temp[56];
  for (int i = 0; i < 28; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // permutation
  for (int i = 0; i < 64; i++)
    temp[i] = temp1[stage3_idx[i]];

  // stage 3
#if CONFIG_USE_LIFTING_FOR_GIVENS
  serial_givens_lifting_tx(temp, temp_out, 20, bd8x8_pp_tj_coords,
                           bd8x8_pp_tj_lifts, bd8x8_pp_tj_idx, n_givens / 4);
  serial_givens_lifting_tx(&temp[20], &temp_out[20], 16, bd8x8_pm_tj_coords,
                           bd8x8_pm_tj_lifts, bd8x8_pm_tj_idx, n_givens / 4);
  serial_givens_lifting_tx(&temp[36], &temp_out[36], 16, bd8x8_mp_tj_coords,
                           bd8x8_mp_tj_lifts, bd8x8_mp_tj_idx, n_givens / 4);
  serial_givens_lifting_tx(&temp[52], &temp_out[52], 12, bd8x8_mm_tj_coords,
                           bd8x8_mm_tj_lifts, bd8x8_mm_tj_idx, n_givens / 4);
#else
  serial_givens_tx(temp, temp_out, 20, bd8x8_pp_tj_coords, bd8x8_pp_tj_angles,
                   bd8x8_pp_tj_idx, n_givens / 4);
  serial_givens_tx(&temp[20], &temp_out[20], 16, bd8x8_pm_tj_coords,
                   bd8x8_pm_tj_angles, bd8x8_pm_tj_idx, n_givens / 4);
  serial_givens_tx(&temp[36], &temp_out[36], 16, bd8x8_mp_tj_coords,
                   bd8x8_mp_tj_angles, bd8x8_mp_tj_idx, n_givens / 4);
  serial_givens_tx(&temp[52], &temp_out[52], 12, bd8x8_mm_tj_coords,
                   bd8x8_mm_tj_angles, bd8x8_mm_tj_idx, n_givens / 4);
#endif
  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i] / 2;

#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
#endif
}

void gft_bd8x8_btf_ptj(const double *input, double *output, int n_layers) {
  double temp[64] = {0.0};
  double temp1[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][28] = {{1,  2,  3,  4,  5,  6,  7,  10, 11, 12, 13, 14, 15, 19,
                      20, 21, 22, 23, 28, 29, 30, 31, 37, 38, 39, 46, 47, 55},
                     {8,  16, 24, 32, 40, 48, 56, 17, 25, 33, 41, 49, 57, 26,
                      34, 42, 50, 58, 35, 43, 51, 59, 44, 52, 60, 53, 61, 62}};
  int inv2[2][28] = {{0,  1,  2,  3,  4,  5,  6,  8,  9,  10, 11, 12, 13, 16,
                      17, 18, 19, 20, 24, 25, 26, 27, 32, 33, 34, 40, 41, 48},
                     {63, 55, 47, 39, 31, 23, 15, 62, 54, 46, 38, 30, 22, 61,
                      53, 45, 37, 29, 60, 52, 44, 36, 59, 51, 43, 58, 50, 57}};

  int stage3_idx[64] = {0,  1,  2,  3,  4,  5,  6,  7,  9,  10, 11, 12, 13,
                        14, 18, 19, 20, 21, 27, 28, 15, 22, 23, 29, 30, 31,
                        36, 37, 38, 39, 45, 46, 47, 54, 55, 63, 8,  16, 17,
                        24, 25, 26, 32, 33, 34, 35, 40, 41, 42, 48, 49, 56,
                        43, 44, 50, 51, 52, 53, 57, 58, 59, 60, 61, 62};
  int output_idx[64] = {0,  3,  6,  8,  12, 15, 18, 22, 26, 29, 33, 37, 40,
                        43, 46, 51, 55, 57, 60, 63, 2,  7,  11, 14, 19, 23,
                        27, 30, 34, 39, 44, 47, 50, 54, 58, 62, 1,  4,  9,
                        13, 16, 21, 25, 28, 32, 36, 41, 45, 48, 53, 56, 61,
                        5,  10, 17, 20, 24, 31, 35, 38, 42, 49, 52, 59};

  // stage 1
  temp[0] = SQRT2 * input[0];
  temp[9] = SQRT2 * input[9];
  temp[18] = SQRT2 * input[18];
  temp[27] = SQRT2 * input[27];
  temp[36] = SQRT2 * input[36];
  temp[45] = SQRT2 * input[45];
  temp[54] = SQRT2 * input[54];
  temp[63] = SQRT2 * input[63];
  for (int i = 0; i < 28; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  temp1[7] = SQRT2 * temp[7];
  temp1[14] = SQRT2 * temp[14];
  temp1[21] = SQRT2 * temp[21];
  temp1[28] = SQRT2 * temp[28];
  temp1[35] = SQRT2 * temp[35];
  temp1[42] = SQRT2 * temp[42];
  temp1[49] = SQRT2 * temp[49];
  temp1[56] = SQRT2 * temp[56];
  for (int i = 0; i < 28; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // permutation
  for (int i = 0; i < 64; i++)
    temp[i] = temp1[stage3_idx[i]];

  // stage 3
#if CONFIG_USE_LIFTING_FOR_GIVENS
  layered_givens_lifting_tx(temp, temp_out, 20, bd8x8_pp_ptj_coords,
                            bd8x8_pp_ptj_lifts, bd8x8_pp_ptj_idx, n_layers);
  layered_givens_lifting_tx(&temp[20], &temp_out[20], 16, bd8x8_pm_ptj_coords,
                            bd8x8_pm_ptj_lifts, bd8x8_pm_ptj_idx, n_layers);
  layered_givens_lifting_tx(&temp[36], &temp_out[36], 16, bd8x8_mp_ptj_coords,
                            bd8x8_mp_ptj_lifts, bd8x8_mp_ptj_idx, n_layers);
  layered_givens_lifting_tx(&temp[52], &temp_out[52], 12, bd8x8_mm_ptj_coords,
                            bd8x8_mm_ptj_lifts, bd8x8_mm_ptj_idx, n_layers);
#else
  layered_givens_tx(temp, temp_out, 20, bd8x8_pp_ptj_coords,
                    bd8x8_pp_ptj_angles, bd8x8_pp_ptj_idx, n_layers);
  layered_givens_tx(&temp[20], &temp_out[20], 16, bd8x8_pm_ptj_coords,
                    bd8x8_pm_ptj_angles, bd8x8_pm_ptj_idx, n_layers);
  layered_givens_tx(&temp[36], &temp_out[36], 16, bd8x8_mp_ptj_coords,
                    bd8x8_mp_ptj_angles, bd8x8_mp_ptj_idx, n_layers);
  layered_givens_tx(&temp[52], &temp_out[52], 12, bd8x8_mm_ptj_coords,
                    bd8x8_mm_ptj_angles, bd8x8_mm_ptj_idx, n_layers);
#endif
  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i] / 2;

#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
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
  mat_times_vec(&temp1[6], &temp_out[6], dct4x4_pmp, 2);
  mat_times_vec(&temp1[8], &temp_out[8], dct4x4_pmp, 2);
  mat_times_vec(&temp1[10], &temp_out[10], dct4x4_pmp, 2);
  mat_times_vec(&temp1[12], &temp_out[12], dct4x4_mmp, 3);
  temp_out[15] = temp1[15] * INVSQRT2 / 2;

  // output
  for (int i = 0; i < 16; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void gft_dct8x8_mat(const double *input, double *output) {
  mat_times_vec(input, output, dct8x8, 64);
#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
#endif
}

void gft_dct8x8_btf(const double *input, double *output) {
  double temp[64] = {0.0};
  double temp1[64] = {0.0};
  double temp2[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][32] = {
      {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
      {56, 57, 58, 59, 60, 61, 62, 63, 48, 49, 50, 51, 52, 53, 54, 55,
       40, 41, 42, 43, 44, 45, 46, 47, 32, 33, 34, 35, 36, 37, 38, 39}};
  int inv2[2][32] = {
      {0,  1,  2,  3,  8,  9,  10, 11, 16, 17, 18, 19, 24, 25, 26, 27,
       32, 33, 34, 35, 40, 41, 42, 43, 48, 49, 50, 51, 56, 57, 58, 59},
      {7,  6,  5,  4,  15, 14, 13, 12, 23, 22, 21, 20, 31, 30, 29, 28,
       39, 38, 37, 36, 47, 46, 45, 44, 55, 54, 53, 52, 63, 62, 61, 60}};
  int inv3[2][30] = {
      {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14,
       15, 32, 40, 48, 56, 33, 41, 49, 57, 37, 38, 39, 46, 47, 55},
      {24, 25, 26, 27, 28, 29, 30, 31, 16, 17, 18, 19, 20, 21, 22,
       23, 35, 43, 51, 59, 34, 42, 50, 58, 44, 52, 60, 53, 61, 62}};
  int inv4[2][16] = {
      {0, 1, 8, 9, 16, 17, 24, 25, 4, 5, 6, 7, 32, 40, 48, 56},
      {3, 2, 11, 10, 19, 18, 27, 26, 12, 13, 14, 15, 33, 41, 49, 57}};
  int inv5[2][7] = {{0, 1, 2, 3, 16, 24, 19}, {8, 9, 10, 11, 17, 25, 26}};

  int stage6_idx1[20] = {0,  1,  8,  9,  4,  5,  6,  7,  12, 13,
                         14, 15, 32, 40, 48, 56, 33, 41, 49, 57};
  int stage6_idx[44] = {2,  3,  10, 11, 16, 24, 17, 25, 18, 19, 27,
                        26, 20, 21, 22, 23, 28, 29, 30, 31, 34, 42,
                        50, 58, 35, 43, 51, 59, 36, 37, 38, 39, 45,
                        46, 47, 54, 55, 63, 44, 52, 53, 60, 61, 62};
  int output_idx[64] = {0,  15,  16,  37,  1, 9, 22, 34, 18, 27, 47, 54, 2,
                        10, 23, 35, 17, 26, 48, 55, 4,  30,  21,  51, 5,  31,
                        20,  52, 8,  36, 60, 38, 7, 13, 28, 33, 43, 46, 57,
                        61, 6, 14, 29, 32, 44, 45, 56, 62, 3, 12, 19, 25,
                        39, 41, 49, 53, 58, 63, 11, 24, 40, 42, 50, 59};

  // stage 1 (left-right symmetry)
  for (int i = 0; i < 32; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2 (up-down symmetry)
  for (int i = 0; i < 32; i++) {
    temp1[inv2[0][i]] = temp[inv2[0][i]] + temp[inv2[1][i]];
    temp1[inv2[1][i]] = temp[inv2[0][i]] - temp[inv2[1][i]];
  }

  // stage 3 (L-R sym for Gpp and Gpm, U-D sym for Gmp, diagonal sym for Gmm)
  temp[36] = SQRT2 * temp1[36];
  temp[45] = SQRT2 * temp1[45];
  temp[54] = SQRT2 * temp1[54];
  temp[63] = SQRT2 * temp1[63];
  for (int i = 0; i < 30; i++) {
    temp[inv3[0][i]] = temp1[inv3[0][i]] + temp1[inv3[1][i]];
    temp[inv3[1][i]] = temp1[inv3[0][i]] - temp1[inv3[1][i]];
  }

  // stage 4 (U-D sym for Gppp, Gppm, and Gmpp, L-R sym for Gpmp)
  for (int i = 0; i < 16; i++) {
    temp1[inv4[0][i]] = temp[inv4[0][i]] + temp[inv4[1][i]];
    temp1[inv4[1][i]] = temp[inv4[0][i]] - temp[inv4[1][i]];
  }

  // stage 5 (L-R sym for Gpppp and Gpppm, U-D sym for Gppmp, and diag sym
  // for Gppmm)
  temp[18] = SQRT2 * temp1[18];
  temp[27] = SQRT2 * temp1[27];
  for (int i = 0; i < 7; i++) {
    temp[inv5[0][i]] = temp1[inv5[0][i]] + temp1[inv5[1][i]];
    temp[inv5[1][i]] = temp1[inv5[0][i]] - temp1[inv5[1][i]];
  }

  // stage 6 (U-D sym for Gppppp)
  temp1[0] = temp[0] + temp[1];
  temp1[1] = temp[0] - temp[1];
  temp1[8] = temp[8] + temp[9];
  temp1[9] = temp[8] - temp[9];

  // permutation
  for (int i = 0; i < 20; i++)
    temp2[i] = temp1[stage6_idx1[i]];
  for (int i = 0; i < 44; i++)
    temp2[20 + i] = temp[stage6_idx[i]];

  // stage 6 (U-D sym for Gppppp, and matrices for others)
  temp_out[0] = temp2[0] / 8;
  temp_out[1] = temp2[1] / 8;
  temp_out[2] = temp2[2] / 8;
  temp_out[3] = temp2[3] / 8;
  mat_times_vec(&temp2[4], &temp_out[4], dct8x8_pmpp, 4);
  mat_times_vec(&temp2[8], &temp_out[8], dct8x8_pmpp, 4);
  mat_times_vec(&temp2[12], &temp_out[12], dct8x8_pmpp, 4);
  mat_times_vec(&temp2[16], &temp_out[16], dct8x8_pmpp, 4);

  mat_times_vec(&temp2[20], &temp_out[20], dct8x8_pppmp, 2);
  mat_times_vec(&temp2[22], &temp_out[22], dct8x8_pppmp, 2);
  mat_times_vec(&temp2[24], &temp_out[24], dct8x8_pppmp, 2);
  mat_times_vec(&temp2[26], &temp_out[26], dct8x8_pppmp, 2);
  mat_times_vec(&temp2[28], &temp_out[28], dct8x8_ppmmp, 3);
  temp_out[31] = temp2[31] / 4 * INVSQRT2;
  mat_times_vec(&temp2[32], &temp_out[32], dct8x8_pmm, 8);
  mat_times_vec(&temp2[40], &temp_out[40], dct8x8_pmm, 8);
  mat_times_vec(&temp2[48], &temp_out[48], dct8x8_mmp, 10);
  mat_times_vec(&temp2[58], &temp_out[58], dct8x8_mmm, 6);

  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 64);
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

void dct8_btf(const double *input, double *output) {
  double temp[8] = {0.0};
  double temp1[8] = {0.0};
  // cos((0:7)*pi/16) / 2
  double cos_factors[8] = {0.5000000000, 0.4903926402, 0.4619397663,
                           0.4157348062, 0.3535533906, 0.2777851165,
                           0.1913417162, 0.0975451610};

  // stage 1
  temp[0] = input[0] + input[7];
  temp[7] = input[0] - input[7];
  temp[1] = input[1] + input[6];
  temp[6] = input[1] - input[6];
  temp[2] = input[2] + input[5];
  temp[5] = input[2] - input[5];
  temp[3] = input[3] + input[4];
  temp[4] = input[3] - input[4];

  // stage 2
  temp1[0] = temp[0] + temp[3];
  temp1[3] = temp[0] - temp[3];
  temp1[1] = temp[1] + temp[2];
  temp1[2] = temp[1] - temp[2];
  temp1[4] = temp[4];
  temp1[5] = (-temp[5] + temp[6]) * INVSQRT2;
  temp1[6] = (temp[5] + temp[6]) * INVSQRT2;
  temp1[7] = temp[7];

  // stage 3
  temp[4] = temp1[4] + temp1[5];
  temp[5] = temp1[4] - temp1[5];
  temp[6] = -temp1[6] + temp1[7];
  temp[7] = temp1[6] + temp1[7];

  // stage 4
  output[0] = (temp1[0] + temp1[1]) * cos_factors[4];
  output[4] = (temp1[0] - temp1[1]) * cos_factors[4];
  output[2] = temp1[2] * cos_factors[6] + temp1[3] * cos_factors[2];
  output[6] = -temp1[2] * cos_factors[2] + temp1[3] * cos_factors[6];
  output[1] = temp[4] * cos_factors[7] + temp[7] * cos_factors[1];
  output[7] = -temp[4] * cos_factors[1] + temp[7] * cos_factors[7];
  output[5] = temp[5] * cos_factors[3] + temp[6] * cos_factors[5];
  output[3] = -temp[5] * cos_factors[5] + temp[6] * cos_factors[3];
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

void gft_dct8x8_sep(const double *input, double *output) {
  double temp[64] = {0.0};
  double temp1[64] = {0.0};

  int flip_idx[64] = {0,  8,  16, 24, 32, 40, 48, 56, 1,  9,  17, 25, 33,
                      41, 49, 57, 2,  10, 18, 26, 34, 42, 50, 58, 3,  11,
                      19, 27, 35, 43, 51, 59, 4,  12, 20, 28, 36, 44, 52,
                      60, 5,  13, 21, 29, 37, 45, 53, 61, 6,  14, 22, 30,
                      38, 46, 54, 62, 7,  15, 23, 31, 39, 47, 55, 63};
  int output_idx[64] = {0,  2,  5,  10, 16, 23, 31, 35, 1,  3,  7,  12, 18,
                        25, 33, 42, 4,  6,  8,  14, 21, 29, 41, 44, 9,  11,
                        13, 19, 27, 38, 46, 50, 15, 17, 20, 26, 37, 48, 52,
                        55, 22, 24, 28, 36, 47, 53, 57, 59, 30, 32, 40, 45,
                        51, 56, 60, 62, 34, 39, 43, 49, 54, 58, 61, 63};

  // column transform
  dct8_btf(input, temp);
  dct8_btf(&input[8], &temp[8]);
  dct8_btf(&input[16], &temp[16]);
  dct8_btf(&input[24], &temp[24]);
  dct8_btf(&input[32], &temp[32]);
  dct8_btf(&input[40], &temp[40]);
  dct8_btf(&input[48], &temp[48]);
  dct8_btf(&input[56], &temp[56]);

  for (int i = 0; i < 64; i++)
    temp1[i] = temp[flip_idx[i]];

  // row transform
  dct8_btf(temp1, temp);
  dct8_btf(&temp1[8], &temp[8]);
  dct8_btf(&temp1[16], &temp[16]);
  dct8_btf(&temp1[24], &temp[24]);
  dct8_btf(&temp1[32], &temp[32]);
  dct8_btf(&temp1[40], &temp[40]);
  dct8_btf(&temp1[48], &temp[48]);
  dct8_btf(&temp1[56], &temp[56]);

  // output (ordered by frequencies)
  for (int i = 0; i < 64; i++)
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

void gft_skeleton25_mat(const double *input, double *output) {
  mat_times_vec(input, output, sk25, 25);
}

void gft_skeleton25_btf(const double *input, double *output) {
  double temp[25] = {0.0};
  double temp_out[25] = {0.0};

  int output_idx[25] = {0,  1,  4, 5, 7,  8,  11, 12, 13, 16, 17, 19, 21,
                        23, 24, 2, 6, 10, 14, 18, 22, 3,  9,  15, 20};

  // stage 1 & permutation
  temp[0] = SQRT2 * input[0];
  temp[1] = SQRT2 * input[1];
  temp[2] = SQRT2 * input[2];
  temp[3] = SQRT2 * input[3];
  temp[12] = SQRT2 * input[20];
  temp[4] = input[8] + input[4];
  temp[5] = input[9] + input[5];
  temp[6] = input[10] + input[6];
  temp[7] = input[11] + input[7];
  temp[8] = input[16] + input[12];
  temp[9] = input[17] + input[13];
  temp[10] = input[18] + input[14];
  temp[11] = input[19] + input[15];
  temp[13] = input[23] + input[21];
  temp[14] = input[24] + input[22];
  temp[15] = input[8] - input[4];
  temp[16] = input[9] - input[5];
  temp[17] = input[10] - input[6];
  temp[18] = input[11] - input[7];
  temp[19] = input[23] - input[21];
  temp[20] = input[24] - input[22];
  temp[21] = input[16] - input[12];
  temp[22] = input[17] - input[13];
  temp[23] = input[18] - input[14];
  temp[24] = input[19] - input[15];

  // stage 2
  mat_times_vec(temp, temp_out, sk25_p, 15);
  mat_times_vec(&temp[15], &temp_out[15], sk25_m1, 6);
  mat_times_vec(&temp[21], &temp_out[21], sk25_m2, 4);

  // output
  for (int i = 0; i < 25; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 25);
#endif
}

void gft_z4x4_mat(const double *input, double *output) {
  mat_times_vec(input, output, z4x4, 16);
}

void gft_z4x4_btf(const double *input, double *output) {
  double temp[16] = {0.0};
  double temp_out[16] = {0.0};

  int output_idx[16] = {0, 2, 4, 7, 8, 11, 12, 14, 1, 3, 5, 6, 9, 10, 13, 15};

  // stage 1
  temp[0] = input[0] + input[15];
  temp[15] = input[0] - input[15];
  temp[1] = input[1] + input[14];
  temp[14] = input[1] - input[14];
  temp[2] = input[2] + input[13];
  temp[13] = input[2] - input[13];
  temp[3] = input[3] + input[12];
  temp[12] = input[3] - input[12];
  temp[4] = input[4] + input[11];
  temp[11] = input[4] - input[11];
  temp[5] = input[5] + input[10];
  temp[10] = input[5] - input[10];
  temp[6] = input[6] + input[9];
  temp[9] = input[6] - input[9];
  temp[7] = input[7] + input[8];
  temp[8] = input[7] - input[8];

  // stage 2
  mat_times_vec(temp, temp_out, z4x4_p, 8);
  mat_times_vec(&temp[8], &temp_out[8], z4x4_m, 8);

  // output
  for (int i = 0; i < 16; i++)
    output[output_idx[i]] = temp_out[i];

#if CONFIG_DEBUG
  show_coefficients(input, output, 16);
#endif
}

void gft_z8x8_mat(const double *input, double *output) {
  mat_times_vec(input, output, z8x8, 64);
}

void gft_z8x8_btf(const double *input, double *output) {
  double temp[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][32] = {
      {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
      {63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48,
       47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32}};
  int output_idx[64] = {0,  2,  5,  7,  8,  11, 12, 14, 16, 18, 20, 22, 24,
                        26, 29, 31, 32, 34, 37, 38, 41, 42, 44, 47, 48, 49,
                        53, 55, 56, 59, 60, 62, 1,  3,  4,  6,  9,  10, 13,
                        15, 17, 19, 21, 23, 25, 27, 28, 30, 33, 35, 36, 39,
                        40, 43, 45, 46, 50, 51, 52, 54, 57, 58, 61, 63};

  // stage 1
  for (int i = 0; i < 32; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
  mat_times_vec(temp, temp_out, z8x8_p, 32);
  mat_times_vec(&temp[32], &temp_out[32], z8x8_m, 32);

  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i] * INVSQRT2;
}

void gft_z8x8_tj(const double *input, double *output, int n_givens) {
#if CONFIG_USE_LIFTING_FOR_GIVENS
  serial_givens_lifting_tx(input, output, 64, z8x8_tj_coords, z8x8_tj_lifts,
                           z8x8_tj_idx, n_givens);
#else
  serial_givens_tx(input, output, 64, z8x8_tj_coords, z8x8_tj_angles,
                   z8x8_tj_idx, n_givens);
#endif
}

void gft_z8x8_ptj(const double *input, double *output, int n_layers) {
#if CONFIG_USE_LIFTING_FOR_GIVENS
  layered_givens_lifting_tx(input, output, 64, z8x8_ptj_coords, z8x8_ptj_lifts,
                            z8x8_ptj_idx, n_layers);
#else
  layered_givens_tx(input, output, 64, z8x8_ptj_coords, z8x8_ptj_angles,
                    z8x8_ptj_idx, n_layers);
#endif
}

void gft_z8x8_btf_tj(const double *input, double *output, int n_givens) {
  double temp[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][32] = {
      {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
      {63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48,
       47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32}};
  int output_idx[64] = {0,  2,  5,  7,  8,  11, 12, 14, 16, 18, 20, 22, 24,
                        26, 29, 31, 32, 34, 37, 38, 41, 42, 44, 47, 48, 49,
                        53, 55, 56, 59, 60, 62, 1,  3,  4,  6,  9,  10, 13,
                        15, 17, 19, 21, 23, 25, 27, 28, 30, 33, 35, 36, 39,
                        40, 43, 45, 46, 50, 51, 52, 54, 57, 58, 61, 63};

  // stage 1
  for (int i = 0; i < 32; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
#if CONFIG_USE_LIFTING_FOR_GIVENS
  serial_givens_lifting_tx(temp, temp_out, 32, z8x8_p_tj_coords,
                           z8x8_p_tj_lifts, z8x8_p_tj_idx, n_givens / 2);
  serial_givens_lifting_tx(&temp[32], &temp_out[32], 32, z8x8_m_tj_coords,
                           z8x8_m_tj_lifts, z8x8_m_tj_idx, n_givens / 2);
#else
  serial_givens_tx(temp, temp_out, 32, z8x8_p_tj_coords, z8x8_p_tj_angles,
                   z8x8_p_tj_idx, n_givens / 2);
  serial_givens_tx(&temp[32], &temp_out[32], 32, z8x8_m_tj_coords,
                   z8x8_m_tj_angles, z8x8_m_tj_idx, n_givens / 2);
#endif
  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i] * INVSQRT2;
}

void gft_z8x8_btf_ptj(const double *input, double *output, int n_layers) {
  double temp[64] = {0.0};
  double temp_out[64] = {0.0};

  int inv1[2][32] = {
      {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
      {63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48,
       47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32}};
  int output_idx[64] = {0,  2,  5,  7,  8,  11, 12, 14, 16, 18, 20, 22, 24,
                        26, 29, 31, 32, 34, 37, 38, 41, 42, 44, 47, 48, 49,
                        53, 55, 56, 59, 60, 62, 1,  3,  4,  6,  9,  10, 13,
                        15, 17, 19, 21, 23, 25, 27, 28, 30, 33, 35, 36, 39,
                        40, 43, 45, 46, 50, 51, 52, 54, 57, 58, 61, 63};

  // stage 1
  for (int i = 0; i < 32; i++) {
    temp[inv1[0][i]] = input[inv1[0][i]] + input[inv1[1][i]];
    temp[inv1[1][i]] = input[inv1[0][i]] - input[inv1[1][i]];
  }

  // stage 2
#if CONFIG_USE_LIFTING_FOR_GIVENS
  layered_givens_lifting_tx(temp, temp_out, 32, z8x8_p_ptj_coords,
                            z8x8_p_ptj_lifts, z8x8_p_ptj_idx, n_layers);
  layered_givens_lifting_tx(&temp[32], &temp_out[32], 32, z8x8_m_ptj_coords,
                            z8x8_m_ptj_lifts, z8x8_m_ptj_idx, n_layers);
#else
  layered_givens_tx(temp, temp_out, 32, z8x8_p_ptj_coords, z8x8_p_ptj_angles,
                    z8x8_p_ptj_idx, n_layers);
  layered_givens_tx(&temp[32], &temp_out[32], 32, z8x8_m_ptj_coords,
                    z8x8_m_ptj_angles, z8x8_m_ptj_idx, n_layers);
#endif
  // output
  for (int i = 0; i < 64; i++)
    output[output_idx[i]] = temp_out[i] * INVSQRT2;
}
