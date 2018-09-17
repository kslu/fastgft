#include "txfm.h"

#define NUM_TJ_GIVENS_TO_RUN 21
#define NUM_BTFTJ_GIVENS_TO_RUN 23
#define NUM_PTJ_LAYERS_TO_RUN 25
#define NUM_BTFPTJ_LAYERS_TO_RUN 21

int main(int argc, char *argv[]) {

  int n = 64;
  int n_inputs;
  int givens_tj[NUM_TJ_GIVENS_TO_RUN] = {
      0,    200,  400,  600,  800,  1000, 1200, 1400, 1600, 1800, 2000,
      2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000};
  int givens_btftj[NUM_BTFTJ_GIVENS_TO_RUN] = {
      0,   100, 200,  300,  400,  500,  600,  700,
      800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
			1900, 2000, 2100, 2200};
  int givens_layers_ptj[NUM_PTJ_LAYERS_TO_RUN] = {
      0,  5,  10, 15, 20, 25, 30, 35,  40,  45,  50,  55, 60,
      65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120};
  int givens_layers_btfptj[NUM_BTFPTJ_LAYERS_TO_RUN] = {
      0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85,
			90, 95, 100};
  double buffer_in[BATCH_SIZE][64];
  double buffer_out_mat[BATCH_SIZE][64];
  double buffer_out_btf[BATCH_SIZE][64];
  double buffer_out_tj[BATCH_SIZE][64];
  double buffer_out_ptj[BATCH_SIZE][64];
  double buffer_out_btftj[BATCH_SIZE][64];
  double buffer_out_btfptj[BATCH_SIZE][64];

  double adiff = 0;
  double acc_error_tj[NUM_TJ_GIVENS_TO_RUN] = {0};
  double acc_error_ptj[NUM_PTJ_LAYERS_TO_RUN] = {0};
  double acc_error_btftj[NUM_BTFTJ_GIVENS_TO_RUN] = {0};
  double acc_error_btfptj[NUM_BTFPTJ_LAYERS_TO_RUN] = {0};

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_mat = 0, t_btf = 0, t_temp = 0;
  clock_t t_tj[NUM_TJ_GIVENS_TO_RUN] = {0};
  clock_t t_ptj[NUM_PTJ_LAYERS_TO_RUN] = {0};
  clock_t t_btftj[NUM_BTFTJ_GIVENS_TO_RUN] = {0};
  clock_t t_btfptj[NUM_BTFPTJ_LAYERS_TO_RUN] = {0};

  for (int b = 0; b < n_batches; b++) {
    fprintf(stderr, "Batch %d/%d ...\n", b + 1, n_batches);
    cur_batch_size = b < n_batches - 1 ? BATCH_SIZE : n_inputs - b * BATCH_SIZE;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < n; j++)
        fscanf(fp_in, "%lf", &buffer_in[i][j]);
    }

    // matrix GFT
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      gft_z8x8_mat(buffer_in[i], buffer_out_mat[i]);
    t_mat += clock() - t_temp;

    // butterfly GFT
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      gft_z8x8_btf(buffer_in[i], buffer_out_btf[i]);
    t_btf += clock() - t_temp;

    // Truncated Jacobi (TJ)-GFT
    for (int k = 0; k < NUM_TJ_GIVENS_TO_RUN; k++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        gft_z8x8_tj(buffer_in[i], buffer_out_tj[i], givens_tj[k]);
      t_tj[k] += clock() - t_temp;
      // compute GFT coefficients error
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < n; j++) {
          adiff = fabs(buffer_out_mat[i][j]) - fabs(buffer_out_tj[i][j]);
          acc_error_tj[k] += adiff * adiff;
        }
      }
    }

    // Parallel Truncated Jacobi (PTJ)-GFT
    for (int k = 0; k < NUM_PTJ_LAYERS_TO_RUN; k++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        gft_z8x8_ptj(buffer_in[i], buffer_out_ptj[i], givens_layers_ptj[k]);
      t_ptj[k] += clock() - t_temp;
      // compute GFT coefficients error
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < n; j++) {
          adiff = fabs(buffer_out_mat[i][j]) - fabs(buffer_out_ptj[i][j]);
          acc_error_ptj[k] += adiff * adiff;
        }
      }
    }

    // butterfly + TJ-GFT
    for (int k = 0; k < NUM_BTFTJ_GIVENS_TO_RUN; k++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        gft_z8x8_btf_tj(buffer_in[i], buffer_out_btftj[i], givens_btftj[k]);
      t_btftj[k] += clock() - t_temp;
      // compute GFT coefficients error
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < n; j++) {
          adiff = fabs(buffer_out_mat[i][j]) - fabs(buffer_out_btftj[i][j]);
          acc_error_btftj[k] += adiff * adiff;
        }
      }
    }

    // butterfly + PTJ-GFT
    for (int k = 0; k < NUM_BTFPTJ_LAYERS_TO_RUN; k++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        gft_z8x8_btf_ptj(buffer_in[i], buffer_out_btfptj[i],
                         givens_layers_btfptj[k]);
      t_btfptj[k] += clock() - t_temp;
      // compute GFT coefficients error
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < n; j++) {
          adiff = fabs(buffer_out_mat[i][j]) - fabs(buffer_out_btfptj[i][j]);
          acc_error_btfptj[k] += adiff * adiff;
        }
      }
    }

#if CONFIG_DEBUG
    // write output GFT coefficients
    for (int i = 0; i < cur_batch_size; i++) {
      fprintf(fp_out, "Input #%d: \n", b * BATCH_SIZE + i);
      fprintf(fp_out, "Matrix GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_mat[i][j]);
      fprintf(fp_out, "\nButterfly GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_btf[i][j]);
      fprintf(fp_out, "\nTJ-GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_tj[i][j]);
      fprintf(fp_out, "\nPTJ-GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_ptj[i][j]);
      fprintf(fp_out, "\nBTFTJ-GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_btftj[i][j]);
      fprintf(fp_out, "\nBTFPTJ-GFT: ");
      for (int j = 0; j < n; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_btfptj[i][j]);
      fprintf(fp_out, "\n");
    }
#endif
  }

  // write runtime
  fprintf(fp_out, "#input = %d\n", n_inputs);
  fprintf(fp_out, "Matrix GFT:    %.8lf\n", ((double)t_mat) / CLOCKS_PER_SEC);
  fprintf(fp_out, "Butterfly GFT: %.8lf\n", ((double)t_btf) / CLOCKS_PER_SEC);

  for (int k = 0; k < NUM_TJ_GIVENS_TO_RUN; k++)
    fprintf(fp_out,
            "TJ-GFT:        %.8lf, (%d Givens rotations, error = %.8lf)\n",
            ((double)t_tj[k]) / CLOCKS_PER_SEC, givens_tj[k],
            acc_error_tj[k] / ((double)n_inputs));

  for (int k = 0; k < NUM_PTJ_LAYERS_TO_RUN; k++)
    fprintf(fp_out, "PTJ-GFT:        %.8lf, (%d layers, error = %.8lf)\n",
            ((double)t_ptj[k]) / CLOCKS_PER_SEC, givens_layers_ptj[k],
            acc_error_ptj[k] / ((double)n_inputs));

  for (int k = 0; k < NUM_BTFTJ_GIVENS_TO_RUN; k++)
    fprintf(fp_out,
            "BTFTJ-GFT:     %.8lf, (%d Givens rotations, error = %.8lf)\n",
            ((double)t_btftj[k]) / CLOCKS_PER_SEC, givens_btftj[k],
            acc_error_btftj[k] / ((double)n_inputs));

  for (int k = 0; k < NUM_BTFPTJ_LAYERS_TO_RUN; k++)
    fprintf(fp_out, "BTFPTJ-GFT:    %.8lf, (%d layers, error = %.8lf)\n",
            ((double)t_btfptj[k]) / CLOCKS_PER_SEC, givens_layers_btfptj[k],
            acc_error_btfptj[k] / ((double)n_inputs));

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
