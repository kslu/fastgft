#include "gftmatrices.h"
#include "txfm.h"

int main(int argc, char *argv[]) {

  /*
  const double test[15] = {.1, .1, .2, .1, .3, .2, 0, .1,
                           .1, .1, .2, .2, 0,  0,  .2};
  double test_out[15] = {0.0};

  gft_skeleton15_mat(test, test_out);
  gft_skeleton15_btf(test, test_out);
  */

  int n = 16;
  int n_inputs;
  double buffer_in[MAX_NUM_DATA][16];
  double buffer_out_mat[MAX_NUM_DATA][16];
  double buffer_out_btf[MAX_NUM_DATA][16];

  // read inputs
  FILE *fp_in;
  fp_in = fopen(argv[1], "r");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  for (int i = 0; i < n_inputs; i++)
    for (int j = 0; j < n; j++)
      fscanf(fp_in, "%lf", &buffer_in[i][j]);

  // matrix GFT
  clock_t t_mat;
  t_mat = clock();
  for (int i = 0; i < n_inputs; i++)
    gft_bd4x4_mat(buffer_in[i], buffer_out_mat[i]);
  t_mat = clock() - t_mat;

  // butterfly GFT
  clock_t t_btf;
  t_btf = clock();
  for (int i = 0; i < n_inputs; i++)
    gft_bd4x4_btf(buffer_in[i], buffer_out_btf[i]);
  t_btf = clock() - t_btf;

  // write results
  double time_mat = ((double)t_mat) / CLOCKS_PER_SEC;
  double time_btf = ((double)t_btf) / CLOCKS_PER_SEC;

  FILE *fp_out_mat;
  fp_out_mat = fopen(argv[2], "w+");
  fprintf(fp_out_mat, "%.8lf\n", time_mat);
  for (int i = 0; i < n_inputs; i++) {
    for (int j = 0; j < n; j++)
      fprintf(fp_out_mat, "%.8lf ", buffer_out_mat[i][j]);
    fprintf(fp_out_mat, "\n");
  }

  FILE *fp_out_btf;
  fp_out_btf = fopen(argv[3], "w+");
  fprintf(fp_out_btf, "%.8lf\n", time_btf);
  for (int i = 0; i < n_inputs; i++) {
    for (int j = 0; j < n; j++)
      fprintf(fp_out_btf, "%.8lf ", buffer_out_btf[i][j]);
    fprintf(fp_out_btf, "\n");
  }

  fclose(fp_in);
  fclose(fp_out_mat);
  fclose(fp_out_btf);
  return 0;
}
