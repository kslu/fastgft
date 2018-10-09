#include "gftmatrices.h"
#include "txfm.h"

int main(int argc, char *argv[]) {

  int n = 25;
  int n_inputs;
  double buffer_in[BATCH_SIZE][25];
  double buffer_out_mat[BATCH_SIZE][25];
  double buffer_out_btf[BATCH_SIZE][25];

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_mat = 0, t_btf = 0, t_temp = 0;

  for (int b = 0; b < n_batches; b++) {
    cur_batch_size = b < n_batches - 1 ? BATCH_SIZE : n_inputs - b * BATCH_SIZE;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < n; j++)
        fscanf(fp_in, "%lf", &buffer_in[i][j]);
    }

    // matrix GFT
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      gft_skeleton25_mat(buffer_in[i], buffer_out_mat[i]);
    t_mat += clock() - t_temp;

    // butterfly GFT
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      gft_skeleton25_btf(buffer_in[i], buffer_out_btf[i]);
    t_btf += clock() - t_temp;

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
      fprintf(fp_out, "\n");
    }
#endif
  }

  // write runtime
  double time_mat = ((double)t_mat) / CLOCKS_PER_SEC;
  double time_btf = ((double)t_btf) / CLOCKS_PER_SEC;
  fprintf(fp_out, "#input = %d\n", n_inputs);
  fprintf(fp_out, "Matrix GFT:    %.8lf\n", time_mat);
  fprintf(fp_out, "Butterfly GFT: %.8lf", time_btf);

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
