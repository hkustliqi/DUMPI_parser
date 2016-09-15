/***                 UNCLASSIFIED//FOR OFFICIAL USE ONLY                   ***/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <fftw3.h>
#include <mpi.h>
//#include <mpp/shmem.h>
#include <math.h>

#include "simple_timer.h"

#define PACKAGE_VERSION "1.0.3"

#define  MIN(a,b) ((a) < (b) ? (a) : (b))

#define SHOW 0

void
initialize(float *_M, uint64_t L, uint64_t l, uint64_t K, uint64_t k,
           uint64_t W, uint64_t id)
{
  float (*M)[k][2] = (float (*)[k][2]) _M;

  uint64_t row_per_pe = l / W;
  uint64_t row_minimum = row_per_pe * id;
  uint64_t row_maximum = MIN(row_minimum + row_per_pe, L);

  uint64_t i;
  uint64_t row, col;

  for (row = row_minimum, i = 0; row < row_maximum; ++row, ++i) {
    for (col = 0; col < K; ++col) {
      M[i][col][0] = row + col * L;
      M[i][col][1] = 0.0;
    }
  }
}

void
show_LK(float *_M, uint64_t L, uint64_t l, uint64_t K, uint64_t k, uint64_t W,
        uint64_t id)
{
  float (*M)[k][2] = (float (*)[k][2]) _M;

  uint64_t row_per_pe = l / W;
  uint64_t row_minimum = row_per_pe * id;
  uint64_t row_maximum = row_minimum + row_per_pe;

  uint64_t i, j;
  uint64_t row, col;

  for (j = 0; j < W; ++j) {
    if (id == j) {
      for (row = row_minimum, i = 0; row < row_maximum; ++row, ++i) {
        printf("%4lu:", id);

        for (col = 0; col < k; ++col)
          printf(" %6.1f,%6.1f", M[i][col][0], M[i][col][1]);

        printf("\n");
      }
    }

    //shmem_barrier_all();
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void
show_KL(float *_M, uint64_t L, uint64_t l, uint64_t K, uint64_t k, uint64_t W,
        uint64_t id)
{
  float (*M)[l][2] = (float (*)[l][2]) _M;

  uint64_t row_per_pe = k / W;
  uint64_t row_minimum = row_per_pe * id;
  uint64_t row_maximum = row_minimum + row_per_pe;

  uint64_t i, j;
  uint64_t row, col;

  for (j = 0; j < W; ++j) {
    if (id == j) {
      for (row = row_minimum, i = 0; row < row_maximum; ++row, ++i) {
        printf("%4lu:", id);

        for (col = 0; col < l; ++col)
          printf(" %6.1f,%6.1f", M[i][col][0], M[i][col][1]);

        printf("\n");
      }
    }

    //shmem_barrier_all();
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void
twiddle(float (*M)[2], uint64_t row, uint64_t K, uint64_t N, double sign)
{
  double delta = 2.0 * M_PI * row / N;
  double real = 1.0;            //cos of zero
  double imag = 0.0;            //sin of zero

  double beta = sin(delta);
  delta = sin(delta / 2.0);
  double alpha = 2.0 * delta * delta;

  uint64_t col;

  for (col = 0; col < K; ++col) {
    double x = M[col][0];
    double y = M[col][1];

    M[col][0] = x * real - sign * y * imag;
    M[col][1] = y * real + sign * x * imag;

    double next_real = real - (alpha * real + beta * imag);
    double next_imag = imag - (alpha * imag - beta * real);

    real = next_real;
    imag = next_imag;
  }
}

void
sfft(float *_M, float *_X, uint64_t L, uint64_t l, uint64_t K, uint64_t k,
     uint64_t W, uint64_t id, fftwf_plan fftw_K, fftwf_plan fftw_L)
{
  // Phase 1: M is L x K -- Transform rows, twiddle, and block transpose 
  // over network into X
  float (*M)[k][2] = (float (*)[k][2]) _M;
  float (*X)[k][2] = (float (*)[k][2]) _X;

  uint64_t L_per_pe = l / W;
  uint64_t L_minimum = L_per_pe * id;
  uint64_t L_maximum = MIN(L_minimum + L_per_pe, L);

  uint64_t K_per_pe = k / W;
  uint64_t K_minimum = K_per_pe * id;
  uint64_t K_maximum = MIN(K_minimum + K_per_pe, K);

#if SHOW
  if (id == 0)
    printf("Input L x K\n");
  show_LK(_M, L, l, K, k, W, id);
#endif

  uint64_t i, pe;
  uint64_t src_row, src_row_block;
  uint64_t src_col, src_col_block;
  uint64_t dst_row, dst_row_block;
  uint64_t dst_col, dst_col_block;

  src_row_block = id;

  for (src_row = 0; src_row < L_per_pe; ++src_row) {
    // Transform Row
    fftwf_execute_dft(fftw_K, (fftwf_complex *) & M[src_row][0][0],
                      (fftwf_complex *) & M[src_row][0][0]);

    // Twiddle Row
    twiddle((float (*)[2]) &M[src_row][0][0], L_minimum + src_row, K, L * K,
            -1);

    // Send row to remote PE
    for (src_col = src_col_block = 0; src_col_block < W;
         ++src_col_block, src_col += K_per_pe) {
      dst_row = src_row;
      dst_row_block = src_col_block;
      dst_col_block = src_row_block;
      dst_col = dst_col_block * K_per_pe;

      //shmem_float_put(&X[dst_row][dst_col][0], &M[src_row][src_col][0],
      //                2 * K_per_pe, dst_row_block);
      MPI_Win win;
      MPI_Win_create(&X[dst_row][dst_col][0], 2 * K_per_pe, sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
      MPI_Win_fence(0, win);
      MPI_Put(&M[src_row][src_col][0], 2 * K_per_pe, MPI_FLOAT, dst_row_block, 0,
              2 * K_per_pe, MPI_FLOAT, win);
      MPI_Win_fence(0, win);
    }
  }

  //shmem_barrier_all();
  MPI_Barrier(MPI_COMM_WORLD);

  // Phase 2: Finish transpose X -> M
  float (*bX)[W][K_per_pe][2] = (float (*)[W][K_per_pe][2]) _X;
  float (*bM)[W][L_per_pe][2] = (float (*)[W][L_per_pe][2]) _M;

  uint64_t block;

  for (block = 0; block < W; ++block) {
    for (src_row = 0; src_row < L_per_pe; ++src_row) {
      for (src_col = 0; src_col < K_per_pe; ++src_col) {
        bM[src_col][block][src_row][0] = bX[src_row][block][src_col][0];
        bM[src_col][block][src_row][1] = bX[src_row][block][src_col][1];
      }
    }
  }

#if SHOW
  if (id == 0)
    printf("Post Transpose K x L\n");
  show_KL(_M, L, l, K, k, W, id);
#endif

  // Phase 3: N is now K x L -- Transform rows
  float (*N)[l][2] = (float (*)[l][2]) _M;

  for (src_row = 0; src_row < K_per_pe; ++src_row) {
    fftwf_execute_dft(fftw_L, (fftwf_complex *) & N[src_row][0][0],
                      (fftwf_complex *) & N[src_row][0][0]);
  }

#if SHOW
  if (id == 0)
    printf("Done K x L\n");
  show_KL(_M, L, l, K, k, W, id);
#endif
}

/**
 * Benchmark an N = K x L long 1-D complex fft using floats.
 *
 * sfft routine assumes input values are load in column major order 
 * into an L x K matric of complex floats.
 *
 * Output is also in column major order.
 **/

int
main(int argc, char *argv[])
{
  if (argc != 5) {
    fprintf(stderr, "%s L K iter threads\n", argv[0]);
    return 1;
  }

  uint64_t L = strtoul(argv[1], 0, 0);
  uint64_t K = strtoul(argv[2], 0, 0);
  uint64_t iter = strtoul(argv[3], 0, 0);
  uint64_t thread_count = strtoul(argv[4], 0, 0);
  simple_timer_t T_init;

  //shmem_init();
  simple_timer(&T_init);
  //start_pes(0);
  MPI_Init(&argc, &argv);

  //shmem_barrier_all();
  MPI_Barrier(MPI_COMM_WORLD);
  int64_t init_usec = simple_timer(&T_init);

  fftwf_init_threads();

  //uint64_t W = shmem_n_pes();
  //uint64_t id = shmem_my_pe();
  int id_int;
  int W_int;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_int);
  MPI_Comm_size(MPI_COMM_WORLD, &W_int);
  uint64_t W = (uint64_t) W_int;
  uint64_t id = (uint64_t) id_int; 

  if (id == 0) {
    printf("shmem-fft version %s\n", PACKAGE_VERSION);
    printf("start_pes in %.4f seconds\n", init_usec / 1e6);
  }
  
  uint64_t lw = (L + W - 1) / W;
  uint64_t kw = (K + W - 1) / W;

  uint64_t l = W * lw;
  uint64_t k = W * kw;
  uint64_t i;

  if (id == 0) {
    printf("l = %lu for L = %lu\n", l, L);
    printf("k = %lu for K = %lu\n", k, K);
  }

  uint64_t entries_per_pe = l * k / W;

  //float *A = shmalloc(entries_per_pe * 2 * sizeof(float));
  //float *X = shmalloc(entries_per_pe * 2 * sizeof(float));
  float *A = malloc(entries_per_pe * 2 * sizeof(float));
  float *X = malloc(entries_per_pe * 2 * sizeof(float));

  if (A == NULL || X == NULL) {
    fprintf(stderr, "ERROR > shmalloc failed\n");
    return 1;
  }

#ifdef ENABLE_WARMUP
  {
    if (id == 0)
      printf("Warming up.\n");
    int dest;
    for (i = 0; i < W; i++) {
      dest = (id + i) % W;
      //shmem_float_put(X, &A[dest], 2, dest);
      MPI_Win win;
      MPI_Win_create(X, 2, sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
      MPI_Win_fence(0, win);
      MPI_Put(&A[dest], 2, MPI_FLOAT, dest, 0,
              2, MPI_FLOAT, win);
      MPI_Win_fence(0, win);
    }
  }
#endif

  // Get FFT ready
  if (thread_count > 1)
    fftwf_plan_with_nthreads(thread_count);
  // some computation
  fftwf_plan fftw_K =
      fftwf_plan_dft_1d(K, (fftwf_complex *) A, (fftwf_complex *) A,
                        FFTW_FORWARD, FFTW_MEASURE);
  fftwf_plan fftw_L =
      fftwf_plan_dft_1d(L, (fftwf_complex *) A, (fftwf_complex *) A,
                        FFTW_FORWARD, FFTW_MEASURE);

  simple_timer_t T;

  simple_timer(&T);
  initialize(A, L, l, K, k, W, id);
  //shmem_barrier_all();
  MPI_Barrier(MPI_COMM_WORLD);

  int64_t usec = simple_timer(&T);

  if (id == 0) {
    printf("initialize in %.4f seconds\n", usec / 1e6);
  }
  
  for (i = 0; i < iter; ++i) {
    sfft(A, X, L, l, K, k, W, id, fftw_K, fftw_L);
    usec = simple_timer(&T);

    if (id == 0)
      printf("FFT in %.4f seconds :: %3.3f BRPS  %3.3f GiRPS\n", usec / 1e6,
             L * K / (usec / 1e6) / 1e9,
             L * K / (usec / 1e6) / (UINT32_C(1) << 30));
  }

  //shfree(X);
  //shfree(A);
  free(X);
  free(A);

  fftwf_cleanup_threads();
#ifdef NEED_SHMEM_FINALIZE
  //shmem_finalize();
  MPI_Finalize();
#endif



  return 0;
}

/***                 UNCLASSIFIED//FOR OFFICIAL USE ONLY                 ***/
