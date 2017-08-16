#include "semi-structured_lib.h"

/*Functions for performing stencils on a GPU. l is courant number. l2 is this squared. g is wall permitivity.*/

__global__ void perform_stencil(struct block *aos, real l2, real l, real g, int swap);

__global__ void perform_stencil_b(struct block *aos, real l2, real l, real g, int swap);

__global__ void perform_IO(real *input_d, real *output_d, real *out_d, real ins, int offset, int t);

__global__ void perform_stencil_structured(real* u, real* u1, char* k_d, real l2, real l, real g, int X, int Y, int Z);
