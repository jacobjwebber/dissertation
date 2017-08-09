#include "stencils.h"

__global__ void perform_IO(real *input_d, real *output_d, real* out_d, real ins, int offset, int t) {
    //Takes two pointers to reals in device memory for input/output locations.
    // These should be calculated from coords elsewhere.
    // sum in source

    *input_d += ins;
    // set output
    out_d[t] = *output_d;
    //*/
}

__global__ void perform_stencil(struct block *aos, real l2, real l, real g, int swap) {
    //Launch with blocksize threads in each dimension.

    int x = threadIdx.x;
    int y = threadIdx.y;
    int z = threadIdx.z;

    int bl = blockIdx.x + 1;

    char k = aos[bl].k[z][y][x];
    __shared__ int left, right, fore, aft, up, down;

    if (x == 0 && y == 0 && z == 0) {
        left = aos[bl].left;
        right = aos[bl].right;
        up = aos[bl].up;
        down = aos[bl].down;
        fore = aos[bl].fore;
        aft = aos[bl].aft;
    }

    //point to arrays in global memory. These should use Cuda broadcast when compiled.

//    bl_array *u1_g, *u_g, *u1_r_g, *u1_l_g, *u1_f_g, *u1_a_g, *u1_u_g, *u1_d_g;

/*   if ( !swap ) {
        u1_g = &(aos[bl].u1);
        u_g  = &(aos[bl].u);

        u1_r_g  = &(aos[aos[bl].right].u1);
        u1_l_g  = &(aos[aos[bl].left ].u1);

        u1_f_g  = &(aos[aos[bl].fore].u1);
        u1_a_g  = &(aos[aos[bl].aft ].u1);

        u1_u_g  = &(aos[aos[bl].up  ].u1);
        u1_d_g  = &(aos[aos[bl].down].u1);
    } else {
        u1_g = &(aos[bl].u);
        u_g  = &(aos[bl].u1);

        u1_r_g  = &(aos[aos[bl].right].u);
        u1_l_g  = &(aos[aos[bl].left ].u);

        u1_f_g  = &(aos[aos[bl].fore].u);
        u1_a_g  = &(aos[aos[bl].aft ].u);

        u1_u_g  = &(aos[aos[bl].up  ].u);
        u1_d_g  = &(aos[aos[bl].down].u);
    }
*/
    __shared__ real u1_s[Bz+2][By+2][Bx+2]; //u1 in shared (L1 cache) memory.

    int i,j,f;
    if (x == 0 && y == 0 && z == 0) {
        for (i = 0; i < Bx; i++) {
            for (j = 0; j < By; j++) {
                for (f = 0; f < Bz; f++) {
                    u1_s[f][j][i] = 0;
                }
            }
        }
    }

    u1_s[z+1][y+1][x+1] = aos[bl].u1[z][y][x];

    if ( x == 0 ) {
        u1_s[z+1][y+1][x+1] = aos[left].u1[Bz-1][y][x];
    } else if( x == Bz-1 ) {
        u1_s[z+1][y+1][x+1] = aos[right].u1[0][y][x];
    }
    if ( y == 0 ) {
        u1_s[z+1][y+1][x+1] = aos[aft].u1[z][By-1][x];
    } if( y == By-1 ) {
        u1_s[z+1][y+1][x+1] = aos[fore].u1[z][0][x];
    }

    if ( z == 0 ) {
        u1_s[z+1][y+1][x+1] = aos[down].u1[z][y][Bx-1];
    } if( z == Bx-1 ) {
        u1_s[z+1][y+1][x+1] = aos[up].u1[z][y][0];
    }

    __syncthreads();

    x++;y++;z++;

    //uzz(IN)=((2-l2*Ki(IN)).*uz(IN) + l2*(uz(iIN+1) + uz(iIN-1) + uz(iIN+Ni) + uz(iIN-Ni) + uz(iIN+Ni*Nj) + uz(iIN-Ni*Nj))+(0.5*l*g*(6-Ki(IN))-1).*uzz(IN))./(1+0.5*l*g*(6-Ki(IN)));

    aos[bl].u[z-1][y-1][x-1] = ((2.0 - l2 * (real) k) * u1_s[z][y][x] +
                                   l2 * ( u1_s[z][y][x+1] +
                                          u1_s[z][y][x-1] +
                                          u1_s[z][y+1][x] +
                                          u1_s[z][y-1][x] +
                                          u1_s[z+1][y][x] +
                                          u1_s[z-1][y][x] ) +
                                          (0.5 * l * g * (6.0 - (real) k) - 1.0) * aos[bl].u[z-1][y-1][x-1])/(1 + 0.5 * l * g * (6 - k));
    if (k == 0) {
        aos[bl].u[z-1][y-1][x-1] = 0.0;
    }

}

__global__ void perform_stencil_structured(real* u, real* u1, char* k_d, real l2, real l, real g, int X, int Y, int Z) {
    // get x,y,z from thread and block Id’s
    int x = blockIdx.x * Bz + threadIdx.x;
    int y = blockIdx.y * By + threadIdx.y;
    int z = blockIdx.z * Bx + threadIdx.z;

    // Test that not at boundary
    if( (x>0) && (x<(X-1))
            && (y>0) && (y<(Y-1))
            && (z>0) && (z<(Z-1)))
    {
        // get linear position
        int cp = z*X*Y+y*X+x;
        char k = k_d[cp];
        u[cp] = ((2 - l2 * k) * u1[cp] +
                l2*(u1[cp-1]+u1[cp+1]+u1[cp-X]+u1[cp+X]+u1[cp-Y*X]+u1[cp+Y*X]) +
                (0.5 * l * g * (6 - k) - 1) * u[cp])/(1 + 0.5 * l * g * (6 - k));
    }
}
