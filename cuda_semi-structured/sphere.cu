#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define TRUE 1
#define FALSE 0

#define real double

//Block size
#define Bl 4
#define Bm 8
#define Bp 8

//Bounding box dimensions
#define L 64
#define M 64
#define P 64

//Sphere radius
#define R 64/2

typedef real bl_array[Bl][Bm][Bp];

struct block {
    int up; //z direction
    int down;
    int left; //x direction
    int right;
    int fore; //y direction.
    int aft;
    bl_array u;
    bl_array u1;
    char k[Bl][Bm][Bp];
};

__global__ void perform_stencil(struct block *aos, real l2, real l, real g, int offset);

__global__ void perform_IO(real *input_d, real *output_d, real *out_d, real ins, int offset, int t);

__global__ void perform_stencil_structured(real* u, real* u1, char* k_d, real l2, real l, real g);

int coords_to_index(int x, int y, int z) {
    int area = M*P;
    return x*area + y*P + x;
}

char is_block_internal(int x, int y, int z, char *array, int xmax, int ymax, int zmax) {
    //This function takes coordinates and returns true if any points within a 
    // blocksize starting on that point are inside the room.
    int i, j, k;
    for (i = x; i < x+Bl; i++) {
        for (j = y; j < y+Bm; j++) {
            for (k = z; k < z+Bp; k++) {
                if ( array[i*ymax*zmax + j*zmax + k]) {
                    //printf("boop , %i %i %i \n", x, y, z);
                    return TRUE;
                 }
            }
        }
    }

    return FALSE;
}

char copy_to_struct(int x, int y, int z, struct block *bl, char *array, int xmax, int ymax, int zmax) {
    int i, j, k;
    for (i = x; i < x+Bl; i++) {
        for (j = y; j < y+Bm; j++) {
            for (k = z; k < z+Bp; k++) {
                (*bl).k[i-x][j-y][k-z] = array[i*ymax*zmax + j*zmax + k];
            }
        }
    }
    return TRUE;
}

real *hann(int big_n) {
    real *hanning = (real *) calloc(big_n, sizeof(real));
    if (hanning == NULL) { return NULL; }

    int i;
    for (i=0; i<big_n; i++) {
        hanning[i] = 0.5 * (1.0 - cos( 2 * M_PI * (i/big_n) ) );
    }

    return hanning;
}

int main() {

    //Coefficients
    real l2 = 1.0/3.0; //courant number and courant number squared.
    real l = sqrt(l2);
    real r = 0.9; //Wall reflection coefficient.
    real g = (1-r)/(1+r);
    real h = 0.1; //grid spacing (m)
    real c = 343; //Speed of sound
    real duration = 1.1; //seconds

    //Error checking - Bounding grid must be divisible by block.
    if(L%Bl != 0 || M%Bm != 0 || P%Bp != 0) {
        printf("Not divisible by block\n");
        return 0;
    }

    //Define array of grid points within sphere.
    char is_in_sphere[L][M][P];
    int i,j,k;
    for (i = 0; i < L; i++) {
        for (j = 0; j < M; j++) {
            for (k = 0; k < P; k++) {
                is_in_sphere[i][j][k] = sqrt( (i-R)*(i-R) + (j-R)*(j-R) + (k-R)*(k-R)) < (float) R;
            }
        }
    }

    //Make is in sphere into record of number of neighbours for each point -- needed later.
    for (i = 0; i < L; i++) {
        for (j = 0; j < M; j++) {
            for (k = 0; k < P; k++) {
                if (is_in_sphere[i][j][k]) {
                    is_in_sphere[i][j][k] = 0;
                    if (i+1 <  L && is_in_sphere[i+1][j][k]) { is_in_sphere[i][j][k]++; } 
                    if (i-1 >= 0 && is_in_sphere[i-1][j][k]) { is_in_sphere[i][j][k]++; } 
                    if (j+1 <  M && is_in_sphere[i][j+1][k]) { is_in_sphere[i][j][k]++; } 
                    if (j-1 >= 0 && is_in_sphere[i][j-1][k]) { is_in_sphere[i][j][k]++; } 
                    if (k+1 <  P && is_in_sphere[i][j][k+1]) { is_in_sphere[i][j][k]++; } 
                    if (k-1 >= 0 && is_in_sphere[i][j][k-1]) { is_in_sphere[i][j][k]++; } 
                }
            }
        }
    }

    // BEGIN DATA PREP SECTION

    int num_blocks_l = L/Bl;
    int num_blocks_m = M/Bm;
    int num_blocks_p = P/Bp;
    int total_blocks = num_blocks_l*num_blocks_m*num_blocks_p;

    //Create an array storing the location of each block.
    int *index_of_struct = (int*) calloc(total_blocks, sizeof(int));

    //total number of internal blocks.
    int blocks_in = 1;

    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
                if ( is_block_internal(i * Bl, j * Bm, k * Bp, &(is_in_sphere[0][0][0]), L,M,P) ) {
                    index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + k] = blocks_in;
                    blocks_in++;
                }
            }
        }
    }

    printf("%i blocks are internal out of %i blocks in total\n", blocks_in, total_blocks);

    printf("Allocating host memory for %i blocks\n", blocks_in);

    //Assign a block for all volumes containing points
    // aos is short for array of structs.
    struct block *aos = (struct block *) calloc(blocks_in, sizeof(struct block));

    if (aos) {
        printf("Memory successfully allocated \n");
    } else {
        printf("Memory allocation error.\n");
        return -1;
    }

    //Copy is in sphere array to k arrrays within blocks.

    int index;
 
    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
                index = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + k];
                
                if (index != 0) {
                    copy_to_struct(i * Bl, j * Bm, k * Bp, &(aos[index]), &(is_in_sphere[0][0][0]), L, M, P);
                }

            }
        }
    }



    // SET LEFT AND RIGHT WITHIN STRUCTS.
 
    struct block *bl;
    //idea - let null neighbour = 0 . Leave 0th block empty.
    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
                
                index = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + k];
                if ( index != 0) {
                    bl = &aos[index];

                    bl->left  = index_of_struct[(i-1)*num_blocks_m*num_blocks_p + j*num_blocks_p + k];
                    bl->right = index_of_struct[(i+1)*num_blocks_m*num_blocks_p + j*num_blocks_p + k];
                    
                    bl->aft   = index_of_struct[i*num_blocks_m*num_blocks_p + (j-1)*num_blocks_p + k];
                    bl->fore  = index_of_struct[i*num_blocks_m*num_blocks_p + (j+1)*num_blocks_p + k];

                    bl->down  = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + (k-1)];
                    bl->up    = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + (k+1)];
                }

            }
        }
    }

    //=======================================================
    //=======================================================
    // END OF DATA PREP SECTION - put this into separate file eventually.

    //Use Hanning curve as input.
    real Ts = h*l / c;
    printf("sample rate=%.1f Hz\n", 1/Ts);
    int Tn = floor(10/l);
    real *usource = hann(Tn);

    int big_n = ceil(duration/Ts);
    printf("there will be %i time steps\n", big_n);
 
    //Set Cuda coefficients,
    dim3 dimsBlocks(blocks_in,1,1);
    dim3 dimsThreads(Bl,Bm,Bp);

    dim3 dimsIO(1,1,1);

    //Allocate device mem.
    struct block *aos_d;
    size_t total_mem =  sizeof(struct block)*blocks_in;
    cudaMalloc((void**) &aos_d, total_mem);
    cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice);

    real *out_d;
    cudaMalloc((void**)&out_d, big_n *sizeof(real));
    cudaMemset(out_d, 0, big_n *sizeof(real));

    float mem_in_KiB = ((float) total_mem) / 1024.0;

    //add error checking for malloc and memcopy.
    printf("Allocated and copied %f KiB of data to device successfully.\n", mem_in_KiB);

    int offset;
    //use this for the pointer swap.
    offset = (int) ( ((char*) &(aos[0].u1[0][0][0])) - ((char*) &(aos[0].u[0][0][0])) );

    printf("offset between arrays in struct is %i chars.\n", offset);
 
    
    //Set input and output locations
    // HARDCODED: calculate from COORDS eventually.
    real *input_d = &(aos_d[blocks_in/2].u[Bl/2][Bm/2][Bp/2]);

    real *output_d = &(aos_d[blocks_in/2].u[Bl/2][Bm/2][Bp/2]);
  
    struct timeval start, end;
    long secs_used,micros_used;

    gettimeofday(&start, NULL);
    //Do the stuff you want to time here
    
    perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 1.0, 0, 0); 
    cudaDeviceSynchronize();

    int t;
    for (t = 1; t < big_n; t++) {
        //do stencil
        perform_stencil<<<dimsBlocks,dimsThreads>>>(aos_d, l2, l, g, offset*(t%2));
        cudaDeviceSynchronize();

        perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 0, offset*(t%2), t); 
        cudaDeviceSynchronize();
    }
    
    gettimeofday(&end, NULL);

    secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
    micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
    printf("For semistructured micros_used: %d\n\n",micros_used);

    real *out = (real *) malloc(big_n*sizeof(real));
    cudaMemcpy(out, out_d, big_n*sizeof(real), cudaMemcpyDeviceToHost);

    printf("first two elements of out_d: %f %f\n", out[0], out[1]);

    cudaFree(aos_d);

    //===================================================
    //Basic version.
    // Set up grid and blocks
    int Gl = L/Bl;
    int Gm = M/Bm;
    int Gp = P/Bp;

    dim3 dimBlockInt(Bl, Bm, Bp);
    dim3 dimGridInt(Gl, Gm, Gp);
    size_t mem_size = L*M*P*sizeof(real);
    real *u_d, *u1_d, *dummy_ptr;
    char *k_d;

    // Initialise memory on device
    cudaMalloc(&u_d, mem_size); 
    cudaMemset(u_d, 0, mem_size);
    cudaMalloc(&u1_d, mem_size); 
    cudaMemset(u1_d, 0, mem_size);
    cudaMalloc(&k_d, L*M*P*sizeof(char));
    cudaMemcpy(k_d, &is_in_sphere, L*M*P*sizeof(char), cudaMemcpyHostToDevice);

    input_d = &(u_d[(L*M*P/2)]);
    output_d = &(u_d[(L*M*P/2)]);
    
    printf("input/output point has %i neighbours.\n", is_in_sphere[L/2][M/2][P/2]);

    gettimeofday(&start, NULL);

    perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 1.0, 0, 0); 
    cudaDeviceSynchronize();

    for (t = 1; t < big_n; t++) {
        
        // update pointers
        dummy_ptr = u1_d;
        u1_d = u_d;
        u_d = dummy_ptr;
        
        input_d = &(u_d[(L*M*P/2)]);
        output_d = &(u_d[(L*M*P/2)]);

        //do stencil
        perform_stencil_structured<<<dimGridInt,dimBlockInt>>>(u_d, u1_d, k_d, l2, l, g);
        cudaDeviceSynchronize();

        perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 0, 0, t); 
        cudaDeviceSynchronize();
    }

    cudaMemcpy(out, out_d, big_n*sizeof(real), cudaMemcpyDeviceToHost);

    gettimeofday(&end, NULL);

    secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
    micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
    printf("For structured: micros_used: %d\n",micros_used);


    printf("first two elements of out_d: %f %f\n", out[0], out[1]);

    cudaFree(u_d);
    cudaFree(u_d);
    cudaFree(u1_d);

    return 0;
}

__global__ void perform_IO(real *input_d, real *output_d, real* out_d, real ins, int offset, int t) {
    //Takes two pointers to reals in device memory for input/output locations.
    // These should be calculated from coords elsewhere.
    // sum in source
    
    input_d  = (real*) ( (char*) input_d + offset);
    output_d  = (real*) ( (char*) output_d + offset);
    
    *input_d += ins;
    // set output
    out_d[t] = *output_d;
}

__global__ void perform_stencil(struct block *aos, real l2, real l, real g, int offset) {
    //Launch with blocksize threads in each dimension.
 
    int x = threadIdx.z; //This is backwards on purpose becauseI named dimensions wrong. :(
    int y = threadIdx.y;
    int z = threadIdx.x;

    int bl = blockIdx.x;

    int k = aos[bl].k[x][y][z];
    //point to arrays in global memory. These should use Cuda broadcast when compiled.

    bl_array *u1_g = (bl_array*) ( (char*) &aos[bl].u1[0][0][0] - offset);
    bl_array *u_g  = (bl_array*) ( (char*) &aos[bl].u[0][0][0] + offset);

    bl_array *u1_r_g  = (bl_array*) ( (char*) &(aos[aos[bl].right].u1[0][0][0]) - offset);
    bl_array *u1_l_g  = (bl_array*) ( (char*) &(aos[aos[bl].left ].u1[0][0][0]) - offset);

    bl_array *u1_f_g  = (bl_array*) ( (char*) &(aos[aos[bl].fore].u1[0][0][0]) - offset);
    bl_array *u1_a_g  = (bl_array*) ( (char*) &(aos[aos[bl].aft ].u1[0][0][0]) - offset);

    bl_array *u1_u_g  = (bl_array*) ( (char*) &(aos[aos[bl].up  ].u1[0][0][0]) - offset);
    bl_array *u1_d_g  = (bl_array*) ( (char*) &(aos[aos[bl].down].u1[0][0][0]) - offset);

    __shared__ real u1_s[Bl+2][Bm+2][Bp+2]; //u1 in shared (L1 cache) memory.

    u1_s[x+1][y+1][z+1] = *u1_g[x][y][z];

    if ( x == 0 ) { 
        u1_s[x][y][z] = *u1_l_g[Bl-1][y][z];
    } else if( x == Bl-1 ) {
        u1_s[x][y][z] = *u1_r_g[0][y][z];
    }
    else if ( y == 0 ) { 
        u1_s[x][y][z] = *u1_a_g[x][Bm-1][z];
    } else if( y == Bm-1 ) {
        u1_s[x][y][z] = *u1_f_g[x][0][z];
    }

    else if ( z == 0 ) {
        u1_s[x][y][z] = *u1_d_g[x][y][Bp-1];
    } else if( z == Bp-1 ) {
        u1_s[x][y][z] = *u1_u_g[x][y][0];
    }
    
    __syncthreads();

    x++;y++;z++;

    *u_g[x-1][y-1][z-1] = ((2 - l2 * k) * u1_s[x][y][z] +
                                   l2 * ( u1_s[x][y][z+1] + 
                                          u1_s[x][y][z-1] + 
                                          u1_s[x][y+1][z] + 
                                          u1_s[x][y-1][z] + 
                                          u1_s[x+1][y][z] +
                                          u1_s[x-1][y][z] ) +
                                          (0.5 * l * g * (6 - k) - 1) * *u_g[x-1][y-1][z-1])/(1 + 0.5 * l * g * (6 - k));

    if (k == 0) {
        *u_g[x-1][y-1][z-1] = 0;
    }
}

__global__ void perform_stencil_structured(real* u, real* u1, char* k_d, real l2, real l, real g) {
    // get x,y,z from thread and block Id’s
    int x = blockIdx.x * Bl + threadIdx.x;
    int y = blockIdx.y * Bm + threadIdx.y;
    int z = blockIdx.z * Bp + threadIdx.z;

    // Test that not at boundary
    if( (x>0) && (x<(L-1))
            && (y>0) && (y<(M-1))
            && (z>0) && (z<(P-1)))
    {
        // get linear position
        int cp = z*M*P+(y*M+x);
        char k = k_d[cp];
        u[cp] = ((2 - l2 * k) * u1[cp] +
                l2*(u1[cp-1]+u1[cp+1]+u1[cp-M]+u1[cp+M]+u1[cp-M*P]+u1[cp+M*P]) +
                (0.5 * l * g * (6 - k) - 1) * u[cp])/(1 + 0.5 * l * g * (6 - k));
    }
}
