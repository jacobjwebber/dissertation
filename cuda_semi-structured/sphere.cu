#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "semi-structured_lib.h"

__global__ void perform_stencil(struct block *aos, real l2, real l, real g, int swap);

__global__ void perform_IO(real *input_d, real *output_d, real *out_d, real ins, int offset, int t);

__global__ void perform_stencil_structured(real* u, real* u1, char* k_d, real l2, real l, real g, int X, int Y, int Z);

void structured_version(int X, int Y, int Z, int big_n, char *is_in_sphere, real l, real l2, real g, int coords[3]);

real *hann(int big_n) {
    real *hanning = (real *) calloc(big_n, sizeof(real));
    if (hanning == NULL) { return NULL; }

    int i;
    for (i=0; i<big_n; i++) {
        hanning[i] = 0.5 * (1.0 - cos( 2 * M_PI * (i/big_n) ) );
    }

    return hanning;
}

char* make_sphere(int diam) {
    int X, Y, Z;
    int ori = diam/2 + 1; //Add one so there is buffer round edge- why not?
    int rad = diam/2;

    X = Y = Z = diam;
   
    char *is_in_sphere = (char*) calloc(X*Y*Z, sizeof(char));

    if(!is_in_sphere) {
        printf("error allocating is_in sphere\n");
        return 0;
    }

    //array[i*ymax*xmax + j*xmax + k];
    int i,j,k;
    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                is_in_sphere[i*Y*X + j*Y + k] = sqrt( (i-ori)*(i-ori) + (j-ori)*(j-ori) + (k-ori)*(k-ori)) < (float) rad;
            }
        }
    }

    //Make is in sphere into record of number of neighbours for each point -- needed later.
    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                if (is_in_sphere[i*Y*X + j*Y + k]) {
                    is_in_sphere[i*Y*X + j*Y + k] = 0;
                    if (i+1 <  Z && is_in_sphere[(i+1)*Y*X + j*X + k]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                    if (i-1 >= 0 && is_in_sphere[(i-1)*Y*X + j*X + k]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                    if (j+1 <  Y && is_in_sphere[i*Y*X + (j+1)*X + k]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                    if (j-1 >= 0 && is_in_sphere[i*Y*X + (j-1)*X + k]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                    if (k+1 <  X && is_in_sphere[i*Y*X + j*X + (k+1)]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                    if (k-1 >= 0 && is_in_sphere[i*Y*X + j*X + (k-1)]) { is_in_sphere[i*Y*X + j*X + k]++; } 
                }
            }
        }
    }

    return is_in_sphere;
}

int main() {

    //Coefficients
    real l2 = 1.0/3.0; //courant number and courant number squared.
    real l = sqrt(l2);
    real r = 0.9; //Wall reflection coefficient.
    real g = (1-r)/(1+r);
    real h = 0.15; //grid spacing (m)
    real c = 343; //Speed of sound
    real duration = 0.1; //seconds

    //Execute circle example.
    real radius = 10.0; // meters
    int diam = ceil(radius/h)+2;
    int X, Y, Z;
    int ori = diam/2 + 1; //Add one so there is buffer round edge- why not?

    X = Y = Z = diam;
    char* is_in_sphere = make_sphere(diam);

    // BEGIN DATA PREP SECTION
    struct block *aos;
    int blocks_in;
    int *index_of_struct;
    printf("***Calling library function to create aos***\n");
    create_aos(X, Y, Z, is_in_sphere, &blocks_in, &aos, &index_of_struct);
    printf("***Function returned. Array contains %i structs***\n", blocks_in);
    //end data prep

    char* k_arr;
    printf("Diam \tInternal \tSize \tSize-structured \tRatio\n");
    for (diam = 8; diam <= 800; diam +=80) {
        k_arr = make_sphere(diam);
        X = Y = Z = diam;
        // BEGIN DATA PREP SECTION
        struct block *aos;
        int blocks_in;
        int *index_of_struct;
        create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);
        size_t bl_siz = ((blocks_in-1)*sizeof(struct block));
        size_t s_siz = (X*Y*Z*(2*sizeof(real) + sizeof(char))) ;
        float ratio = ((float) bl_siz) / ((float) s_siz);
        printf("%i \t%i \t%f \t%f \t%f \n", diam, blocks_in, (float) bl_siz/(float)(1024*1024),(float) s_siz/ (float)(1024*1024), ratio);
        //end data prep
        //Set input and output locations
    }
    
    int arrindx;
    int arrindy;
    int arrindz;
    int io_block_ind;
    int arrindices[3];
    int input_coords[3];

    input_coords[0] = ori;
    input_coords[1] = ori;
    input_coords[2] = ori;

    get_coords(input_coords, X, Y, Z, index_of_struct, &io_block_ind, arrindices);

    arrindx = arrindices[0];
    arrindy = arrindices[1];
    arrindz = arrindices[2];


    //Use Hanning curve as input.
    real Ts = h*l / c;
    printf("sample rate=%.1f Hz\n", 1/Ts);
    int Tn = floor(10/l);
    real *usource = hann(Tn);

    int big_n = ceil(duration/Ts);
    printf("there will be %i time steps\n", big_n);
 
    //Set Cuda coefficients,
    dim3 dimsBlocks(blocks_in,1,1);
    dim3 dimsThreads(Bx,By,Bz);

    dim3 dimsIO(1,1,1);

    //Allocate device mem.
    struct block *aos_d;
    size_t total_mem =  sizeof(struct block)*(blocks_in);
    float mem_in_KiB = ((float) total_mem) / 1024.0;
    printf("Allocating %f MiB CUDA memory.\n", mem_in_KiB/1024.0);
    CUCALL( cudaMalloc((void**) &aos_d, total_mem) );
    printf("Copying data from host to device\n");
    CUCALL( cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice) );

    real *out_d;
    CUCALL( cudaMalloc((void**)&out_d, big_n *sizeof(real)) );
    CUCALL( cudaMemset(out_d, 0, big_n *sizeof(real)) );


    //add error checking for malloc and memcopy.
    printf("Allocated and copied %f KiB of data to device successfully.\n", mem_in_KiB);

    CUCALL( cudaGetLastError());
    real *input_d = &(aos_d[io_block_ind].u[arrindz][arrindy][arrindx]);
    real *output_d = &(aos_d[io_block_ind].u[arrindz][arrindy][arrindx]);
  
    printf("input/output point has %i neighbours.\n", aos[io_block_ind].k[arrindz][arrindy][arrindx]);

    struct timeval start, end;
    long secs_used,micros_used;

    CUCALL( cudaGetLastError());
    cudaDeviceSynchronize();
    gettimeofday(&start, NULL);
    
    perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 1.0, 0, 0); //ZERO
    cudaDeviceSynchronize();
    int t;
    big_n = 1;
    for (t = 1; t < big_n; t++) {
        //do stencil
        perform_stencil<<<dimsBlocks,dimsThreads>>>(aos_d, l2, l, g, t%2);
        cudaDeviceSynchronize();
        if ( cudaGetLastError() != cudaSuccess) {
            printf("PJG: Error!\n");
            CUCALL( cudaGetLastError());
            break;
        }

        perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 0, 0, t); 
        cudaDeviceSynchronize();
    }
    
    CUCALL( cudaGetLastError());
    gettimeofday(&end, NULL);

    secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
    micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
    cudaDeviceSynchronize();
    printf("For semistructured micros_used: %d\n\n",micros_used);

    real *out = (real *) malloc(big_n*sizeof(real));
    cudaMemcpy(out, out_d, big_n*sizeof(real), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    printf("first two elements of out_d: %f %f %f %f\n", out[0], out[1], out[2], out[3]);

    cudaError_t er = cudaMemcpy(aos, aos_d, total_mem, cudaMemcpyDeviceToHost);
    
    if ( er != cudaSuccess ) {
        printf("balls\n");
    }
    cudaDeviceSynchronize();
    printf("Element in aos is %f\n", aos[io_block_ind].u[arrindz][arrindy][arrindx]);
    printf("Element in aos is %f\n", aos[io_block_ind+1].u[arrindz][arrindy][arrindx]);

    cudaFree(aos_d);
    cudaFree(out_d);
    free(aos);
    printf("freed cuda and host mem\n");

    //int coordsy[3];
    //coordsy[0] = coordsy[1] = coordsy[2] = ori;

    //structured_version( X, Y, Z, big_n, is_in_sphere, l, l2, g, coordsy);
    
    return 0;
}

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

    int bl = blockIdx.x;

    char k = aos[bl].k[z][y][x];
    //point to arrays in global memory. These should use Cuda broadcast when compiled.
    
    bl_array *u1_g, *u_g, *u1_r_g, *u1_l_g, *u1_f_g, *u1_a_g, *u1_u_g, *u1_d_g;

    if ( !swap ) {
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
 
    __shared__ real u1_s[Bz+2][By+2][Bx+2]; //u1 in shared (L1 cache) memory.

    u1_s[z+1][y+1][x+1] = (*u1_g)[z][y][x];

    if ( x == 0 ) { 
        u1_s[z][y][x] = (*u1_l_g)[Bz-1][y][x];
    } else if( x == Bz-1 ) {
        u1_s[z][y][x] = (*u1_r_g)[0][y][x];
    }
    else if ( y == 0 ) { 
        u1_s[z][y][x] = (*u1_a_g)[z][By-1][x];
    } else if( y == By-1 ) {
        u1_s[z][y][x] = (*u1_f_g)[z][0][x];
    }

    else if ( z == 0 ) {
        u1_s[z][y][x] = (*u1_d_g)[z][y][Bx-1];
    } else if( z == Bx-1 ) {
        u1_s[z][y][x] = (*u1_u_g)[z][y][0];
    }
    
    __syncthreads();

    x++;y++;z++;

    //uzz(IN)=((2-l2*Ki(IN)).*uz(IN) + l2*(uz(iIN+1) + uz(iIN-1) + uz(iIN+Ni) + uz(iIN-Ni) + uz(iIN+Ni*Nj) + uz(iIN-Ni*Nj))+(0.5*l*g*(6-Ki(IN))-1).*uzz(IN))./(1+0.5*l*g*(6-Ki(IN)));
    
    (*u_g)[z-1][y-1][x-1] = ((2.0 - l2 * (real) k) * u1_s[z][y][x] +
                                   l2 * ( u1_s[z][y][x+1] + 
                                          u1_s[z][y][x-1] + 
                                          u1_s[z][y+1][x] + 
                                          u1_s[z][y-1][x] + 
                                          u1_s[z+1][y][x] +
                                          u1_s[z-1][y][x] ) +
                                          (0.5 * l * g * (6.0 - (real) k) - 1.0) * (*u_g)[z-1][y-1][x-1])/(1 + 0.5 * l * g * (6 - k));

    if (k == 0) {
        (*u_g)[z-1][y-1][x-1] = 0.0;
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

void structured_version(int X, int Y, int Z, int big_n, char *is_in_sphere, real l, real l2, real g, int coords[3]) {
    //===================================================
    // Set up grid and blocks
    printf("running basic version of experiment\n");
    int Gl = X/Bz;
    int Gm = Y/By;
    int Gp = Z/Bx;

    dim3 dimBzockInt(Bz, By, Bx);
    dim3 dimGridInt(Gl, Gm, Gp);
    dim3 dimsIO(1,1,1);

    size_t mem_size = X*Y*Z*sizeof(real);
    real *u_d, *u1_d, *dummy_ptr;
    char *k_d;

    // Initialise memory on device
    printf("Allocating device memory.\n");
    cudaMalloc(&u_d, mem_size); 
    cudaMemset(u_d, 0, mem_size);
    cudaMalloc(&u1_d, mem_size); 
    cudaMemset(u1_d, 0, mem_size);
    cudaMalloc(&k_d, X*Y*Z*sizeof(char));
    printf("Copying data to device.\n");

    cudaMemcpy(k_d, &(is_in_sphere[0]), X*Y*Z*sizeof(char), cudaMemcpyHostToDevice);

    printf("Copied data to device.\n");

    real *out_d;
    CUCALL( cudaMalloc((void**)&out_d, big_n *sizeof(real)) );
    CUCALL( cudaMemset(out_d, 0, big_n *sizeof(real)) );


    real *input_d, *output_d;
    input_d = &(u_d[coords[2]*X*Y + coords[1]*X + coords[0]]);
    output_d = &(u_d[coords[2]*X*Y + coords[1]*X + coords[0]]);
    
    printf("input/output point has %i neighbours.\n", is_in_sphere[coords[2]*X*Y + coords[1]*X + coords[0]]);
 
    struct timeval start, end;
    long secs_used,micros_used;
   
    int t;
    gettimeofday(&start, NULL);

    perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 1.0, 0, 0); 
    cudaDeviceSynchronize();
    
    for (t = 1; t < big_n; t++) {
        
        // update pointers
        dummy_ptr = u1_d;
        u1_d = u_d;
        u_d = dummy_ptr;
        
        input_d = &(u_d[(X*Y*Z/2)]);
        output_d = &(u_d[(X*Y*Z/2)]);

        //do stencil
        perform_stencil_structured<<<dimGridInt,dimBzockInt>>>(u_d, u1_d, k_d, l2, l, g, X, Y, Z);
        cudaDeviceSynchronize();

        perform_IO<<<dimsIO,dimsIO>>> (input_d, output_d, out_d, 0, 0, t); 
        cudaDeviceSynchronize();
        if (t%1000 == 0)
            printf("#");
    }

    printf("\n");

    real *out = (real *) malloc(big_n*sizeof(real));
    cudaMemcpy(out, out_d, big_n*sizeof(real), cudaMemcpyDeviceToHost);

    gettimeofday(&end, NULL);

    secs_used=(end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
    micros_used= ((secs_used*1000000) + end.tv_usec) - (start.tv_usec);
    printf("For structured: micros_used: %d\n",micros_used);


    printf("first two elements of out_d: %f %f %f %f\n", out[0], out[1], out[2], out[3]);

    cudaFree(u_d);
    cudaFree(u1_d);
    
    CUCALL( cudaGetLastError());

}
