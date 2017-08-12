#include <stdio.h>
#include <math.h>

#include "semi-structured_lib.h"
#include "stencils.h"
#include "make_rooms.h"

#define T 1000

int time_sphere(int diam, int big_t);
int time_cube(int diam, int big_t);
float time_room_ss(int X, int Y, int Z, int big_t, char* k_arr);
float time_room_s(int X, int Y, int Z, int big_t, char* k_arr);

int main() {
    printf("TIMING SPHERE\n");
    if( !time_sphere(32,2)) {
        printf("SUCCESS\n");
    }

    printf("TIMING CUBE\n");
    if( !time_cube(32,1)) {
        printf("SUCCESS\n");
    }
    return 0;
}

int time_sphere(int diam, int big_t) {

    char* k_arr;
    k_arr = make_sphere(diam);
    int X,Y,Z;
    X = Y = Z = diam;

    time_room_ss(X,Y,Z, big_t, k_arr);

    return 0;
}

int time_cube(int diam, int big_t) {

    char* k_arr;
    k_arr = make_cube(diam);
    int X,Y,Z;
    X = Y = Z = diam;

    //time_room_ss(X,Y,Z, big_t, k_arr);
    time_room_ss(X,Y,Z, big_t, k_arr);

    return 0;
}

float time_room_ss(int X, int Y, int Z, int big_t, char* k_arr) {

    struct block *aos;
	struct block *aos_processed;
    int blocks_in;
    int *index_of_struct;
    ss_t data;
    data = create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);

    int block_ind;
    int arrind[3];
    int coords[3];
    coords[0] = X/2;
    coords[1] = Y/2;
    coords[2] = Z/2;
    get_coords(coords, X, Y, Z, index_of_struct, &block_ind, arrind);
    aos[block_ind].u[arrind[2]][arrind[1]][arrind[0]] = 1.0;

	struct block *aos_d;
	size_t total_mem = sizeof(struct block) * blocks_in;

	CUCALL(cudaMalloc((void** ) &aos_d, total_mem));
    cudaDeviceSynchronize();
	CUCALL(cudaGetLastError());
	printf("Copying data from host to device\n");
	CUCALL(cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();
	CUCALL(cudaGetLastError());
    printf("done\n");
    dim3 dims(Bx,By,Bz);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    int t;
    for (t=0;t<(big_t/2);t++) {
        perform_stencil<<<blocks_in,dims>>>(aos_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
        cudaDeviceSynchronize();
        perform_stencil_b<<<blocks_in,dims>>>(aos_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
        cudaDeviceSynchronize();
    }
    cudaEventRecord(stop);

	CUCALL(cudaGetLastError());
    
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

	aos_processed = (struct block*) calloc(blocks_in, sizeof(struct block));
	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[block_ind].u1[arrind[2]][arrind[1]][arrind[0]]);
    printf("%f\n", aos_processed[block_ind].u[arrind[2]][arrind[1]][arrind[0]]);

    free_ss(data);
    cudaDeviceSynchronize();
    CUCALL(cudaFree(aos_d));
    
    int n_voxels = Bx*By*Bz*blocks_in;
    printf("Elapsed time: %f ms\n", milliseconds);
    printf("%f MiB of data processed %d times in %f seconds\n", ((float) sizeof(struct block)*blocks_in)/( (float) 1024*1024), big_t, milliseconds/1000.0);
    printf("%f voxels/s achieved\n",((float) n_voxels*big_t)/(milliseconds/1000.0));
 
    return milliseconds;
}

float time_room_s(int X, int Y, int Z, int big_t, char* k_arr) {

    real *u1;
    real *u;
	real *u_processed;
    real *u1_processed;

    int coords[3];
    coords[0] = X/2;
    coords[1] = Y/2;
    coords[2] = Z/2;

	real *u_d;
	real *u1_d;
    char*k_d;

	size_t total_mem = sizeof(real) * X*Y*Z;
	size_t total_mem_k = sizeof(real) * X*Y*Z;

    u1 = (real*) calloc(total_mem,1);
    u = (real*) calloc(total_mem,1);

    u1_processed = (real*) calloc(total_mem,1);
    u_processed = (real*) calloc(total_mem,1);

	CUCALL(cudaMalloc((void** ) &u_d, total_mem));
	CUCALL(cudaMalloc((void** ) &u1_d, total_mem));
	CUCALL(cudaMalloc((void** ) &k_d, total_mem_k));
	CUCALL(cudaGetLastError());

	printf("Copying data from host to device\n");
	CUCALL(cudaMemcpy(u1_d, u1, total_mem, cudaMemcpyHostToDevice));
	CUCALL(cudaMemcpy(u_d, u, total_mem, cudaMemcpyHostToDevice));
	CUCALL(cudaGetLastError());
    printf("done\n");
    dim3 dims(Bx,By,Bz);
    dim3 blocks_dims(X/Bx,Y/By,Z/Bz);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    int t;
    for (t=0;t<(big_t/2);t++) {
        perform_stencil_structured<<<blocks_dims,dims>>>(u_d, u1_d, k_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, X,Y,Z);
        cudaDeviceSynchronize();
        perform_stencil_structured<<<blocks_dims,dims>>>(u1_d, u_d, k_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, X,Y,Z);
        cudaDeviceSynchronize();
    }
    cudaEventRecord(stop);

	CUCALL(cudaGetLastError());
    
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

	CUCALL(cudaMemcpy(u1_processed, u1_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaMemcpy(u_processed, u_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", u_processed[coords[2]*X*Y + coords[1]*X + coords[0]]); 
    printf("%f\n", u1_processed[coords[2]*X*Y + coords[1]*X + coords[0]]); 

    free(u1);
    free(u);
    free(u1_processed);
    free(u_processed);
    cudaFree(u1_d);
    cudaFree(u_d);
    cudaFree(k_d);
    
    int n_voxels = X*Y*Z;
    printf("Elapsed time: %f ms\n", milliseconds);
    printf("%f MiB of data processed %d times in %f seconds\n", ((float) sizeof(real)*X*Y*Z)/( (float) 1024*1024), big_t, milliseconds/1000.0);
    printf("%f voxels/s achieved\n",((float) n_voxels*big_t)/(milliseconds/1000.0));
 
    return milliseconds;
}
