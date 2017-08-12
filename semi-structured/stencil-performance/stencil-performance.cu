#include <stdio.h>
#include <math.h>

#include "semi-structured_lib.h"
#include "stencils.h"
#include "make_rooms.h"

#define T 1000

int test1();
int time_sphere(int diam, int big_t);

int main() {
    printf("TIMING SPHERE\n");
    if( !time_sphere(304,500)) {
        printf("SUCCESS\n");
    }

    return 0;
}

int time_sphere(int diam, int big_t) {
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);



    int X,Y,Z;
    struct block *aos;

	struct block *aos_processed;
    
    int blocks_in;
    int *index_of_struct;
    ss_t data;
    char* k_arr;
    k_arr = make_sphere(diam);
    X = Y = Z = diam;
    data = create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);

    int block_ind;
    int arrind[3];
    int coords[3];
    coords[0] = 5;
    coords[1] = 5;
    coords[2] = 5;
    get_coords(coords, X, Y, Z, index_of_struct, &block_ind, arrind);

	struct block *aos_d;
	size_t total_mem = sizeof(struct block) * blocks_in;

	CUCALL(cudaMalloc((void** ) &aos_d, total_mem));
	CUCALL(cudaGetLastError());
	printf("Copying data from host to device\n");
	CUCALL(cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice));
	CUCALL(cudaGetLastError());
    printf("done\n");
    dim3 dims(Bx,By,Bz);

    cudaEventRecord(start);
    int t;
    for (t=0;t<big_t;t++) {
        perform_stencil<<<1,dims>>>(aos_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
        cudaDeviceSynchronize();
        perform_stencil_b<<<1,dims>>>(aos_d, 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
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
    
    printf("Elapsed time: %f ms\n", milliseconds);
    return 0;
}
