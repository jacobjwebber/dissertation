#include <stdio.h>
#include <math.h>

#include "semi-structured_lib.h"
#include "stencils.h"
#include "make_rooms.h"

#define T 1000

int test1();
int test2();

int main() {
    printf("RUNNING TEST 1\n");
    if( !test1()) {
        printf("PASSED\n");
    }
    printf("RUNNING TEST 2\n");
    if( !test2()) {
        printf("PASSED\n");
    }

    return 0;
}

int test2() {
    int X,Y,Z;
    int diam = 21;
    struct block *aos;
    int blocks_in;
    int *index_of_struct;
    ss_t data;
    char* k_arr;
    k_arr = make_sphere(diam);
    X = Y = Z = diam;
    data = create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);
    
    int block_ind;
    int arrind[3];
	int i, j, k;
    int coords[3];
	for (i = 0; i < X; i++) {
		for (j = 0; j < Y; j++) {
			for (k = 0; k < Z; k++) {
                coords[0] = i;
                coords[1] = j;
                coords[2] = k;
                get_coords(coords, X, Y, Z, index_of_struct, &block_ind, arrind);
                if (k_arr[ coords[2]*X*Y + coords[1]*X + coords[0] ] 
                    != aos[block_ind].k[arrind[2]][arrind[1]][arrind[0]]) {
                    printf("Test Failed");
                    return -1;
                }
			}
		}
	}

    free_ss(data);
     
    return 0;
}

int test1() {
	struct block *aos;
	struct block *aos_processed;

    printf("Testing stencil. Should print 2.2 1.4 then 1.4 2.2 if dims are all 4\n");
    printf("Bx By and Bz = %d %d and %d\n", Bx, By, Bz);
	aos = (struct block*) calloc(2, sizeof(struct block));
	aos_processed = (struct block*) calloc(2, sizeof(struct block));

	int i, j, k;
	for (i = 0; i < Bx; i++) {
		for (j = 0; j < By; j++) {
			for (k = 0; k < Bz; k++) {
				aos[1].u[k][j][i] = 3.0;
				aos[1].u1[k][j][i] = 2.2;
                aos[0].u1[k][j][i] = 6.0;
				aos[1].k[k][j][i] = 6;
			}
		}
	}


	struct block *aos_d;
	size_t total_mem = sizeof(struct block) * 2;
	CUCALL(cudaMalloc((void** ) &aos_d, total_mem));
	CUCALL(cudaGetLastError());
	printf("Copying data from host to device\n");
	CUCALL(cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice));
	CUCALL(cudaGetLastError());
    printf("done\n");
    dim3 dims(Bx,By,Bz);


    printf("testing forward\n");
    perform_stencil<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
    cudaDeviceSynchronize();

	CUCALL(cudaGetLastError());
	cudaDeviceSynchronize();

	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[1].u1[2][2][2]);
    printf("%f\n", aos_processed[1].u[2][2][2]);


	for (i = 0; i < Bx; i++) {
		for (j = 0; j < Bx; j++) {
			for (k = 0; k < Bx; k++) {
				aos[1].u1[k][j][i] = 3.0;
				aos[1].u[k][j][i] = 2.2;
                aos[0].u[k][j][i] = 6.0;
				aos[1].k[k][j][i] = 6;
			}
		}
	}


	printf("Copying data from host to device\n");
	CUCALL(cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice));
	CUCALL(cudaGetLastError());
    printf("done\n");
 
    printf("testing backward\n");
    perform_stencil_b<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
    cudaDeviceSynchronize();

	CUCALL(cudaGetLastError());

	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[1].u1[2][2][2]);
    printf("%f\n", aos_processed[1].u[2][2][2]);

    printf("Running %d times to ensure stability. Should print sensible numbers (not very high)\n",T);
    int t;
    for (t=0;t<T;t++) {
        perform_stencil<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
        cudaDeviceSynchronize();
        perform_stencil_b<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
        cudaDeviceSynchronize();
    }

	CUCALL(cudaGetLastError());

	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[1].u1[2][2][2]);
    printf("%f\n", aos_processed[1].u[2][2][2]);




    free(aos_processed);
    free(aos);
    cudaFree(aos_d);


	return 0;
}
