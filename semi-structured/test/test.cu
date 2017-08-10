#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "semi-structured_lib.h"
#include "stencils.h"
#include "make_rooms.h"

int main() {
	struct block *aos;
	struct block *aos_processed;

	aos = (struct block*) calloc(2, sizeof(struct block));
	aos_processed = (struct block*) calloc(2, sizeof(struct block));

	int i, j, k;
	for (i = 0; i < Bx; i++) {
		for (j = 0; j < Bx; j++) {
			for (k = 0; k < Bx; k++) {
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


    int t;

    perform_stencil<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
    cudaDeviceSynchronize();
    //perform_stencil<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
    cudaDeviceSynchronize();

	CUCALL(cudaGetLastError());
	cudaDeviceSynchronize();

	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[1].u1[2][2][2]);
    printf("%f\n", aos_processed[1].u[2][2][2]);

    free(aos_processed);
    free(aos);
    cudaFree(aos_d);


	return 0;
}
