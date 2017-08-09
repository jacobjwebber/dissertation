#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "semi-structured_lib.h"
#include "stencils.h"
#include "make_rooms.h"

void structured_version(int X, int Y, int Z, int big_n, char *is_in_sphere,
		real l, real l2, real g, int coords[3]);

int main() {
	struct block *aos;
	struct block *aos_processed;

	aos = (struct block*) calloc(2, sizeof(struct block));
	aos_processed = (struct block*) calloc(2, sizeof(struct block));

	int i, j, k;
	for (i = 0; i < Bx; i++) {
		for (j = 0; j < Bx; j++) {
			for (k = 0; k < Bx; k++) {
				aos[1].u[k][j][i] = 2.0;
				aos[1].u1[k][j][i] = 8.0;
                aos[0].u1[k][j][i] = 4.0;
				aos[1].k[k][j][i] = 6;
			}
		}
	}

    aos[1].u1[2][2][2] = 1.0;

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
    perform_stencil<<<1,dims>>>(&(aos_d[0]), 1.0/3.0, sqrt(1.0/3.0), 0.1/1.9, 0);
    cudaDeviceSynchronize();

	CUCALL(cudaGetLastError());
	cudaDeviceSynchronize();

	CUCALL(cudaMemcpy(aos_processed, aos_d, total_mem, cudaMemcpyDeviceToHost));
	CUCALL(cudaGetLastError());

    printf("%f\n", aos_processed[1].u1[2][2][2]);
    printf("%f\n", aos_processed[1].u[2][2][2]);


	return 0;
	/*
	 //Coefficients
	 real l2 = 1.0/3.0; //courant number and courant number squared.
	 real l = sqrt(l2);
	 real r = 0.9; //Wall reflection coefficient.
	 real g = (1-r)/(1+r);
	 real h = 0.08; //grid spacing (m)
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
	 real *usource = hanning_window(Tn);

	 int big_n = ceil(duration/Ts);
	 printf("there will be %i time steps\n", big_n);

	 //Set Cuda coefficients,
	 dim3 dimsBlocks(blocks_in,1,1);
	 dim3 dimsThreads(Bx,By,Bz);

	 dim3 dimsIO(1,1,1);

	 //Allocate device mem.
	 struct block *aos_d;
	 size_t total_mem =  sizeof(struct block)*blocks_in;
	 float mem_in_KiB = ((float) total_mem) / 1024.0;
	 printf("Allocating %f MiB CUDA memory.\n", mem_in_KiB/1024.0);
	 CUCALL( cudaMalloc((void**) &aos_d, total_mem) );
	 printf("Copying data from host to device\n");
	 CUCALL( cudaMemcpy(aos_d, aos, total_mem, cudaMemcpyHostToDevice) );
	 CUCALL( cudaGetLastError());

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
	 printf("Error\n");
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
	 */
	//return 0;
}

void structured_version(int X, int Y, int Z, int big_n, char *is_in_sphere,
		real l, real l2, real g, int coords[3]) {
	//===================================================
	// Set up grid and blocks
	printf("running basic version of experiment\n");
	int Gl = X / Bz;
	int Gm = Y / By;
	int Gp = Z / Bx;

	dim3 dimBzockInt(Bz, By, Bx);
	dim3 dimGridInt(Gl, Gm, Gp);
	dim3 dimsIO(1, 1, 1);

	size_t mem_size = X * Y * Z * sizeof(real);
	real *u_d, *u1_d, *dummy_ptr;
	char *k_d;

	// Initialise memory on device
	printf("Allocating device memory.\n");
	cudaMalloc(&u_d, mem_size);
	cudaMemset(u_d, 0, mem_size);
	cudaMalloc(&u1_d, mem_size);
	cudaMemset(u1_d, 0, mem_size);
	cudaMalloc(&k_d, X * Y * Z * sizeof(char));
	printf("Copying data to device.\n");

	cudaMemcpy(k_d, &(is_in_sphere[0]), X * Y * Z * sizeof(char),
			cudaMemcpyHostToDevice);

	printf("Copied data to device.\n");

	real *out_d;
	CUCALL(cudaMalloc((void**)&out_d, big_n *sizeof(real)));
	CUCALL(cudaMemset(out_d, 0, big_n *sizeof(real)));

	real *input_d, *output_d;
	input_d = &(u_d[coords[2] * X * Y + coords[1] * X + coords[0]]);
	output_d = &(u_d[coords[2] * X * Y + coords[1] * X + coords[0]]);

	printf("input/output point has %i neighbours.\n",
			is_in_sphere[coords[2] * X * Y + coords[1] * X + coords[0]]);

	struct timeval start, end;
	long secs_used, micros_used;

	int t;
	gettimeofday(&start, NULL);

	perform_IO<<<dimsIO, dimsIO>>>(input_d, output_d, out_d, 1.0, 0, 0);
	cudaDeviceSynchronize();

	for (t = 1; t < big_n; t++) {

		// update pointers
		dummy_ptr = u1_d;
		u1_d = u_d;
		u_d = dummy_ptr;

		input_d = &(u_d[(X * Y * Z / 2)]);
		output_d = &(u_d[(X * Y * Z / 2)]);

		//do stencil
		perform_stencil_structured<<<dimGridInt, dimBzockInt>>>(u_d, u1_d, k_d,
				l2, l, g, X, Y, Z);
		cudaDeviceSynchronize();

		perform_IO<<<dimsIO, dimsIO>>>(input_d, output_d, out_d, 0, 0, t);
		cudaDeviceSynchronize();
		if (t % 1000 == 0)
			printf("#");
	}

	printf("\n");

	real *out = (real *) malloc(big_n * sizeof(real));
	cudaMemcpy(out, out_d, big_n * sizeof(real), cudaMemcpyDeviceToHost);

	gettimeofday(&end, NULL);

	secs_used = (end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
	micros_used = ((secs_used * 1000000) + end.tv_usec) - (start.tv_usec);
	printf("For structured: micros_used: %d\n", micros_used);

	printf("first two elements of out_d: %f %f %f %f\n", out[0], out[1], out[2],
			out[3]);

	cudaFree(u_d);
	cudaFree(u1_d);

	CUCALL(cudaGetLastError());

}
