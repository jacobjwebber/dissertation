#define ReaL double
//#include "CJW Audio.h"
//#include "CJW Cuda.h"
#include <stdio.h>
// Define Grid dims, modulo Thread block
#define Nl 256
#define Nm 296
#define Np 208
#define area (Nl*Np)
// Define Thread block size
#define Bl 32
#define Bm 4
#define Bp 2
// Define Source and Read
#define Sl 100
#define Sm 80
#define Sp 70
#define Rl 100
#define Rm 140
#define Rp 70
#define PI 3.141
// kernel methods
__global__ void UpDateScheme(ReaL *u, const ReaL * __restrict__ u1, ReaL l2);
__global__ void inout(ReaL *u,ReaL *out,ReaL ins,int n);
//
int main(){

    // Simulation parameters
    ReaL SR = 44100.0;
    int NF = 4410;
    ReaL c = 344.0;
    ReaL k = 1/SR;
    ReaL h = sqrt(3.0)*c*k;
    ReaL l2 = 1.0/3.0;
    //
    // Initialise input
    int n;
    size_t pr_size = sizeof(ReaL);
    int dur = 20;
    ReaL *si_h = (ReaL *)calloc(NF,pr_size);
    for(n=0;n<dur;n++){
        si_h[n] = 0.5*(1.0-cos(2.0*PI*n/(ReaL)dur));
    }
    //
    // Set up grid and blocks
    int Gl = Nl/Bl;
    int Gm = Nm/Bm;
    int Gp = Np/Bp;
    dim3 dimBlockInt(Bl, Bm, Bp);
    dim3 dimGridInt(Gl, Gm, Gp);
    dim3 dimBlockIO(1, 1, 1);
    dim3 dimGridIO(1, 1, 1);
    size_t mem_size = area*Np*pr_size;
    ReaL *out_d, *u_d, *u1_d, *dummy_ptr;
    ReaL ins;
    //
    // Initialise memory on device
    cudaMalloc(&u_d, mem_size); 
    cudaMemset(u_d, 0, mem_size);
    cudaMalloc(&u1_d, mem_size); 
    cudaMemset(u1_d, 0, mem_size);
    cudaMalloc(&u_d, NF*pr_size); 
    cudaMemset(u_d, 0, NF*pr_size);
    //
    // initialise memory on host
    ReaL *out_h = (ReaL *)calloc(NF,pr_size);
    if((out_h == NULL)){
        printf("\nout h memory alloc failed...\n");
        exit(EXIT_FAILURE);
    }
    //
    // Compute scheme
    //
    cudaDeviceSynchronize();
    time_t start = time(NULL);
    for(n=0;n<NF;n++)
    {
        UpDateScheme<<<dimGridInt,dimBlockInt>>>(u_d,u1_d,l2);
        // perform read in out
        ins = si_h[n];
        inout<<<dimGridIO,dimBlockIO>>>(u_d,u_d,ins,n);
        // update pointers
        dummy_ptr = u1_d;
        u1_d = u_d;
        u_d = dummy_ptr;
    }
    // print process time
    //checkLastCUDAError("Kernel");
    cudaDeviceSynchronize();
    time_t end = time(NULL);
    printf("\nProcess time : %ld seconds\n", (end-start));
    // copy result back from device
    cudaMemcpy(out_h, u_d, NF*pr_size, cudaMemcpyDeviceToHost);
    // print last samples, and write output file
    //printLastSamples(out_h, NF, 5);
    //
    // Free memory
    free(si_h);
    free(out_h);

    cudaFree(u_d);
    cudaFree(u_d);
    cudaFree(u1_d);

    exit(EXIT_SUCCESS);
}

//
// Standard 3D update scheme
__global__ void UpDateScheme(ReaL *u, const ReaL * __restrict__ u1, ReaL l2)
{
    // get L,M,P from thread and block Id’s
    int L = blockIdx.x * Bl + threadIdx.x;
    int M = blockIdx.y * Bm + threadIdx.y;
    int P = blockIdx.z * Bp + threadIdx.z;
    // Test that not at boundary
    if( (L>0) && (L<(Nl-1))
            && (M>0) && (M<(Nm-1))
            && (P>0) && (P<(Np-1)))
    {
        // get linear position
        int cp = P*area+(M*Nl+L);
        u[cp] = l2*(u1[cp-1]+u1[cp+1]+u1[cp-Nl]+u1[cp+Nl]+u1[cp-area]+u1[cp+area])-u[cp];
    }
}
//
// read output and sum in input
__global__ void inout(ReaL *u,ReaL *out,ReaL ins,int n)
{
    // sum in source
    u[(Sp*area)+(Sm*Nl+Sl)] += ins;
    // noninterp read out
        out[n] = u[(Rp*area)+(Rm*Nl+Rl)];
}
