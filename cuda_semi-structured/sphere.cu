#include <stdio.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define real double

//Block size
#define Bl 8
#define Bm 8
#define Bp 4

//Bounding box dimensions
#define L 64
#define M 64
#define P 64

//Sphere radius
#define R 64/2

struct block {
    int up; //z direction
    int down;
    int left; //x direction
    int right;
    int fore; //y direction.
    int aft;
    float u[Bl][Bm][Bp];
    float u1[Bl][Bm][Bp];
    char k[Bl][Bm][Bp];
};

__global__ void perform_stencil_internal(struct block *u, struct block *u1, real l2, real l, real g);

//DUPLICATE
char is_block_in(int x, int y, int z, char is_in_sphere[L][M][P]) {

    int i,j,k;
    for (i = 0; i < Bl; i++) {
        for (j = 0; j < Bm; j++) {
            for (k = 0; k < Bp; k++) {
                if (is_in_sphere[i][j][k]) {return TRUE;}
            }
        }
    }
    return FALSE;
}

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
    real duration = 0.1; //seconds

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

    /*
    //Print slice in middle
    for (i = 0; i < L; i++) {
        for (j = 0; j < M; j++) {
            printf("%i ", is_in_sphere[i][j][R]);
        }
        printf("\n");
    }
    // */

    // BEGIN DATA PREP SECTION

    int num_blocks_l = L/Bl;
    int num_blocks_m = M/Bm;
    int num_blocks_p = P/Bp;
    int total_blocks = num_blocks_l*num_blocks_m*num_blocks_p;

    //Create an array storing the location of each block.
    int *index_of_struct = (int*) calloc(total_blocks, sizeof(int));

    for (i = 0; i < total_blocks; i++) {
        index_of_struct[i] = -1;
    }


    //total number of internal blocks.
    int blocks_in = 0;

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
                
                if (index != -1) {
                    copy_to_struct(i * Bl, j * Bm, k * Bp, &(aos[index]), &(is_in_sphere[0][0][0]), L, M, P);
                }

            }
        }
    }



    // SET LEFT AND RIGHT WITHIN STRUCTS.
 
    struct block *bl;

    for (i = 0; i < num_blocks_l; i++) {
        for (j = 0; j < num_blocks_m; j++) {
            for (k = 0; k < num_blocks_p; k++) {
                
                index = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + k];
                if ( index != -1) {
                    bl = &aos[index];

                    if (i==0) {bl->left = -1;}
                    else {bl->left  = index_of_struct[(i-1)*num_blocks_m*num_blocks_p + j*num_blocks_p + k];}

                    if (i==num_blocks_l-1) {bl->right = -1;} 
                    else {bl->right = index_of_struct[(i+1)*num_blocks_m*num_blocks_p + j*num_blocks_p + k];}
                    
                    if (j==0) {bl->aft = -1;} 
                    else {bl->aft   = index_of_struct[i*num_blocks_m*num_blocks_p + (j-1)*num_blocks_p + k];}

                    if (j == num_blocks_m-1) {bl->fore = -1;} 
                    else {bl->fore  = index_of_struct[i*num_blocks_m*num_blocks_p + (j+1)*num_blocks_p + k];}

                    if (k==0) {bl->down = -1;} 
                    else {bl->down  = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + (k-1)];}
                    
                    if (j==num_blocks_p-1) {bl->up = -1;} 
                    else {bl->up    = index_of_struct[i*num_blocks_m*num_blocks_p + j*num_blocks_p + (k+1)];}
                }

            }
        }
    }

   /*
    //Print slice in middle
    printf("Printing block blocks_in/2\n");
    for (k=0; k< Bp; k++) {
        for (i = 0; i < Bl; i++) {
            for (j = 0; j < Bm; j++) {
                printf("%i ", array[blocks_in/2].k[i][j][k]);
            }
            printf("\n");
        }
    }
    // */

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
    
    int t;
    for (t = 0; t < big_n; t++) {
        //do stencil
    }

    return 0;
}

__global__ void perform_stencil(struct block *aos, real l2, real l, real g) {
    //Launch with blocksize+2 in each dimension.
    int bl = blockIdx.x;

    int bl_r = aos[bl].right;
    int bl_l = aos[bl].left;
    
    int bl_f = aos[bl].fore;
    int bl_a = aos[bl].aft;

    int bl_u = aos[bl].up;
    int bl_d = aos[bl].down;
    
    int x = threadIdx.x;
    int y = threadIdx.y;
    int z = threadIdx.z;

    __shared__ real arr[Bl+2][Bm+2][Bp+2];

    arr[x][y][z] = 0;

    if ( x == 0 && bl_l != -1) { // this means doing it twice eg x==y==0
        arr[x][y][z] = aos[bl_l].u1[Bl-1][y][z];
    } else if( x == Bl+1 && bl_r != -1) {
        arr[x][y][z] = aos[bl_r].u1[0][y][z];
    }
    else if ( y == 0 && bl_a != -1) { // this means doing it twice eg x==y==0
        arr[x][y][z] = aos[bl_a].u1[x][Bm-1][z];
    } else if( y == Bm+1 && bl_f != -1) {
        arr[x][y][z] = aos[bl_f].u1[x][0][z];
    }

    else if ( z == 0 && bl_d != -1) { // this means doing it twice eg x==y==0
        arr[x][y][z] = aos[bl_d].u1[Bp-1][y][z];
    } else if( z == Bp+1 && bl_u != -1) {
        arr[x][y][z] = aos[bl_u].u1[x][y][0];
    }
    
    else {
        arr[x+1][y+1][z+1] = aos[bl].u1[x][y][z];
    }

    __syncthreads();

    if (x<Bl && y<Bm && z<Bp) {

        int k = aos[bl].k[x][y][z];

        aos[bl].u[x][y][z] = ((2 - l2 * k) * aos[bl].u1[x][y][z] +
                        l2 * ( arr[x][y][z+1] + 
                               arr[x][y][z-1] + 
                               arr[x][y+1][z] + 
                               arr[x][y-1][z] + 
                               arr[x+1][y][z] +
                               arr[x-1][y][z] ) +
                               (0.5 * l * g * (6 - k) - 1) * arr[x][y][z])/(1 + 0.5 * l * g * (6 - k));

        if (k == 0) {
            aos[bl].u[x][y][z] = 0;
        }
    }    
}


__global__ void perform_stencil_internal(struct block *aos, real l2, real l, real g) {
    //Do stencil.
    //internal elements
    //Block size should be dimensions of struct blocks -2 in each direction.

    //The following statement effectively caches the block.
    //FIX THIS using SHARED.???
    int bl = blockIdx.x;
    
    int x = threadIdx.x + 1;
    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;
    int k = aos[bl].k[x][y][z];

    aos[bl].u[x][y][z] = ((2 - l2 * k) * aos[bl].u1[x][y][z] +
                    l2 * ( aos[bl].u1[x][y][z+1] + 
                           aos[bl].u1[x][y][z-1] + 
                           aos[bl].u1[x][y+1][z] + 
                           aos[bl].u1[x][y-1][z] + 
                           aos[bl].u1[x+1][y][z] +
                           aos[bl].u1[x-1][y][z] ) +
                           (0.5 * l * g * (6 - k) - 1) * aos[bl].u[x][y][z])/(1 + 0.5 * l * g * (6 - k));

    if (k == 0) {
        aos[bl].u[x][y][z] = 0;
    }

    
}

__global__ void perform_stencil_surfaces_LR(struct block *u, real l2, real l, real g) {
    //Do stencil.
    //Edge elements left and right.
    int bl = blockIdx.x;
    int bl_r = u[bl].right;
    int bl_l = u[bl].left;
    
    real left, right;

    int y = threadIdx.y + 1;
    int z = threadIdx.z + 1;

    int is_right = threadIdx.x;
   
    if (is_right) {
        
        left = u[bl].u1[Bl-2][y][z];
        if( bl_r != -1 ) {
            right = u[bl_r].u1[0][y][z];
        } else {
            right = 0;
        }

    } else {
        
        right = u[bl].u1[1][y][z];
        if ( bl_l != -1 ) {
            left = u[bl_l].u1[Bl-1][y][z];
        } else {
            left = 0;
        }

    }
    
    
    int k = u[bl].k[0][y][z];
    
    
    u[bl].u[0][y][z] = ((2 - l2 * k) * u[bl].u1[0 + is_right*(Bl-1)][y][z] +
                    l2 * ( u[bl].u1[0 + is_right*(Bl-1)][y][z+1] + 
                           u[bl].u1[0 + is_right*(Bl-1)][y][z-1] + 
                           u[bl].u1[0 + is_right*(Bl-1)][y+1][z] + 
                           u[bl].u1[0 + is_right*(Bl-1)][y-1][z] + 
                           right +
                           left ) +
                           (0.5 * l * g * (6 - k) - 1) * u[bl].u[0 + is_right*(Bl-1)][y][z])/(1 + 0.5 * l * g * (6 - k));
    
   

    if (k == 0) {
        u[bl].u[0][y][z] = 0;
    }
    
}
