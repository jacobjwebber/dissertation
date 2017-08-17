#include <stdio.h>
#include "semi-structured_lib.h"

char is_block_internal(int x, int y, int z, char *array, int xmax, int ymax, int zmax) {
    //This function takes coordinates and returns true if any points within a 
    // blocksize starting on that point are inside the room.
    int i, j, k;
    for (i = x; i < x+Bz; i++) {
        for (j = y; j < y+By; j++) {
            for (k = z; k < z+Bx; k++) {
                if (  i<zmax && j<ymax && k<xmax && array[i*ymax*xmax + j*xmax + k]) {
                    return TRUE;
                 }
            }
        }
    }

    return FALSE;
}

char copy_to_struct(int x, int y, int z, struct block *bl, char *array, int xmax, int ymax, int zmax) {
    int i, j, k;
    for (i = x; i < x+Bz; i++) {
        for (j = y; j < y+By; j++) {
            for (k = z; k < z+Bx; k++) {
                if ( i<zmax && j<ymax && k<xmax ) {
                    bl->k[i-x][j-y][k-z] = array[i*ymax*xmax + j*xmax + k];
                } else {
                    bl->k[i-x][j-y][k-z] = 0;
                }
            }
        }
    }
    return TRUE;
}

int free_ss(ss_t data) {
    //frees all data in ss_t.
    free(data.aos);
    free(data.index_of_struct);
    return 0;
}

ss_t create_aos(int X, int Y, int Z, char* is_in_sphere, int* blocks_in, struct block** aos_pp, int** index_of_struct) {
    //Creates an array of structs for stencil operation.

    ss_t data;

    int num_blocks_z = (Z + Bz - 1)/Bz; //Round up in case end block is half populated.
    int num_blocks_y = (Y + By - 1)/By;
    int num_blocks_x = (X + Bx - 1)/Bx;
    int total_blocks = num_blocks_z*num_blocks_y*num_blocks_x;

    //printf("allocating index of struct\n");
    //Create an array storing the location of each block.
    int* index_of_struct_s = (int*) calloc(total_blocks, sizeof(int));

    //total number of internal blocks.
    int blocks_in_s = 1;

    //printf("allocating index of struct %i\n", index_of_struct_s[1]);

    int i,j,k;
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                if ( is_block_internal(i * Bz, j * By, k * Bx, &(is_in_sphere[0]), X,Y,Z) ) {
                    index_of_struct_s[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k] = blocks_in_s;
                    blocks_in_s++;
                }
            }
        }
    }

    //Assign a block for all volumes containing points
    // aos is short for array of structs.
    struct block *aos = (struct block *) calloc(blocks_in_s, sizeof(struct block));

    if (aos) {
        //printf("Memory successfully allocated \n");
    } else {
        printf("Memory allocation error.\n");
        data.error = -1;
    }

    //Copy is in sphere array to k arrrays within blocks.
    int index;
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                index = index_of_struct_s[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k];
                
                if (index != 0) {
                    copy_to_struct(i * Bz, j * By, k * Bx, &(aos[index]), &(is_in_sphere[0]), X, Y, Z);
                }

            }
        }
    }

    // SET NEIGHBOURS WITHIN STRUCTS.
    struct block *bl;
    //idea - let null neighbour = 0 . Leave 0th block empty.
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                
                index = index_of_struct_s[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k];
                if ( index != 0) {
                    bl = &(aos[index]);

                    if (i>0) {bl->down  = index_of_struct_s[(i-1)*num_blocks_y*num_blocks_x + j*num_blocks_x + k];}
                    if (i<(num_blocks_z-1)) {bl->up = index_of_struct_s[(i+1)*num_blocks_y*num_blocks_x + j*num_blocks_x + k];}
                    
                    if (j>0) {bl->aft   = index_of_struct_s[i*num_blocks_y*num_blocks_x + (j-1)*num_blocks_x + k];}
                    if (j>num_blocks_y-1) { bl->fore  = index_of_struct_s[i*num_blocks_y*num_blocks_x + (j+1)*num_blocks_x + k];}

                    if (k>0) { bl->left  = index_of_struct_s[i*num_blocks_y*num_blocks_x + j*num_blocks_x + (k-1)];}
                    if (k<num_blocks_x) {bl->right    = index_of_struct_s[i*num_blocks_y*num_blocks_x + j*num_blocks_x + (k+1)];}
                }

            }
        }
    }

    *index_of_struct = index_of_struct_s;
    *blocks_in = blocks_in_s;
    *aos_pp = aos;

    data.X = X;
    data.Y = Y;
    data.Z = Z;
    data.blocks_in = blocks_in_s;
    data.aos = aos;
    data.index_of_struct = index_of_struct_s;
    data.k = is_in_sphere;

    return data;

}

int get_coords(int coords[3], int X, int Y, int Z, int* index_of_struct, int* io_block_ind, int arrind[3]) {
    //Translate from structured array index to index of, and within, struct.
    //coords array is input, io_block_ind and arrind array are outputs.

    //int num_blocks_z = (Z + Bz - 1)/Bz; //Round up in case end block is half populated.
    int num_blocks_y = (Y + By - 1)/By;
    int num_blocks_x = (X + Bx - 1)/Bx;

    //coords array is {x,y,z}

    int blockindx = coords[0]/Bx;
    int blockindy = coords[1]/By;
    int blockindz = coords[2]/Bz;

    *io_block_ind = index_of_struct[blockindz *num_blocks_x*num_blocks_y + blockindy *num_blocks_x + blockindx];

    arrind[0] = coords[0]%Bx;
    arrind[1] = coords[1]%By;
    arrind[2] = coords[2]%Bz;

    return TRUE;
}

real *hanning_window(int big_n) {
    //Creates a hanning window for use as an input signal.
    real *hanning = (real *) calloc(big_n, sizeof(real));
    if (hanning == NULL) { return NULL; }

    int i;
    for (i=0; i<big_n; i++) {
        hanning[i] = 0.5 * (1.0 - cos( 2 * M_PI * (i/big_n) ) );
    }

    return hanning;
}
