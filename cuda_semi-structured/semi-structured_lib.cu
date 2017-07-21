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


int create_aos(int X, int Y, int Z, char* is_in_sphere, int* blocks_in, struct block** aos_pp, int** index_of_struct) {

    int num_blocks_z = (Z + Bz - 1)/Bz; //Round up in case end block is half populated.
    int num_blocks_y = (Y + By - 1)/By;
    int num_blocks_x = (X + Bx - 1)/Bx;
    int total_blocks = num_blocks_z*num_blocks_y*num_blocks_x;

    printf("allocating index of struct\n");
    //Create an array storing the location of each block.
    *index_of_struct = (int*) calloc( 2, sizeof(int));//total_blocks, sizeof(int));

    //total number of internal blocks.
    *blocks_in = 1;

    printf("allocating index of struct %i\n", *index_of_struct[1]);

    int i,j,k;
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                if ( is_block_internal(i * Bz, j * By, k * Bx, &(is_in_sphere[0]), X,Y,Z) ) {
                    printf("%i %i %i\n",i,j,k);
                    *index_of_struct[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k] = *blocks_in;
                    *blocks_in++;
                }
            }
        }
    }


    printf("%i blocks are internal out of %i blocks in total\n", blocks_in, total_blocks);

    printf("Allocating host memory for %i blocks\n", blocks_in);

    //Assign a block for all volumes containing points
    // aos is short for array of structs.
    struct block *aos = (struct block *) calloc(*blocks_in, sizeof(struct block));

    if (aos) {
        printf("Memory successfully allocated \n");
    } else {
        printf("Memory allocation error.\n");
        return -1;
    }

    //Copy is in sphere array to k arrrays within blocks.
    printf("Copying data to blocks\n");
    int index;
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                index = *index_of_struct[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k];
                
                if (index != 0) {
                    copy_to_struct(i * Bz, j * By, k * Bx, &(aos[index]), &(is_in_sphere[0]), X, Y, Z);
                }

            }
        }
    }



    // SET LEFT AND RIGHT WITHIN STRUCTS.
 
    printf("Assigning block neighbours\n");
    struct block *bl;
    //idea - let null neighbour = 0 . Leave 0th block empty.
    for (i = 0; i < num_blocks_z; i++) {
        for (j = 0; j < num_blocks_y; j++) {
            for (k = 0; k < num_blocks_x; k++) {
                
                index = *index_of_struct[i*num_blocks_y*num_blocks_x + j*num_blocks_x + k];
                if ( index != 0) {
                    bl = &(aos[index]);

                    if (i>0) {bl->left  = *index_of_struct[(i-1)*num_blocks_y*num_blocks_x + j*num_blocks_x + k];}
                    if (i<(num_blocks_z-1)) {bl->right = *index_of_struct[(i+1)*num_blocks_y*num_blocks_x + j*num_blocks_x + k];}
                    
                    if (j>0) {bl->aft   = *index_of_struct[i*num_blocks_y*num_blocks_x + (j-1)*num_blocks_x + k];}
                    if (j>num_blocks_y-1) { bl->fore  = *index_of_struct[i*num_blocks_y*num_blocks_x + (j+1)*num_blocks_x + k];}

                    if (k>0) { bl->down  = *index_of_struct[i*num_blocks_y*num_blocks_x + j*num_blocks_x + (k-1)];}
                    if (k<num_blocks_x) {bl->up    = *index_of_struct[i*num_blocks_y*num_blocks_x + j*num_blocks_x + (k+1)];}
                }

            }
        }
    }


    return TRUE;

}

int get_coords(int coords[3], int X, int Y, int Z, int** index_of_struct, int* io_block_ind, int arrind[3]) {
    //Translate from structured array index to index of, and within, struct.

    int num_blocks_z = (Z + Bz - 1)/Bz; //Round up in case end block is half populated.
    int num_blocks_y = (Y + By - 1)/By;
    int num_blocks_x = (X + Bx - 1)/Bx;

    //coords array is {x,y,z}

    int blockindx = coords[0]/Bx;
    int blockindy = coords[1]/By;
    int blockindz = coords[2]/Bz;

    *io_block_ind = *index_of_struct[blockindz *num_blocks_x*num_blocks_y + blockindy *num_blocks_x + blockindx];

    arrind[0] = coords[0]%Bx;
    arrind[1] = coords[1]%By;
    arrind[2] = coords[2]%Bz;

    return TRUE;
}
