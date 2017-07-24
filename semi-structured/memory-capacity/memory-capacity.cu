#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "semi-structured_lib.h"

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
char* make_cube(int diam) {
    int X, Y, Z;

    X = Y = Z = diam;
   
    char *cube = (char*) calloc(X*Y*Z, sizeof(char));

    if(!cube) {
        printf("error allocating is_in sphere\n");
        return 0;
    }

    //array[i*ymax*xmax + j*xmax + k];
    int i,j,k;
    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                cube[i*Y*X + j*Y + k] = 1;
            }
        }
    }

    //Make is in sphere into record of number of neighbours for each point -- needed later.
    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                if (cube[i*Y*X + j*Y + k]) {
                    cube[i*Y*X + j*Y + k] = 0;
                    if (i+1 <  Z && cube[(i+1)*Y*X + j*X + k]) { cube[i*Y*X + j*X + k]++; } 
                    if (i-1 >= 0 && cube[(i-1)*Y*X + j*X + k]) { cube[i*Y*X + j*X + k]++; } 
                    if (j+1 <  Y && cube[i*Y*X + (j+1)*X + k]) { cube[i*Y*X + j*X + k]++; } 
                    if (j-1 >= 0 && cube[i*Y*X + (j-1)*X + k]) { cube[i*Y*X + j*X + k]++; } 
                    if (k+1 <  X && cube[i*Y*X + j*X + (k+1)]) { cube[i*Y*X + j*X + k]++; } 
                    if (k-1 >= 0 && cube[i*Y*X + j*X + (k-1)]) { cube[i*Y*X + j*X + k]++; } 
                }
            }
        }
    }

    return cube;
}
int main() {

    //Coefficients
    real l2 = 1.0/3.0; //courant number and courant number squared.
    real l = sqrt(l2);
    real h = 0.15; //grid spacing (m)

    //Execute circle example.
    real radius = 10.0; // meters
    int diam = ceil(radius/h)+2;
    int X, Y, Z;

    // BEGIN DATA PREP SECTION
    struct block *aos;
    int blocks_in;
    int *index_of_struct;
    ss_t data;
    char* k_arr;
    printf("Memory Capacity Experiment\n\n");
    printf("SPHERE\n");
    printf("Diam \tInternal Blocks\tSize semi-Structured (MiB) \tSize - Structured (MiB) \tRatio\n");

    for (diam = 8; diam <= 800; diam +=40) {
        k_arr = make_sphere(diam);
        X = Y = Z = diam;
        // BEGIN DATA PREP SECTION
        data = create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);
        size_t bl_siz = ((blocks_in-1)*sizeof(struct block));
        size_t s_siz = (X*Y*Z*(2*sizeof(real) + sizeof(char))) ;
        float ratio = ((float) bl_siz) / ((float) s_siz);
        printf("%i \t%i \t%f \t%f \t%f \n", diam, blocks_in-1, (float) bl_siz/(float)(1024*1024),(float) s_siz/ (float)(1024*1024), ratio);
        free_ss(data);
        //end data prep
        //Set input and output locations
    }
 
    printf("\nCUBE\n");
    printf("Diam \tInternal Blocks \tSize semi-Structured (MiB) \tSize - Structured (MiB) \tRatio\n");

    for (diam = 8; diam <= 800; diam +=40) {
        k_arr = make_cube(diam);
        X = Y = Z = diam;
        // BEGIN DATA PREP SECTION
        data = create_aos(X, Y, Z, k_arr, &blocks_in, &aos, &index_of_struct);
        size_t bl_siz = ((blocks_in-1)*sizeof(struct block));
        size_t s_siz = (X*Y*Z*(2*sizeof(real) + sizeof(char))) ;
        float ratio = ((float) bl_siz) / ((float) s_siz);
        printf("%i \t%i \t%f \t%f \t%f \n", diam, blocks_in-1, (float) bl_siz/(float)(1024*1024),(float) s_siz/ (float)(1024*1024), ratio);
        free_ss(data);
        //end data prep
        //Set input and output locations
    }
   
    return 0;
}

