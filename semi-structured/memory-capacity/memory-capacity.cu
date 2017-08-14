#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "semi-structured_lib.h"
#include "make_rooms.h"


int main() {

    //Coefficients
    real l2 = 1.0/3.0; //courant number and courant number squared.
    real l = sqrt(l2);
    real h = 0.15; //grid spacing (m)

    //Execute circle example.
    real radius = 10.0; // meters
    int diam = ceil(radius/h)+2;
    int X, Y, Z;

    diam = 8;
   
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
 
    printf("\nCROSS\n");
    printf("Diam \tInternal Blocks\tSize semi-Structured (MiB) \tSize - Structured (MiB) \tRatio\n");
    for (diam = 8; diam <= 500; diam +=40) {
        k_arr = make_cross(diam);
        Z = diam;
        X = Y = 3*diam;
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

