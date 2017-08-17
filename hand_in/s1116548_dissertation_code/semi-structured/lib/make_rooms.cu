/*Generates K arrays for the room shapes used in the experiments. Function names are self explanatory.*/
#include "make_rooms.h"
#include <stdio.h>

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

char* make_cross(int diam) {
    // Generates a cross made of 5 cubes. Each has side diam.
    int X, Y, Z;
    Z = diam;
    X = Y = 3*diam;

    char *cross = (char*) calloc(X*Y*Z, sizeof(char));

    int i,j,k;
    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                if (k >= diam && k < 2*diam
                        || j >=diam && j < 2*diam) {
                    cross[i*Y*X + j*X + k] = 1;
                }
            }
        }
    }

    for (i = 0; i < Z; i++) {
        for (j = 0; j < Y; j++) {
            for (k = 0; k < X; k++) {
                if (cross[i*Y*X + j*Y + k]) {
                    cross[i*Y*X + j*Y + k] = 0;
                    if (i+1 <  Z && cross[(i+1)*Y*X + j*X + k]) { cross[i*Y*X + j*X + k]++; } 
                    if (i-1 >= 0 && cross[(i-1)*Y*X + j*X + k]) { cross[i*Y*X + j*X + k]++; } 
                    if (j+1 <  Y && cross[i*Y*X + (j+1)*X + k]) { cross[i*Y*X + j*X + k]++; } 
                    if (j-1 >= 0 && cross[i*Y*X + (j-1)*X + k]) { cross[i*Y*X + j*X + k]++; } 
                    if (k+1 <  X && cross[i*Y*X + j*X + (k+1)]) { cross[i*Y*X + j*X + k]++; } 
                    if (k-1 >= 0 && cross[i*Y*X + j*X + (k-1)]) { cross[i*Y*X + j*X + k]++; } 
                }
            }
        }
    }

    return cross;
}
