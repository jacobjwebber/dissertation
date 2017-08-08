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
