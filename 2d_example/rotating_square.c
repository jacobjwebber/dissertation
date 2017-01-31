#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define POINTS_PER_UNIT 36

#define THETA M_PI/4 


typedef struct {
    float value; // floating point value (from physics model)
    char inside; //Boolean char for whether point is inside model square
    int coords[2];
} point;


void array_position_to_cart(int j, int i, float* x, float *y)
{
    *x = ( (float) j )/(POINTS_PER_UNIT) - (0.5* (cos(THETA) + sin(THETA)) );
    *y = ( (float) i )/(POINTS_PER_UNIT) - (0.5* (cos(THETA) + sin(THETA)) );
}

char is_in_square(float theta, float x, float y)
{
    //Tests if point (x,y) is within unit square rotated by theta.
    //Does this by applying rotation matrix to points coords.

    float x_prime, y_prime;

    x_prime = x*cos(theta) + y*sin(theta);
    y_prime = x*sin(theta) - y*cos(theta);

    if (x_prime <= 0.5 && x_prime >= -0.5
            && y_prime <= 0.5 && y_prime >= -0.5)
    {
        return TRUE;
    } else 
    {
        return FALSE;
    }
}

int main()
{
    int i, j, x_dim, y_dim, count;
   
    //For reasons of rotational symmetry it is not necessary to support 
    // values of Theta > pi/2

    if (THETA > M_PI/2 || THETA < 0)
    {
        printf("Theta out of bounds\n");
        return 0;
    }
    
    x_dim = ceil( POINTS_PER_UNIT * (cos(THETA) + sin(THETA)));
    y_dim = x_dim;

    printf("allocating memory\n");

    printf("theta = %f, \n", THETA);
    printf("Unit square is %ix%i points\n", POINTS_PER_UNIT, POINTS_PER_UNIT);
    printf("bounding box is: %i by %i\n", x_dim, y_dim);

    point *bounding_array = (point *)malloc(x_dim * y_dim * sizeof(point *));
    
    printf("memory allocated\n");
    
    float x, y;

    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            array_position_to_cart(j, i, &x, &y);
            printf("x= %f, y = %f\n", x, y);
            bounding_array[i*y_dim + j].inside = is_in_square(THETA,x,y);
        }
    }
    
    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            printf("%i ", bounding_array[i*y_dim + j].inside);
        }
        printf("\n");
    }

    return 0;
}
