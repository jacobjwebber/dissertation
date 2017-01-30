#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define POINTS_PER_UNIT 20

#define THETA M_PI/15 


typedef struct {
    float value; // floating point value (from physics model)
    char inside; //Boolean char for whether point is inside model square
    int coords[2];
} point;

void alloc_point_array(int x_dim, int y_dim, point **point_array);

void array_position_to_cart(int j, int i, float* x, float *y)
{
    *x = ( (float) j )/(POINTS_PER_UNIT) - 0.5;
    *y = ( (float) i )/(POINTS_PER_UNIT) - 0.5;
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
   
    //For reasons of rotational symmetry it is not necessary to support values of Theta > pi/2
    if (THETA > M_PI/2 || THETA < 0)
    {
        printf("Theta out of bounds\n");
        return 0;
    }
    
    int bounding_box_x = ceil( POINTS_PER_UNIT * (cos(THETA) + sin(THETA)));
    int bounding_box_y = bounding_box_x;

    printf("allocating memory\n");
    printf("theta = %f, \n", THETA);
    printf("Unit square is %ix%i points\n", POINTS_PER_UNIT, POINTS_PER_UNIT);
    printf("bounding box is: %i by %i\n", bounding_box_x, bounding_box_y);

    //ALLOCATE ARRAY OF POINTS - MAKE INTO FUCTION.
    x_dim = bounding_box_x;
    y_dim = bounding_box_y;

    point **bounding_array = (point **)malloc(x_dim * sizeof(point *));
    
    for (i=0; i < x_dim; i++)
    {
        bounding_array[i] = (point *)malloc(y_dim * sizeof(point));
    }   
    //END FUCTION.

    point *bounding_array2 = (point *)malloc(x_dim * y_dim * sizeof(point *));
    
    printf("memory allocated\n");
    
    float x, y;

    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            array_position_to_cart(j, i, &x, &y);
            bounding_array[i][j].inside = is_in_square(THETA,x,y);
        }
    }
    
    
    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            printf("%i ", bounding_array[i][j].inside);
        }
        printf("\n");
    }

    printf("Printing method 2\n");
    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            array_position_to_cart(j, i, &x, &y);
            bounding_array2 [i*y_dim + j].inside = is_in_square(THETA,x,y);
        }
    }
    
    
    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            printf("%i ", bounding_array2[i*y_dim + j].inside);
        }
        printf("\n");
    }


    return 0;
}

void alloc_point_array(int x_dim, int y_dim, point **point_array)
{
    //point **point_array = (point *)malloc(x_dim * sizeof(point));
    
    int i,j;
    *point_array = malloc(x_dim * sizeof(point *));
    
    for (i=0; i < x_dim; i++)
    {
        point_array[i] = (point *)malloc(y_dim * sizeof(point));
    }
}
    

