#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define POINTS_PER_UNIT 25

#define THETA M_PI/4


typedef struct {
    float value; // floating point value (from physics model)
    char inside; //Boolean char for whether point is inside model square
    int coords[2];
} point;

typedef struct linked_point {
    struct linked_point *up;
    struct linked_point *down;
    struct linked_point *left;
    struct linked_point *right;
    float value;
    char inside;
    int coords[2];
} linked_point_t;

void print_array_of_structs(point *array, int x_dim, int y_dim);

void array_position_to_cart(int j, int i, float* x, float *y, float scale)
{
    *x = ( (float) j )/(POINTS_PER_UNIT) - (0.5*scale );
    *y = ( (float) i )/(POINTS_PER_UNIT) - (0.5*scale );
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
    
    float bounding_scale = cos(THETA) + sin(THETA);

    x_dim = ceil( POINTS_PER_UNIT * bounding_scale );
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
            array_position_to_cart(j, i, &x, &y, bounding_scale);
            bounding_array[i*y_dim + j].inside = is_in_square(THETA,x,y);
        }
    }

    print_array_of_structs(bounding_array, x_dim, y_dim);

    float ratio = ( (float) sizeof(point) )/( (float) sizeof(linked_point_t) );
    printf("ratio of point size to linked point size is %f\n", ratio);


    //Find how many points are inside square.
    int total_inside = 0;
    for (i = 0; i < (x_dim * y_dim); i++)
    {
        if ( bounding_array[i].inside )
        {
            total_inside++;
        }
    }

    printf("There are %i points inside square\n", total_inside);

    return 0;
}

void print_array_of_structs(point *array, int x_dim, int y_dim)
{
    int i, j;
    for (i = 0; i <  y_dim; i++)
    {
        for (j = 0; j < x_dim; j++)
        {
            printf("%i ", array[i*y_dim + j].inside);
        }
        printf("\n");
    }
}
