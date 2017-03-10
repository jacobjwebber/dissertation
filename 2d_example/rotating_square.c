#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define POINTS_PER_UNIT 25

#define THETA 0 // M_PI/4


struct point {
    float value; // floating point value (from physics model)
    float old_value;
    char inside; //Boolean char for whether point is inside model square
};

struct linked_point {
    struct linked_point *up;
    struct linked_point *down;
    struct linked_point *left;
    struct linked_point *right;
    float value;
    float old_value;
    char inside;
    int coords[2];
};

struct linked_point_small {
    struct linked_point *up;
    struct linked_point *down;
    struct linked_point *left;
    struct linked_point *right;
    float value;
    float old_value;
    char inside;
};
    
void print_array_of_structs(struct point *array, int x_dim, int y_dim);
void print_array_values(struct linked_point *array, int N);

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
            && y_prime <= 0.5 && y_prime >= -0.5) {
        return TRUE;
    } else {
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
        printf("Theta out of bounds. Select a value between 0 and PI/2\n");
        return 0;
    }
    
    float bounding_scale = cos(THETA) + sin(THETA);

    x_dim = ceil( POINTS_PER_UNIT * bounding_scale );
    y_dim = x_dim;

    printf("allocating memory\n");

    printf("theta = %f, \n", THETA);
    printf("Unit square is %ix%i points\n", POINTS_PER_UNIT, POINTS_PER_UNIT);
    printf("bounding box is: %i by %i\n", x_dim, y_dim);

    size_t bound_array_siz = x_dim * y_dim * sizeof(struct point);
    struct point *bounding_array = (struct point *)malloc(bound_array_siz);

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

    float ratio = ( (float) sizeof(struct point) )/( (float) sizeof(struct linked_point) );
    printf("ratio of point size to linked point size is %f\n", ratio);


    //Find how many points are inside square.
    float counter = 0.0;
    int total_inside = 0;
    for (i = 0; i < (x_dim * y_dim); i++)
    {
        if ( bounding_array[i].inside )
        {
            total_inside++;
            bounding_array[i].value = counter;
            counter += 1.0;
        }
    }

    printf("There are %i points inside square\n", total_inside);

    //Allocate memory for new 'compressed' data structure.
    size_t linked_array_siz = total_inside * sizeof(struct linked_point);
    struct linked_point *linked_array = (struct linked_point *)malloc(linked_array_siz);

    j =0;
    for (i = 0; i < (x_dim * y_dim); i++)
    {
        if ( bounding_array[i].inside )
        {
            linked_array[j].value  = bounding_array[i].value;
            linked_array[j].coords[0] = i % (x_dim * y_dim); //bounding_array[i].coords[0];
            linked_array[j].coords[1] = i / (x_dim * y_dim); //bounding_array[i].coords[1];
            j++;
        }
    }

    printf("Copied %i elements to linked array\n", j);

    if ( j != total_inside )
    {
        printf("ERROR, %i elements copied,  %i inside", j, total_inside);
        return -1;
    }

    float siz_percent = 100*( (float)linked_array_siz )/( (float)bound_array_siz );
    printf("linked array uses %f%% of space of bound array.\n", siz_percent);

    //Find Neighbors.
    for (i = 0; i < total_inside; i++)
    {
        linked_array[i].up = NULL;
        linked_array[i].down = NULL;
        linked_array[i].left = NULL;
        linked_array[i].right = NULL;

        for (j = 0; j < total_inside; j++)
        {
            if ( linked_array[j].coords[0] == (linked_array[i].coords[0] + 1)
                    && linked_array[j].coords[1] == (linked_array[i].coords[1] ) )
            {
                linked_array[i].down = &linked_array[j];
            }
            
            if ( linked_array[j].coords[0] == (linked_array[i].coords[0] - 1)
                    && linked_array[j].coords[1] == (linked_array[i].coords[1] ) )
            {
                linked_array[i].up = &linked_array[j];
            }
            
            if ( linked_array[j].coords[0] == (linked_array[i].coords[0] )
                    && linked_array[j].coords[1] == (linked_array[i].coords[1] - 1) )
            {
                linked_array[i].right = &linked_array[j];
            } 

            if ( linked_array[j].coords[0] == (linked_array[i].coords[0] + 1)
                    && linked_array[j].coords[1] == (linked_array[i].coords[1] + 1) )
            {
                linked_array[i].left = &linked_array[j];
            }
        }
    }

    //Use following as test - should print consecutive numbers - write unit test later
    //print_array_values(linked_array, total_inside);


    //Replace linked_array in place with structure that does not include coords.

    size_t linked_array_small_siz = total_inside * sizeof(struct linked_point_small);
    float siz_percent_small = 100*( (float)linked_array_small_siz )/( (float)bound_array_siz );
    printf("linked array: 'small' uses %f%% of space of bound array.\n", siz_percent_small);

    struct linked_point_small *smol_array = (struct linked_point_small *)malloc(linked_array_small_siz);

    //loop that provides small array that does not include coords.
    for (i = 0; i++; i < total_inside) {
        smol_array[i].up = linked_array[i].up;
        smol_array[i].down = linked_array[i].down;
        smol_array[i].left = linked_array[i].left;
        smol_array[i].right = linked_array[i].right;

        smol_array[i].value = linked_array[i].value;
        smol_array[i].old_value = linked_array[i].old_value;
    }

    return 0;
}

void print_array_of_structs(struct point *array, int x_dim, int y_dim)
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

void print_array_values(struct linked_point *array, int N)
{
    int i;
    printf("printing all internal values\n");
    for (i = 0; i < N; i++)
    {
        printf("%f\n", array[i].value);
    }
}
