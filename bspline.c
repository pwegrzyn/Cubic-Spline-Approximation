#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.14157
#define NUMBER_OF_NODES 10
#define NUMBER_OF_BASE_FUNCTIONS 5
#define DEGREE 3
#define PLOT_GRANULITY 100
#define DOMAIN_SIZE 1
#define DOMAIN_START 0

double coeffs[NUMBER_OF_BASE_FUNCTIONS];
double approx_nodes_x[NUMBER_OF_NODES];
double approx_nodes_y[NUMBER_OF_NODES];
double parameters[NUMBER_OF_NODES];
double knots[NUMBER_OF_NODES + DEGREE + 1];

void gen_approx_nodes(double (*f)(double));
double bspline(int i, int j, double t);
void spline_approximation(double (*f)(double));
void approximate_points(double (*f)(double), FILE *data, FILE *file_int);
double eval_approx_spline(double x);
double f1(double x);
double f2(double x);
void gaussian_elimination(double *a, double *b, double *x, int n);
void matrix_mul_optimized(double *matrix_1, int a, int b, double *matrix_2, int c, int d, double *res);
void print_results();
void gen_plot_data(FILE *file, FILE *file_int, double (*f)(double));
void matrix_vector_mul(double *matrix, int a, int b, double *vect, int n, double *res);
double euclid_dist(double x1, double y1, double x2, double y2);
double calc_knot_val(int j);

double euclid_dist(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double calc_knot_val(int j) {
    double d = NUMBER_OF_NODES / (NUMBER_OF_BASE_FUNCTIONS - 3);
    int i = floor(j * d);
    double alpha = j * d - i;
    return (1 - alpha) * parameters[i - 1] + alpha * parameters[i];
}

void gen_approx_nodes(double (*f)(double))
{
    double h = ((double)DOMAIN_SIZE) / (NUMBER_OF_NODES - 1);
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        approx_nodes_x[i] = ((double)i) * h + DOMAIN_START;
        approx_nodes_y[i] = f(approx_nodes_x[i]);
    }

    // *************************************************
    // IMPORTANT

    // Calculate the parameters vector using normalized accumulated chord parametrization method
    // double total_chord_length = 0.0;
    // for(int i = 1; i < NUMBER_OF_NODES; i++) {
    //     total_chord_length += euclid_dist(approx_nodes_x[i], approx_nodes_y[i], approx_nodes_x[i - 1], approx_nodes_y[i - 1]);
    // }
    // parameters[0] = 0;
    // parameters[NUMBER_OF_NODES - 1] = 1;
    // for(int i = 1; i < NUMBER_OF_NODES - 1; i++) {
    //     parameters[i] = parameters[i - 1] + euclid_dist(approx_nodes_x[i], approx_nodes_y[i], approx_nodes_x[i - 1], approx_nodes_y[i - 1]) / total_chord_length;
    // }
    // Using equidistance points
    // parameters[0] = 0;
    // parameters[NUMBER_OF_NODES - 1] = 1;
    // for(int i = 1; i < NUMBER_OF_NODES - 1; i++) {
    //     parameters[i] = ((double)i) / (NUMBER_OF_NODES - 1);
    // }

    // // // Calculte the knot vector
    // knots[0] = 0;
    // knots[1] = 0;
    // knots[2] = 0;
    // knots[3] = 0;
    // for(int i = 4; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
    //     knots[i] = calc_knot_val(i);
    // }
    // knots[NUMBER_OF_NODES] = 1;
    // knots[NUMBER_OF_NODES + 1] = 1;
    // knots[NUMBER_OF_NODES + 2] = 1;
    // knots[NUMBER_OF_NODES + 3] = 1;

    knots[0] = approx_nodes_x[0];
    knots[1] = approx_nodes_x[0];
    knots[2] = approx_nodes_x[0];
    knots[3] = approx_nodes_x[0];
    for(int i = 4; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        knots[i] = calc_knot_val(i);
    }
    knots[NUMBER_OF_NODES] = 1;
    knots[NUMBER_OF_NODES + 1] = 1;
    knots[NUMBER_OF_NODES + 2] = 1;
    knots[NUMBER_OF_NODES + 3] = 1;

    // double gap = ((double)DOMAIN_SIZE) / (NUMBER_OF_BASE_FUNCTIONS - 1);
    // knots[0] = DOMAIN_START;
    // knots[1] = DOMAIN_START;
    // for (int i = 2; i < NUMBER_OF_BASE_FUNCTIONS + DEGREE + 1 - 2; i++)
    // {
    //     knots[i] = ((double)i-2) * gap + DOMAIN_START;
    // }
    // knots[NUMBER_OF_BASE_FUNCTIONS + DEGREE - 1] = DOMAIN_START + DOMAIN_SIZE;
    // knots[NUMBER_OF_BASE_FUNCTIONS + DEGREE] = DOMAIN_START + DOMAIN_SIZE;
    // Simple
    // knots[0] = 0;
    // knots[NUMBER_OF_BASE_FUNCTIONS + DEGREE] = 1;
    // for(int i = 1; i < NUMBER_OF_BASE_FUNCTIONS + DEGREE; i++) {
    //     knots[i] = ((double)i) / (NUMBER_OF_BASE_FUNCTIONS + DEGREE);
    // }
    // *************************************************
}

double bspline(int i, int j, double t) {
    if(j == 0) {
        if(t >= knots[i] && t < knots[i+1]) {
            return 1;
        } else {
            return 0;
        }
    } else {
        double frac_1 = knots[i + j] - knots[i] != 0 ? (t - knots[i]) / (knots[i + j] - knots[i]) * bspline(i, j - 1, t) : 0;
        double frac_2 = knots[i + j + 1] - knots[i + 1] != 0 ? (knots[i + j + 1] - t) / (knots[i + j + 1] - knots[i + 1]) * bspline(i + 1, j - 1, t) : 0;
        //printf("Bspline for %d %d %f: %f\n", i, j, t, frac_1 + frac_2);
        return frac_1 + frac_2;
    }
}

void spline_approximation(double (*f)(double))
{
    double *a = (double *)malloc(NUMBER_OF_BASE_FUNCTIONS * NUMBER_OF_NODES * sizeof(double));
    double *a_transpose = (double *)malloc(NUMBER_OF_BASE_FUNCTIONS * NUMBER_OF_NODES * sizeof(double));
    double *a_mul_res = (double *)malloc(NUMBER_OF_BASE_FUNCTIONS * NUMBER_OF_BASE_FUNCTIONS * sizeof(double));
    double *b = (double *)malloc(NUMBER_OF_NODES * sizeof(double));
    double *b_mul_res = (double *)malloc(NUMBER_OF_BASE_FUNCTIONS * sizeof(double));
    if (!a || !b || !a_transpose || !a_mul_res || !b_mul_res)
    {
        perror("Error while allocating memory\n");
        exit(1);
    }
    // Generate A matrix
    printf("A matrix:\n");
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        for (int j = 0; j < NUMBER_OF_BASE_FUNCTIONS; j++)
        {
            a[i * NUMBER_OF_BASE_FUNCTIONS + j] = bspline(j, 3, approx_nodes_x[i]);
            printf("%f ", a[i * NUMBER_OF_BASE_FUNCTIONS + j]);
        }
        printf("\n");
    }
    // A transposed
    printf("A matrix transposed:\n");
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        for(int j = 0; j < NUMBER_OF_NODES; j++) {
            a_transpose[i * NUMBER_OF_NODES + j] = a[j * NUMBER_OF_BASE_FUNCTIONS + i];
            printf("%f ", a_transpose[i * NUMBER_OF_NODES + j]);
        }
        printf("\n");
    }
    // Generate b vector
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        b[i] = approx_nodes_y[i];
    }

    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        coeffs[i] = 0;
    }

    matrix_mul_optimized(a_transpose, NUMBER_OF_BASE_FUNCTIONS, NUMBER_OF_NODES, a, NUMBER_OF_NODES, NUMBER_OF_BASE_FUNCTIONS, a_mul_res);

    printf("A mul res matrix:\n");
    for(int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++) {
        for(int j = 0; j < NUMBER_OF_BASE_FUNCTIONS; j++) {
            printf("%f ", a_mul_res[i * NUMBER_OF_BASE_FUNCTIONS + j]);
        }
        printf("\n");
    }

    matrix_vector_mul(a_transpose, NUMBER_OF_BASE_FUNCTIONS, NUMBER_OF_NODES, b, NUMBER_OF_NODES, b_mul_res);

    printf("B mul vector:\n");
    for (int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++)
    {
        printf("%f ", b_mul_res[i]);
    }
    printf("\n");

    gaussian_elimination(a_mul_res, b_mul_res, coeffs, NUMBER_OF_BASE_FUNCTIONS);
    // The coeffs[] array now holds the approximation coefficients
}

void approximate_points(double (*f)(double), FILE *data, FILE *file_int)
{

    gen_approx_nodes(f);
    spline_approximation(f);
    print_results();
    printf("done\n");
    gen_plot_data(data, file_int, f);
}

double eval_approx_spline(double x)
{
    double res = 0;
    for (int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++)
    {
        res += coeffs[i] * bspline(i, 3, x);
    }
    return res;
}

// (0, 4.5)
double f1(double x)
{
    return x / (2 + pow(x, 2));
}

// (0.1, 1)
double f2(double x)
{
    return x * sin(PI / x);
}

void gaussian_elimination(double *a, double *b, double *x, int n)
{

    int column, row, diagonal, max_pivot_row, j;
    double max_pivot, tmp;
    for (diagonal = 0; diagonal < n; diagonal++)
    {
        max_pivot_row = diagonal;
        max_pivot = *(a + (diagonal * n + diagonal)); // i,ith element of the matrix
        for (row = diagonal + 1; row < n; row++)
        {
            tmp = fabs(*(a + (row * n + diagonal)));
            if (tmp > max_pivot)
            {
                max_pivot_row = row;
                max_pivot = tmp;
            }
        }

        if (diagonal != max_pivot_row)
        {
            for (int k = 0; k < n; k++)
            {
                double *tmp_pointer1 = a + (diagonal * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
            tmp = b[diagonal];
            b[diagonal] = b[max_pivot_row];
            b[max_pivot_row] = tmp;
        }

        for (row = diagonal + 1; row < n; row++)
        {
            tmp = *(a + (row * n + diagonal)) / *(a + (diagonal * n + diagonal));
            for (column = diagonal + 1; column < n; column++)
            {
                *(a + (row * n + column)) -= tmp * *(a + (diagonal * n + column));
            }
            *(a + (row * n + diagonal)) = 0;
            b[row] -= tmp * b[diagonal];
        }
    }

    for (row = n - 1; row >= 0; row--)
    {
        tmp = b[row];
        for (j = n - 1; j > row; j--)
        {
            tmp -= x[j] * *(a + (row * n + j));
        }
        x[row] = tmp / *(a + (row * n + row));
    }
}

void matrix_vector_mul(double *matrix, int a, int b, double *vect, int n, double *res) {

    if(b != n) {
        perror("Shapes do not match");
        exit(1);
    }
    for(int i = 0; i < a; i++) {
        res[i] = 0;
        for(int j = 0; j < b; j++) {
            res[i] += matrix[i * b + j] * vect[j];
        }
    }

}

void matrix_mul_optimized(double *matrix_1, int a, int b, double *matrix_2, int c, int d, double *res)
{

    int i, j, k;
    double sum = 0.0;

    for (i = 0; i < a; i++)
    {
        for (j = 0; j < d; j++)
        {
            for (k = 0; k < c; k++)
            {
                sum += matrix_1[i * b + k] * matrix_2[k * d + j];
            }
            res[i * b + j] = sum;
            sum = 0;
        }
    }
}

void print_results()
{
    printf("Apporximation nodes:\n");
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        printf("    x = %10f  |  y = %10f\n", approx_nodes_x[i], approx_nodes_y[i]);
    }
    printf("Knot vector:\n");
    for (int i = 0; i < NUMBER_OF_BASE_FUNCTIONS + DEGREE + 1; i++)
    {
        printf("    x = %10f  |\n", knots[i]);
    }
    printf("Approximation spline coefficients:\n");
    for (int i = 0; i < NUMBER_OF_BASE_FUNCTIONS; i++)
    {
        printf("    a%d = %10f |\n", i, coeffs[i]);
    }
    printf("----------------------------------------------------\n\n");
}

void gen_plot_data(FILE *file, FILE *file_int, double (*f)(double))
{
    double h = ((double)DOMAIN_SIZE) / (PLOT_GRANULITY - 1);
    char buffer[100];
    for (int i = 0; i < PLOT_GRANULITY; i++)
    {
        double x = ((double)i) * h + DOMAIN_START;
        double y = eval_approx_spline(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file_int);

        x = ((double)i) * h + DOMAIN_START;
        y = f(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file);
    }
}

int main(void) {

    FILE *data;
    FILE *approx_data;
    data = fopen("data.csv", "w+");
    approx_data = fopen("approx_data.csv", "w+");

    approximate_points(f1, data, approx_data);

}