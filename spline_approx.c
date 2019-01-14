#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.14157
#define NUMBER_OF_NODES 20
#define DEGREE 3
#define PLOT_GRANULITY 1000
#define DOMAIN_SIZE 0.9
#define DOMAIN_START 0.1

double approx_nodes_x[NUMBER_OF_NODES];
double approx_nodes_y[NUMBER_OF_NODES];
double knots[NUMBER_OF_NODES + DEGREE + 1];
double weights[NUMBER_OF_NODES];
double intermediate_points_x[NUMBER_OF_NODES - 1];
double intermediate_points_y[NUMBER_OF_NODES - 1];
double approximated_points_y[NUMBER_OF_NODES - 1];
double approximation_errors[NUMBER_OF_NODES - 1];

void gen_approx_nodes(double (*f)(double));
void bspline(double t, double *res);
void approximate_points(double (*f)(double), FILE *data, FILE *file_int);
double f1(double x);
double f2(double x);
void print_results();
void gen_plot_data(FILE *file, FILE *file_int, double (*f)(double));

void gen_approx_nodes(double (*f)(double))
{
    double h = ((double)DOMAIN_SIZE) / (NUMBER_OF_NODES - 1);
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        approx_nodes_x[i] = ((double)i) * h + DOMAIN_START;
        approx_nodes_y[i] = f(approx_nodes_x[i]);
    }

    int k;
    for(k = 0; k < DEGREE + 1; k++) {
        knots[k] = 0;
    }
    for(; k < NUMBER_OF_NODES; k++) {
        knots[k] = k - DEGREE;
    }
    int w = k;
    for(; k < NUMBER_OF_NODES + DEGREE + 1; k++) {
        knots[k] = w - DEGREE;
    }

    double start = h / 2;
    for (int i = 0; i < NUMBER_OF_NODES - 1; i++)
    {
        intermediate_points_x[i] = start + ((double)i) * h + DOMAIN_START;
        intermediate_points_y[i] = f(intermediate_points_x[i]);
    }

    for(int i = 0; i < NUMBER_OF_NODES; i++)
    {
        weights[i] = 1;
    }
}

void bspline(double t, double *res)
{
    double low = knots[DEGREE];
    double high = knots[NUMBER_OF_NODES];

    double t_new = t * (high - low) + low;
    if(t_new < low || t_new > high) {
        fprintf(stderr, "Not in bounds!\n");
        exit(1);
    }

    int s;
    for(s = DEGREE; s < NUMBER_OF_NODES; s++) {
        if(t_new >= knots[s] && t_new <= knots[s+1]) {
            break;
        }
    }

    double v[NUMBER_OF_NODES][3];
    for(int i = 0; i < NUMBER_OF_NODES; i++) {
        v[i][0] = approx_nodes_x[i] * weights[i];
        v[i][1] = approx_nodes_y[i] * weights[i];
        v[i][2] = weights[i];
    }

    double alpha;
    for(int l = 1; l <= DEGREE + 1; l++) {
        for(int i = s; i > s - DEGREE - 1 + l; i--) {
            alpha = (t_new - knots[i]) / (knots[i + DEGREE + 1 - l] - knots[i]);
            for (int j = 0; j < 2 + 1; j++) {
                v[i][j] = (1 - alpha) * v[i - 1][j] + alpha * v[i][j];
            }
        }
    }

    for(int i = 0; i < 2; i++) {
        res[i] = v[s][i] / v[s][2];
    }
}

void approximate_points(double (*f)(double), FILE *data, FILE *file_int)
{
    gen_approx_nodes(f);
    print_results();
    printf("done\n");
    gen_plot_data(data, file_int, f);
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

void print_results()
{
    printf("Apporximation nodes:\n");
    for (int i = 0; i < NUMBER_OF_NODES; i++)
    {
        printf("    x = %10f  |  y = %10f\n", approx_nodes_x[i], approx_nodes_y[i]);
    }
    printf("Knot vector:\n");
    for (int i = 0; i < NUMBER_OF_NODES + DEGREE + 1; i++)
    {
        printf("    x = %10f  |\n", knots[i]);
    }
    printf("----------------------------------------------------\n\n");
}

void gen_plot_data(FILE *file, FILE *file_int, double (*f)(double))
{
    double h, res[2], x, y;
    char buffer[100];
    for (int i = 0; i < PLOT_GRANULITY; i++)
    {
        h = ((double)1) / (PLOT_GRANULITY - 1);
        x = ((double)i) * h;
        bspline(x, res);
        sprintf(buffer, "%f, %f\n", res[0], res[1]);
        fwrite(buffer, sizeof(char), strlen(buffer), file_int);

        h = ((double)DOMAIN_SIZE) / (PLOT_GRANULITY - 1);
        x = ((double)i) * h + DOMAIN_START;
        y = f(x);
        sprintf(buffer, "%f, %f\n", x, y);
        fwrite(buffer, sizeof(char), strlen(buffer), file);
    }
}

int main(void)
{
    FILE *data;
    FILE *approx_data;
    data = fopen("data.csv", "w+");
    approx_data = fopen("approx_data.csv", "w+");

    approximate_points(f2, data, approx_data);
}