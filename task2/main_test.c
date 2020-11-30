#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "pso.h"


//==============================================================
//                      TEST FUNCTIONS
//==============================================================

double pso_Rastrigin(double *vec, int dim) {
    double sum = 0;
    int i;
    for (i=0; i<dim; i++)
        sum += 10*dim + (pow(vec[i], 2) - 10 * cos(2 * M_PI * vec[i])) + (pow(vec[i], 2) - 10 * cos(2 * M_PI * vec[i]));
    return sum;
}

double pso_Rosenbrock(double *vec, int dim) {
    double sum = 0;
    int i;
    for (i=0; i<dim-1; i++)
        sum += 100 * pow((vec[i+1] - pow(vec[i], 2)), 2) + pow((1 - vec[i]), 2);
    return sum;
}

//==============================================================

int main(int argc, char **argv) {

    pso_settings_t *settings = NULL;
    pso_obj_fun_t obj_fun = NULL;

    // parse command line argument (function name)
    if (argc == 6) {
        if (strcmp(argv[1], "Rastrigin") == 0) {
            obj_fun = pso_Rastrigin;
            settings = pso_settings_new(strcmp(argv[2], "RING"), atoi(argv[3]), -5.12, 5.12, atoi(argv[4]), atoi(argv[5])/atoi(argv[4])-1);
            printf("Optimizing function: Rastrigin (dim=%d, swarm size=%d)\n", settings->dim, settings->size);

        } else if (strcmp(argv[1], "Rosenbrock") == 0) {
            obj_fun = pso_Rosenbrock;
            settings = pso_settings_new(strcmp(argv[2], "RING"), atoi(argv[3]), -30, 30, atoi(argv[4]), atoi(argv[5])/atoi(argv[4])-1);
            printf("Optimizing function: Rosenbrock (dim=%d, swarm size=%d)\n", settings->dim, settings->size);

        } else {
            printf("Unsupported objective function: %s", argv[1]);
            return 1;
        }
    } else {
        printf("Usage: ./pso [PROBLEM] [TOPOLOGY] [DIM] [SIZE] [NUM_EVAL]\n \
        where:\t[PROBLEM] is optional with values [Rastrigin|Rosenbrock|Eggholder|Ackley]\n \
        \t[TOPOLOGY] is neighborhood toplogy with [STAR|RING]\n \
        \t[DIM] is  number of variables\n \
        \t[SIZE] is number of particles\n \
        \t[NUM_EVAL] is max number of evaluation caller");
        return 1;
    }

    pso_result_t solution;

    int i;
    for (i=0; i<10; i++) {
        solution.gbest = (double *)malloc(settings->dim * sizeof(double));
        pso_solve(obj_fun, &solution, settings, i);
        free(solution.gbest);
    }

    pso_settings_free(settings);
    return 0;
}