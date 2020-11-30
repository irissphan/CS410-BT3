#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "pso.h"
#include "topologies.h"

#define RNG_UNIFORM() (rand()/(double)RAND_MAX)
#define RNG_UNIFORM_INT(s) (rand()%s)

typedef void (*inform_fun_t)(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, pso_settings_t *settings);
typedef double (*inertia_fun_t)(int step, pso_settings_t *settings);

pso_settings_t *pso_settings_new(int Topology, int dim, double range_lo, double range_hi, int paticle_size, int gen) {
    pso_settings_t *settings = (pso_settings_t *)malloc(sizeof(pso_settings_t));
    if (settings == NULL) { 
        return NULL; 
    }
    settings->dim = dim;
    settings->range_lo = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_lo == NULL) { 
        free(settings); 
        return NULL; 
    }
    settings->range_hi = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_hi == NULL) {
        free(settings); 
        free(settings->range_lo); 
        return NULL; 
    }
    for (int i=0; i<settings->dim; i++) {
        settings->range_lo[i] = range_lo;
        settings->range_hi[i] = range_hi;
    }
    settings->size              = paticle_size;
    settings->steps             = gen;
    settings->c1                = PSO_COGNITIVE;
    settings->c2                = PSO_SOCIAL;
    settings->nhood_strategy    = Topology;
    return settings;
}

void pso_settings_free(pso_settings_t *settings) {
    free(settings->range_lo);
    free(settings->range_hi);
    free(settings);
}

double **pso_matrix_new(int size, int dim) {
    double **m = (double **)malloc(size * sizeof(double *));
    for (int i=0; i<size; i++) {
        m[i] = (double *)malloc(dim * sizeof(double));
    }
    return m;
}

void pso_matrix_free(double **m, int size) {
    for (int i=0; i<size; i++) {
        free(m[i]);
    }
    free(m);
}

void pso_solve(pso_obj_fun_t obj_fun, pso_result_t *solution, pso_settings_t *settings) {
    double **pos    = pso_matrix_new(settings->size, settings->dim);
    double **vel    = pso_matrix_new(settings->size, settings->dim);
    double **pos_b  = pso_matrix_new(settings->size, settings->dim); 
    double *fit     = (double *)malloc(settings->size * sizeof(double));
    double *fit_b   = (double *)malloc(settings->size * sizeof(double));
    double **pos_nb = pso_matrix_new(settings->size, settings->dim); 
    int *comm       = (int *)malloc(settings->size * settings->size * sizeof(int));

    int i, d, step;
    double a, b;
    double rho1, rho2; 
    double w = PSO_INERTIA;
    inform_fun_t inform_fun = NULL;


    srand(19520166);
    switch (settings->nhood_strategy) {
        case 0:
            inform_fun = inform_global;
            break;
        case 1:
            init_comm_ring(comm, settings);
            inform_fun = inform_ring;
            break;
    }
    solution->error = DBL_MAX;
    for (i=0; i<settings->size; i++) {
        for (d=0; d<settings->dim; d++) {
            a = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) * RNG_UNIFORM();
            b = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) * RNG_UNIFORM();
            pos[i][d] = a;
            pos_b[i][d] = a;
            vel[i][d] = (a-b) / 2.;
        }
        fit[i] = obj_fun(pos[i], settings->dim);
        fit_b[i] = fit[i]; 
        if (fit[i] < solution->error) {
            solution->error = fit[i];
            memmove((void *)solution->gbest, (void *)pos[i], sizeof(double) * settings->dim);
        }
    }

    for (step=0; step<settings->steps; step++) {
        settings->step = step;
        inform_fun(comm, (double **)pos_nb, (double **)pos_b, fit_b, solution->gbest, settings);
        for (i=0; i<settings->size; i++) {
            for (d=0; d<settings->dim; d++) {
                rho1 = settings->c1 * RNG_UNIFORM();
                rho2 = settings->c2 * RNG_UNIFORM();
                vel[i][d] = w * vel[i][d] +	rho1 * (pos_b[i][d] - pos[i][d]) + rho2 * (pos_nb[i][d] - pos[i][d]);
                pos[i][d] += vel[i][d];
                if (pos[i][d] < settings->range_lo[d]) {
                    pos[i][d] = settings->range_lo[d];
                    vel[i][d] = 0;
                }
                else if (pos[i][d] > settings->range_hi[d]) {
                    pos[i][d] = settings->range_hi[d];
                    vel[i][d] = 0;
                }
            }
            fit[i] = obj_fun(pos[i], settings->dim);
            if (fit[i] < fit_b[i]) {
                fit_b[i] = fit[i];
                memmove((void *)pos_b[i], (void *)pos[i], sizeof(double) * settings->dim);
            }
            if (fit[i] < solution->error) {
                solution->error = fit[i];
                memmove((void *)solution->gbest, (void *)pos[i], sizeof(double) * settings->dim);
            }
            printf("%.5f %.5f ", pos[i][0], pos[i][1]);
        }
        printf("\n");
        // printf("Gen %d (w=%.4f) :: min err=%.5e\n", step, w, solution->error);
    }
    printf("min err=%.5f\n", solution->error);
    printf("gbest=(%.5f, %.5f)", solution->gbest[0], solution->gbest[1]);

    pso_matrix_free(pos, settings->size);
    pso_matrix_free(vel, settings->size);
    pso_matrix_free(pos_b, settings->size);
    pso_matrix_free(pos_nb, settings->size);
    free(comm);
    free(fit);
    free(fit_b);
}