#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "topologies.h"

void inform_global(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, pso_settings_t *settings) {
    int i;
    for (i=0; i<settings->size; i++)
        memmove((void *)pos_nb[i], (void *)gbest, sizeof(double) * settings->dim);
}

void inform(int *comm, double **pos_nb, double **pos_b, double *fit_b, pso_settings_t * settings) {
    int i, j;
    int b_n;
    for (j=0; j<settings->size; j++) {
        b_n = j;
        for (i=0; i<settings->size; i++)
            if (comm[i*settings->size + j] && fit_b[i] < fit_b[b_n])
                b_n = i;
        memmove((void *)pos_nb[j], (void *)pos_b[b_n], sizeof(double) * settings->dim);
    }
}

void init_comm_ring(int *comm, pso_settings_t * settings) {
    int i;
    memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);
    for (i=0; i<settings->size; i++) {
        comm[i*settings->size+i] = 1;
        if (i==0) {
            comm[i*settings->size+i+1] = 1;
            comm[(i+1)*settings->size-1] = 1;
        } else if (i == settings->size-1) {
            comm[i*settings->size] = 1;
            comm[i*settings->size+i-1] = 1;
        } else {
            comm[i*settings->size+i+1] = 1;
            comm[i*settings->size+i-1] = 1;
        }
    }
}

void inform_ring(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, pso_settings_t * settings) {
    inform(comm, pos_nb, pos_b, fit_b, settings);
}