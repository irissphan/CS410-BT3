#include "pso.h"
void inform_global(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, pso_settings_t *settings);
void inform(int *comm, double **pos_nb, double **pos_b, double *fit_b, pso_settings_t * settings);
void init_comm_ring(int *comm, pso_settings_t * settings);
void inform_ring(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, pso_settings_t * settings);