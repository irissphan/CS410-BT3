#ifndef PSO_H_
#define PSO_H_

#define PSO_INERTIA 0.7298
#define PSO_COGNITIVE 1.49618
#define PSO_SOCIAL 1.49618

typedef struct {

    double error;
    double *gbest;

} pso_result_t;

typedef double (*pso_obj_fun_t)(double *, int);

typedef struct {
    int dim;
    double *range_lo;
    double *range_hi;
    int size;
    int steps;
    int step;
    double c1;
    double c2;
    int nhood_strategy;
} pso_settings_t;

pso_settings_t *pso_settings_new(int Topology, int dim, double range_lo, double range_hi, int paticle_size, int gen);
void pso_settings_free(pso_settings_t *settings);

void pso_solve(pso_obj_fun_t obj_fun, pso_result_t *solution, pso_settings_t *settings);

#endif 