#include <stdio.h>
#include <stdlib.h>

#include "riemann/riemann_exact.h"

#define N_CELLS 100
#define TIME_END 0.2f

#define BOUNDARY_TYPE 0 // periodic=0, reflective=1

#define GAMMA 5.f / 3.f

struct Cell {
    float mass;
    float momentum;
    float energy;
    float density;
    float velocity;
    float pressure;
    
    float volume;
    float surface_area;

    float boundary_x;
    float boundary_type;

    struct Cell *ngb_left;
};

void check_conserved(struct Cell cell) {
    if (cell.mass < 0) {
        printf("Mass is less than zero: %f\n", cell.mass);
        exit(-1);
    }
    if (cell.energy < 0) {
        printf("Energy is less than zero: %f\n", cell.energy);
        exit(-1);
    }
}

void conserved_to_primitive(struct Cell *cell) {
    cell->density = cell->mass / cell->volume;
    cell->velocity = cell->momentum / cell->mass;
    cell->pressure = (GAMMA - 1) * 
                                (cell->energy / cell->volume - 
                                 0.5f * cell->density * cell->velocity * cell->velocity);
}

void update(struct Cell *cell, float time_step) {
    // Initialise the data arrays for the left and right states
    float state_left[5];
    struct Cell *cell_ngb_left = cell->ngb_left;
    state_left[0] = cell_ngb_left->density;
    state_left[1] = cell_ngb_left->velocity; // velocity_x
    state_left[2] = 0.f; // velocity_y
    state_left[3] = 0.f; // velocity_z
    state_left[4] = cell_ngb_left->pressure;

    float state_right[5];
    state_right[0] = cell->density;
    state_right[1] = cell->velocity; // velocity_x
    state_right[2] = 0.f; // velocity_y
    state_right[3] = 0.f; // velocity_z
    state_right[4] = cell->pressure;

    // Data array for the solved state
    float state_solved[5];

    // Needed for the SWIFT riemann solver; the unit vector for the face
    float n_unit[3];
    n_unit[0] = 1.f;
    n_unit[1] = 0.f;
    n_unit[2] = 0.f;

    riemann_solver_solve(&state_left[0], &state_right[0], &state_solved[0],
                                        &n_unit[0]);

    float density = state_solved[0];
    float velocity = state_solved[1];
    float pressure = state_solved[4];

    // printf("Solved state: %f, %f, %f\n", density, velocity, pressure);

    // Convert the primitive quantities of the solved state 
    // into fluxes of the conserved qunatities
    float flux_mass = density * velocity;
    float flux_momentum = density * velocity * velocity + pressure;
    float flux_energy = (pressure * GAMMA / (GAMMA - 1) + 
                                    0.5f * density * velocity * velocity) * velocity;

    // printf("Fluxes: %f, %f, %f\n", flux_mass, flux_momentum, flux_energy);

    printf("Before: %f, %f\n", cell_ngb_left->mass, flux_mass);
    printf("Interm: %f, %f\n", cell_ngb_left->mass, cell_ngb_left->surface_area);
    cell_ngb_left->mass -= flux_mass * cell_ngb_left->surface_area * time_step;
    cell_ngb_left->momentum -= flux_momentum * cell_ngb_left->surface_area * time_step;
    cell_ngb_left->energy -= flux_energy * cell_ngb_left->surface_area * time_step;
    printf("After: %f, %f\n", cell_ngb_left->mass, flux_mass);

    cell->mass += flux_mass * cell->surface_area * time_step;
    cell->momentum += flux_momentum * cell->surface_area * time_step;
    cell->energy += flux_energy * cell->surface_area * time_step;
}

void write_to_file(struct Cell *cells, char mode[1]) {
    FILE *fp = fopen("output.txt", mode);
    if (fp == NULL){
        exit(1);
    }

    for (int i = 0; i < N_CELLS; i++){
        if (i == 0){
            fprintf(fp, "%f", cells[i].density);
        }else if (i == N_CELLS - 1){
            fprintf(fp, ", %f\n", cells[i].density);
        } else {
            fprintf(fp, ", %f", cells[i].density);
        }
    }
    fclose(fp);
}

int main(){
    // Generate structure in which we store our cells
    struct Cell *cells = malloc(N_CELLS * sizeof(struct Cell));
    
    // Generate initial conditions
    for (int i = 0; i < N_CELLS; i++) {
        cells[i].volume = 1.f / N_CELLS;
        cells[i].surface_area = 1.f;

        if (i < N_CELLS/2) {
            cells[i].mass = 1.f / N_CELLS;
            cells[i].energy = 1.f / (GAMMA - 1) / N_CELLS;
        } else {
            cells[i].mass = 0.125f / N_CELLS;
            cells[i].energy = 0.1f / (GAMMA - 1) / N_CELLS;
        }

        if (i > 0) {
            cells[i].ngb_left = &cells[i - 1];
        }
    }
    cells[0].ngb_left = &cells[N_CELLS - 1];
    write_to_file(cells, "w");

    int n_iter = 0;
    float time_step = 0.001f;
    float time_current = 0.f;
    // Main time integration loop
    while (time_current < TIME_END){
        for (int i = 0; i < N_CELLS; i++) {
            check_conserved(cells[i]);
            conserved_to_primitive(&cells[i]);
        }
        for (int i = 0; i < N_CELLS; i++) {
              update(&cells[i], time_step);
        }

        n_iter++;
        time_current += time_step;
        printf("Numer of iterations: %d\n", n_iter);
        printf("Current time: %f\n", time_current);
        
        write_to_file(cells, "a");
    }

    for (int i = 0; i < N_CELLS; i++) {
        printf("%f\n", cells[i].density);
    }
}

