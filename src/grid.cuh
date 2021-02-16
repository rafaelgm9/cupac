//
// Created by rafa on 2/5/21.
//

#ifndef CUPAC_GRID_CUH
#define CUPAC_GRID_CUH

#include "define.cuh"

Grid *getGrid(double boxsize, int nside, long long npart, const double *positions);

void gridToOrderedArray(Grid *grid, int nside, double *orderedPositions, long long *numParticlesInGrid, long long *offset);

#endif //CUPAC_GRID_CUH
