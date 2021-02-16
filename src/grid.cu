//
// Created by rafa on 2/5/21.
//

#include "grid.cuh"


Grid *getGrid(double boxsize, int nside, long long npart, const double *positions) {
    Grid *grid;

    grid = (Grid *) malloc(nside * nside * nside * sizeof(Grid));

    for (int i = 0; i < nside * nside * nside; i++) {
        grid[i].np = 0;
    }

    for (long long ii = 0; ii < npart; ii++) {
        int i, j, k;
        int s;

        i = (int) (positions[3 * ii + 0] / boxsize * nside);
        j = (int) (positions[3 * ii + 1] / boxsize * nside);
        k = (int) (positions[3 * ii + 2] / boxsize * nside);

        i %= nside;
        j %= nside;
        k %= nside;

        s = nside * nside * k + nside * j + i;

        grid[s].np++;
    }

    for (int i = 0; i < nside * nside * nside; i++) {
        if (grid[i].np > 0) {
            grid[i].pos = (double *) malloc(3 * grid[i].np * sizeof(double));
        }
        grid[i].np = 0;
    }

    for (long long ii = 0; ii < npart; ii++) {
        int i, j, k;
        int s;
        long long offset;

        i = (int) (positions[3 * ii + 0] / boxsize * nside);
        j = (int) (positions[3 * ii + 1] / boxsize * nside);
        k = (int) (positions[3 * ii + 2] / boxsize * nside);

        i %= nside;
        j %= nside;
        k %= nside;

        s = nside * nside * k + nside * j + i;

        offset = 3 * grid[s].np;
        grid[s].pos[offset + 0] = positions[3 * ii + 0];
        grid[s].pos[offset + 1] = positions[3 * ii + 1];
        grid[s].pos[offset + 2] = positions[3 * ii + 2];
        grid[s].np++;
    }

    return grid;
}

void
gridToOrderedArray(Grid *grid, int nside, double *orderedPositions, long long *numParticlesInGrid, long long *offset) {

    for (int s = 0; s < nside * nside * nside; s++) {
        numParticlesInGrid[s] = grid[s].np;
    }

    offset[0] = 0;
    for (int s = 0; s < nside * nside * nside - 1; s++) {
        offset[s + 1] = grid[s].np;
    }

    for (int s = 0; s < nside * nside * nside - 1; s++) {
        offset[s + 1] += offset[s];
    }

    for (int s = 0; s < nside * nside * nside; s++) {
        for (long long i = 0; i < grid[s].np; i++) {
            orderedPositions[3 * offset[s] + 3 * i] = grid[s].pos[3 * i];
            orderedPositions[3 * offset[s] + 3 * i + 1] = grid[s].pos[3 * i + 1];
            orderedPositions[3 * offset[s] + 3 * i + 2] = grid[s].pos[3 * i + 2];
        }
    }
}