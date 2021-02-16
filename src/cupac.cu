//
// Created by rafa on 2/5/21.
//

#include "cupac.cuh"
#include <iostream>


unsigned long long *
autoPairCount(long long npart, const double *positions, double boxsize, int nside, double minsep, double maxsep, int nbins) {

    // Build Grid struct from positions
    Grid *grid = getGrid(boxsize, nside, npart, positions);

    // Get orderedPositions from Grid
    double *orderedPositions;
    long long *numParticlesInGrid, *offset;

    cudaMallocManaged(&orderedPositions, 3 * npart * sizeof(double));
    cudaMallocManaged(&numParticlesInGrid, nside * nside * nside * sizeof(long long));
    cudaMallocManaged(&offset, nside * nside * nside * sizeof(long long));

    gridToOrderedArray(grid, nside, orderedPositions, numParticlesInGrid, offset);

    // Free grid
    for (int i = 0; i < nside * nside * nside; i++) {
        if (grid[i].np > 0) {
            free(grid[i].pos);
        }
    }
    free(grid);

    // Set up rbinsSquared
    double *rbins, *rbinsSquared;

    cudaMallocManaged(&rbins, nbins * sizeof(double));
    cudaMallocManaged(&rbinsSquared, nbins * sizeof(double));

    for (int k = 0; k <= nbins; k++) {
        rbins[k] = (double) (pow(10, k * (log10(maxsep) - log10(minsep)) / nbins + log10(minsep)));
        rbinsSquared[k] = rbins[k] * rbins[k];
    }

    // Precompute key and iwrap for boundary conditions
    double agrid = boxsize / nside;
    int index_max = (int) (maxsep / agrid) + 1;
    int keysize = nside + 2 * index_max;
    int *key, *iwrap;

    cudaMallocManaged(&key, keysize * sizeof(int));
    cudaMallocManaged(&iwrap, keysize * sizeof(int));

    for (int ii = -index_max; ii <= nside + index_max; ii++) {
        if (ii < 0) {
            key[ii + index_max] = ii + nside;
            iwrap[ii + index_max] = -1;
        } else if (ii >= nside) {
            key[ii + index_max] = ii - nside;
            iwrap[ii + index_max] = 1;
        } else {
            key[ii + index_max] = ii;
            iwrap[ii + index_max] = 0;
        }
    }

    // Set up histogram
    unsigned long long *paircounts;

    cudaMallocManaged(&paircounts, nbins * sizeof(unsigned long long));

    for (int k = 0; k < nbins; k++) {
        paircounts[k] = 0;
    }

    // Get maximum potential gridSize and blockSize
    int gridSize, blockSize;
    cudaOccupancyMaxPotentialBlockSize(&gridSize, &blockSize, doPairCount, 0, 0);
    doPairCount<<<gridSize, blockSize>>>(npart, orderedPositions, numParticlesInGrid, offset, boxsize, nside, minsep,
                                         maxsep, nbins, rbinsSquared, paircounts, key, iwrap);

    cudaDeviceSynchronize();

    cudaFree(orderedPositions);
    cudaFree(numParticlesInGrid);
    cudaFree(offset);
    cudaFree(rbins);
    cudaFree(rbinsSquared);
    cudaFree(key);
    cudaFree(iwrap);

    return paircounts;
}