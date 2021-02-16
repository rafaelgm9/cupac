//
// Created by rafa on 2/5/21.
//

#include "counter.cuh"

__global__
void
doPairCount(long long np, const double *array, const long long *particlesInGrid, const long long *offset,
            double boxsize, int nside, double minsep, double maxsep, const int nbins, const double *rbins2,
            unsigned long long *paircounts, const int *key, const int *iwrap) {
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    long long stride = blockDim.x * gridDim.x;

    double agrid = boxsize / nside;
    int index_max = (int) (maxsep / agrid) + 1;

    long long i, j;
    int i1, j1, k1;
    double x1, y1, z1;
    double x11, y11, z11;
    double dx, dy, dz;
    double r2;
    int iwrapx, iwrapy, iwrapz;
    int i2, j2, k2;
    int s;
    double factor = nside / boxsize;

//    __shared__  unsigned long long private_histogram[20];
//    if (threadIdx.x < 20) {
//        private_histogram[threadIdx.x] = 0;
//    }
//    __syncthreads();

    for (i = index; i < np; i += stride) {
        x1 = array[3 * i + 0];
        y1 = array[3 * i + 1];
        z1 = array[3 * i + 2];

        i1 = (int) (x1 * factor);
        j1 = (int) (y1 * factor);
        k1 = (int) (z1 * factor);

        i1 %= nside;
        j1 %= nside;
        k1 %= nside;

        for (int idx = -index_max; idx <= index_max; idx++) {
            i2 = i1 + idx;
            iwrapx = iwrap[i2 + index_max];
            i2 = key[i2 + index_max];

            for (int idy = -index_max; idy <= index_max; idy++) {
                j2 = j1 + idy;
                iwrapy = iwrap[j2 + index_max];
                j2 = key[j2 + index_max];

                for (int idz = -index_max; idz <= index_max; idz++) {
                    k2 = k1 + idz;
                    iwrapz = iwrap[k2 + index_max];
                    k2 = key[k2 + index_max];

                    s = nside * nside * k2 + nside * j2 + i2;

                    x11 = x1 - boxsize * iwrapx;
                    y11 = y1 - boxsize * iwrapy;
                    z11 = z1 - boxsize * iwrapz;

                    for (j = 0; j < particlesInGrid[s]; j += 1) {
                        dx = x11 - array[3 * offset[s] + 3 * j + 0];
                        dy = y11 - array[3 * offset[s] + 3 * j + 1];
                        dz = z11 - array[3 * offset[s] + 3 * j + 2];

                        r2 = dx * dx + dy * dy + dz * dz;

                        if (r2 > rbins2[nbins] || r2 < rbins2[0]) continue;

                        for (int k = nbins - 1; k >= 0; k--) {
                            if (r2 >= rbins2[k]) {
                                atomicAdd(&(paircounts[k]), 1);
//                                atomicAdd(&(private_histogram[k]), 1);
//                                paircounts[k] += 1;
                                break;
                            }
                        }
                    } //loop over neighbor grid2 particles
                } //loop over z neighbors
            } //loop over y neighbors
        } //loop over x neighbors
    } //loop over catalog1 particles

//    __syncthreads();
//
//    if (threadIdx.x < 20) {
//        atomicAdd(&(paircounts[threadIdx.x]), private_histogram[threadIdx.x]);
//    }

//    return paircounts;
}