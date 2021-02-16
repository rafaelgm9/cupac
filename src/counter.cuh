//
// Created by rafa on 2/5/21.
//

#ifndef CUPAC_COUNTER_CUH
#define CUPAC_COUNTER_CUH

__global__
void doPairCount(long long np, const double *array, const long long *particlesInGrid, const long long *offset,
                 double boxsize, int nside, double minsep, double maxsep, int nbins, const double *rbins2,
                 unsigned long long *paircounts, const int *key, const int *iwrap);

#endif //CUPAC_COUNTER_CUH
