//
// Created by rafa on 2/5/21.
//

#ifndef CUPAC_CUPAC_CUH
#define CUPAC_CUPAC_CUH

#include "define.cuh"
#include "grid.cuh"
#include "counter.cuh"


extern "C" {
unsigned long long *
autoPairCount(long long npart, const double *positions, double boxsize, int nside, double minsep, double maxsep,
              int nbins);
}

#endif //CUPAC_CUPAC_CUH
