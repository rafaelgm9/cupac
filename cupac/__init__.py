import numpy as np
import ctypes
import os
import glob


def pair_count(positions1, boxsize, nside, minsep, maxsep, nbins, positions2=None):
    lib_path = glob.glob(os.path.join(os.path.dirname(__file__), '../libcupac*.so'))[0]
    cupac_lib = ctypes.cdll.LoadLibrary(lib_path)

    cdouble_ptr = ctypes.POINTER(ctypes.c_double)
    cull_ptr = ctypes.POINTER(ctypes.c_ulonglong)

    if positions2 is not None:
        pair_count_run = cupac_lib.crossPairCount
        pair_count_run.argtypes = [ctypes.c_longlong, cdouble_ptr, ctypes.c_longlong, cdouble_ptr, ctypes.c_double,
                                   ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]
        pair_count_run.restype = cdouble_ptr

        positions1 = np.ascontiguousarray(positions1)
        positions2 = np.ascontiguousarray(positions2)
        positions1_ptr = positions1.ctypes.data_as(cdouble_ptr)
        positions2_ptr = positions2.ctypes.data_as(cdouble_ptr)

        print("doing pair count")
        paircounts_ptr = pair_count_run(positions1.shape[0], positions1_ptr, positions2.shape[0], positions2_ptr,
                                        boxsize, nside, minsep, maxsep, nbins)
        paircounts = np.ctypeslib.as_array(paircounts_ptr, shape=(nbins,))
    else:
        pair_count_run = cupac_lib.autoPairCount
        pair_count_run.argtypes = [ctypes.c_longlong, cdouble_ptr, ctypes.c_double, ctypes.c_int, ctypes.c_double,
                                   ctypes.c_double, ctypes.c_int]
        pair_count_run.restype = cull_ptr

        positions = np.ascontiguousarray(positions1)
        positions_ptr = positions.ctypes.data_as(cdouble_ptr)

        print("doing pair count")
        paircounts_ptr = pair_count_run(positions.shape[0], positions_ptr, boxsize, nside, minsep, maxsep, nbins)
        paircounts = np.ctypeslib.as_array(paircounts_ptr, shape=(nbins,))

    return paircounts
