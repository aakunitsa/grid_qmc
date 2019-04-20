#ifndef POISSON_FORT_H
#define POISSON_FORT_H


extern "C" void construct_potential_(double *rho, double *pot);
extern "C" void initialize_poisson_(int *nrad, int *nang, int *iat); // requires the number of radial and angular points as well as atomic number
extern "C" void finalize_poisson_(void);
extern "C" double eri_fortran_(double *ci, double *cj, double *ck, double *cl);

#endif
