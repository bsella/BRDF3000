#ifndef BRDFREADER_H
#define BRDFREADER_H

#include "stdlib.h"
#include "math.h"
#include <cstdio>
#include <iostream>

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define RED_SCALE (1.0/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE (1.66/1500.0)
//#define M_PI	3.1415926535897932384626433832795

class BRDFReader
{
public:
    BRDFReader();
    bool read_brdf(const char *filename, double* &brdf);

private:

};

#endif // BRDFREADER_H
