#include <stdlib.h>
#include <iostream>
#include <cstdio>

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360


// Read BRDF data
void read_brdf(const char *filename, double* &brdf, unsigned int num_coefficientsNeeded)
{
    FILE *f = fopen(filename, "rb");
    if (!f)
        exit(EXIT_FAILURE);

    int dims[3];
    fread(dims, sizeof(int), 3, f);
    unsigned int n = 3 * dims[0] * dims[1] * dims[2];
    if (n != num_coefficientsNeeded)
    {
        fprintf(stderr, "Dimensions don't match\n");
        fclose(f);
        exit(EXIT_FAILURE);
    }

    brdf = (double*) malloc (sizeof(double)*3*n);
    fread(brdf, sizeof(double), 3*n, f);

    fclose(f);
}


void write_brdf(const char *filePath, const double *brdf, unsigned int num_coefficientsNeeded)
{
    FILE *f = fopen(filePath, "w");
    if (!f)
        exit(EXIT_FAILURE);

    for (unsigned int i = 0; i < num_coefficientsNeeded; ++i) {
        fprintf(f, "%.20f ",brdf[i]);
    }

    fclose(f);
}


void write_setBRDFs(const char *filePath, double *brdf[], unsigned int num_brdfs, unsigned int num_coefficientsToWrite) {
    FILE *f = fopen(filePath, "w");
    if (!f)
        exit(EXIT_FAILURE);

    for (unsigned int index_coefficient = 0; index_coefficient < num_coefficientsToWrite; ++index_coefficient) {
        for (unsigned int index_brdf = 0; index_brdf < num_brdfs; ++index_brdf) {
            fprintf(f, "%.20f ",brdf[index_brdf][index_coefficient]);
        }
        fprintf(f, "\n");
    }

    fclose(f);
}


int main() {
    double *brdf;
    const unsigned int num_coefficientsToWrite = 100;
    const unsigned int num_coefficientsBRDF = 3 * BRDF_SAMPLING_RES_THETA_H *
                                            BRDF_SAMPLING_RES_THETA_D *
                                            BRDF_SAMPLING_RES_PHI_D / 2;

    std::cout << num_coefficientsBRDF << std::endl;

    read_brdf("brdfs/cherry-235.binary", brdf, num_coefficientsBRDF);
    write_brdf("brdf_output.txt", brdf, num_coefficientsToWrite);

    // To test the createZ function
    const unsigned int num_brdfs = 3;
    const std::string setBRDFpaths[3] = {"brdfs/cherry-235.binary", "brdfs/natural-209.binary", "brdfs/chrome.binary"};
    double *setBRDFs[3];

    for (unsigned int i = 0; i < num_brdfs; ++i) {
        read_brdf(setBRDFpaths[i].c_str(), setBRDFs[i], num_coefficientsBRDF);
    }

    write_setBRDFs("setBRDF_output.txt", setBRDFs, num_brdfs, num_coefficientsToWrite);

    /* Free */

    free(brdf);

    for (unsigned int i = 0; i < num_brdfs; ++i) {
        free(setBRDFs[i]);
    }

    exit(0);
}