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


void write_brdf(const char *filename, const double * brdf, unsigned int num_coefficientsNeeded)
{
	FILE *f = fopen(filename, "w");
	if (!f)
		exit(EXIT_FAILURE);

	for (unsigned int i = 0; i < num_coefficientsNeeded; ++i) {
		fprintf(f, "%.20f ",brdf[i]);
	}

	fclose(f);
}


int main()
{
	double *brdf;
	unsigned int num_coefficientsNeeded = 3 * BRDF_SAMPLING_RES_THETA_H * 
										BRDF_SAMPLING_RES_THETA_D * 
										BRDF_SAMPLING_RES_PHI_D / 2;

	std::cout << num_coefficientsNeeded << std::endl;

	read_brdf("brdfs/cherry-235.binary", brdf, num_coefficientsNeeded);
	write_brdf("brdf_output.txt", brdf, 100);

	free(brdf);
	exit(0);
}