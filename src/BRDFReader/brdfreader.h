#ifndef BRDFREADER_H
#define BRDFREADER_H

#include "stdlib.h"
#include "math.h"
#include <cstdio>
#include <iostream>

#include "../types.h"

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define RED_SCALE (1.0/1500.0)
#define GREEN_SCALE (1.15/1500.0)
#define BLUE_SCALE (1.66/1500.0)
//#define M_PI	3.1415926535897932384626433832795

/**
* @brief This class is used to read all the BRDF references and to sample them
*/   
template <typename Scalar>
class BRDFReader
{
public:
    BRDFReader();

    /**
    * @brief Read all the brdfs stored in a given directory
    * @param fileDirectory the path of the directory where all the brdfs are stored
    * @return The set of brdfs, the Z matrix
    */
    Matrix<Scalar> create_brdfSet(const char *fileDirectory);

    /**
    * @brief Samples a BRDF
    * @param brdf the BRDF to be sampled
    * @param num_sampling the number of possible values for the angles that parametrizes the retained BRDF values
    * @return The sampled BRDF
    */
    Vector<Scalar> sample_brdf(const Vector<Scalar> &brdf, unsigned int num_sampling);

private:

	/**
	* @brief Read a BRDF from a file
	* @param filepath the path of the BRDF's file 
	* @return All the coefficients of a BRDF as a vector of scalars
	*
	* If the file is not found, returns an error
	*/
    Vector<Scalar> read_brdf(const char *filepath);

};

#endif // BRDFREADER_H
