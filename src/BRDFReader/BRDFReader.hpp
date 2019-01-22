#include "BRDFReader.h"
/*
using namespace std;

namespace ChefDevr
{
	template <typename Scalar>
	BRDFReader<Scalar>::BRDFReader()
	{

	}

	template <typename Scalar>
	Matrix<Scalar> BRDFReader<Scalar>::create_brdfSet(const char *fileDirectory)
	{

	}

	template <typename Scalar>
	Vector<Scalar> BRDFReader<Scalar>::sample_brdf(const Vector<Scalar> &brdf, unsigned int num_sampling)
	{

	}

	template <typename Scalar>
	Vector<Scalar> BRDFReader<Scalar>::read_brdf(const char *filePath)
	{
		const FILE *file = fopen(filePath, "rb");
		if (!file) {
			throw BRDFReaderError{string{"The file "} + filePath + " could not have been opened"};
		}

		int dims[3];
		const unsigned int num_valuesNeeded = samplingResolution_thetaH *
								 samplingResolution_thetaD * samplingResolution_phiD * 0.5;

		fread(dims, sizeof(int), 3, file);
		const unsigned int num_values = dims[0] * dims[1] * dims[2];
		if (num_values != num_valuesNeeded){
				fclose(f);
				throw BRDFReaderError{"Dimensions don't match : " + num_values +
									  " is not equal to " + num_valuesNeeded};
		}

		const unsigned int num_coefficents = 3 * num_values;
		Vector<Scalar> brdf{num_coefficents};
		fread(brdf.data(), sizeof(Scalar), num_coefficents, file);

		fclose(file);

		return brdf;
	}

} // namespace ChefDevr

*/
