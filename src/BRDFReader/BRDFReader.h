#ifndef BRDFREADER_H
#define BRDFREADER_H

#include "stdlib.h"
#include "math.h"
#include <cstdio>
#include <iostream>

#include "../Parametrisation/types.h"

namespace ChefDevr
{

	/**
	* @brief This class is used to read all the BRDF references and to sample them
	*/   
	template <typename Scalar>
	class BRDFReader
	{
	public:
	    BRDFReader();

	    ~BRDFReader(){}

	    /**
	    * @brief Read all the brdfs stored in a given directory
	    * @param fileDirectory the path of the directory where all the brdfs are stored
	    * @return The set of brdfs, the Z matrix
	    *
	    * Initializes the list of BRDF paths in the order in which they were read
	    */
	    Matrix<Scalar> create_brdfSet(const char *fileDirectory);

	    /**
	    * @brief Samples a BRDF
	    * @param brdf the BRDF to be sampled
	    * @param num_sampling the number of possible values for the angles that parametrizes the retained BRDF values
	    * @return The sampled BRDF
	    */
	    Vector<Scalar> sample_brdf(const Vector<Scalar> &brdf, unsigned int num_sampling);

	    /**
	     * @return the list of BRDF paths in the order in which they were read
	     */
	    inline const std::vector<std::string>& getBRDFPaths () const {
	    	return brdf_filePaths;
	    }

	    class BRDFReaderError : public std::runtime_error {
	    public:
			explicit BRDFReaderError(const std::string &message_error) :
	    		runtime_error{message_error} {}
	    };

	private:

		/* ------------*/
		/* Attributes */
		/* ------------*/

		const static unsigned int samplingResolution_thetaH = 90;
		const static unsigned int samplingResolution_thetaD = 90;
		const static unsigned int samplingResolution_phiD = 360;

		constexpr static double red_scale = 1.0 / 1500.0;
		constexpr static double green_scale = 1.15 / 1500.0;
		constexpr static double blue_scale = 1.66 / 1500.0;

		/**
		 * @brief the list of BRDF paths in the order in which they were read
		 */
		std::vector<std::string> brdf_filePaths;

		/* ------------*/
		/* Functions */
		/* ------------*/

		/**
		* @brief Read a BRDF from a file
		* @param filePath the path of the BRDF's file 
		* @return All the coefficients of a BRDF as a vector of scalars
		*
		* If the file is not found, returns an error
		*/
	    Vector<Scalar> read_brdf(const char *filePath);

		/**
		 * @brief Extracts a color in a BRDF from a pair of incoming and outgoing angles
		 * @param[in] brdf the BRDF from which the color is extracted
		 * @param[in] theta_in incoming angle of theta
		 * @param[in] phi_in incoming angle of phi
		 * @param[in] theta_out outgoing angle of theta
		 * @param[in] phi_out outgoing angle of phi
		 * @param[out] red_value red channel of the extracted color
		 * @param[out] green_value green channel of the extracted color
		 * @param[out] blue_value blue channel of the extracted color
		 */
		void lookup_brdf_val(const Vector<Scalar> &brdf, double theta_in, double phi_in,
							 double theta_out, double phi_out,
							 double& red_value,double& green_value,double& blue_value);

		/**
		 * @brief Converts standard coordinates to half vector/difference vector coordinates
		 * @param[in] theta_in incoming angle of theta
		 * @param[in] phi_in incoming angle of phi
		 * @param[in] theta_out outgoing angle of theta
		 * @param[in] phi_out outgoing angle of phi
		 * @param[out] theta_half theta in half vector coordinates
		 * @param[out] phi_half phi in half vector coordinates
		 * @param[out] theta_diff theta in difference vector coordinates
		 * @param[out] phi_diff phi in difference vector coordinates
		 */
		void std_coords_to_half_diff_coords(double theta_in, double phi_in, double theta_out, double phi_out,
											double& theta_half,double& phi_half,double& theta_diff,double& phi_diff);

	};
} // namespace ChefDevr

#include "BRDFReader.hpp"

#endif // BRDFREADER_H
