#ifndef BRDFREADER_H
#define BRDFREADER_H

/**
 * @file BRDFReader.h
 */

#include "stdlib.h"
#include "math.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <stxxl/vector>

#include "../Parametrisation/types.h"

namespace ChefDevr {

    using namespace Eigen;
    using namespace std;
    
    /**
    * @brief This class is used to read all the BRDF references and to sample them
    */

    class BRDFReader {
    public:
        BRDFReader();
        ~BRDFReader(){}

        /**
        * @brief Read all the BRDFs stored in a given directory
        * @param fileDirectory the path of the directory where all the brdfs are stored
        * @return Z BRDFs data matrix where each column represents a BRDF
        *
        * Initializes the list of BRDFs paths in the order in which they were read.
        *
        * As the set of BRDFs can be heavy, we use the stxxl library
        * Thus, a problem is not likely to occur if the RAM is too small compared to the set of BRDFs
        * Indeed, in this case, the set of BRDFs is stored inside the disk
        */
        template<typename Scalar>
        Matrix<Scalar> createZ(const char *fileDirectory);

        /**
        * @brief Samples a BRDF
        * @param brdf the BRDF to be sampled
        * @param num_sampling the number of possible values for the angles
        * that parametrizes the retained BRDF values
        */
        template<typename Scalar>
        Vector<Scalar> sampleBRDF(const Vector<Scalar> &brdf, unsigned int num_sampling);

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
        *
        * As the set of BRDFs can be heavy, we use the stxxl library
        * Thus, a problem is not likely to occur if the RAM is too small compared to the set of BRDFs
        * Indeed, in this case, the set of BRDFs is stored inside the disk
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
        void lookup_brdf_val(const Vector<Scalar> &brdf, double theta_in, double phi_in, double theta_out, double phi_out, double& red_value, double& green_value, double& blue_value);

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
        void std_coords_to_half_diff_coords(double theta_in, double phi_in, double theta_out, double phi_out, double& theta_half, double& phi_half, double& theta_diff, double& phi_diff);

        /**
         * @brief Rotates a vector along an axis
         * @param vector the vector to be rotated
         * @param axis the axis along which the vector is rotated
         * @param angle the angle of the rotation
         * @return the rotated vector
         */
        Eigen::Vector3d rotate_vector(const Eigen::Vector3d &vector, const Eigen::Vector3d &axis, double angle);


        /**
         * @brief Computes a direction from angles
         * @param theta the theta angle
         * @param phi the phi angle
         * @return the direction
         */
        Eigen::Vector3d compute_direction (double theta, double phi);


        /**
         * @brief Lookup theta_half index
         * @param theta_half the angle corresponding to the index
         * @return the index
         * @pre theta_half is between 0 and PI / 2
         * @post the result is between 0 and 89
         *
         * This is a non-linear mapping!
         */
        unsigned int theta_half_index(double theta_half);

        /**
         * @brief Lookup theta_diff index
         * @param theta_diff the angle corresponding to the index
         * @return the index
         * @pre theta_diff is between 0 and PI / 2
         * @post the result is between 0 and 89
         */
        unsigned int theta_diff_index(double theta_diff);

        /**
         * @brief Lookup phi_diff index
         * @param phi_diff the angle corresponding to the index
         * @return the index
         * @pre phi_diff is between 0 and pi
         * @post the result is between 0 and 179
         */
        unsigned int phi_diff_index(double phi_diff);
    };
} // namespace ChefDevr

#include "BRDFReader.hpp"

#endif // BRDFREADER_H
