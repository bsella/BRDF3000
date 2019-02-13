#ifndef BRDFREADER_H
#define BRDFREADER_H

/**
 * @file BRDFReader.h
 */

#include <cstdio>
#include <vector>
#include "Parametrisation/types.h"
#include "../tests/BRDFReaderTest.h"


namespace ChefDevr {
    
    /**
    * @brief This class is used to read all the BRDF references and to sample them
    */
    class BRDFReader {
    public:

        const static int samplingResolution_thetaH = 90;
        const static int samplingResolution_thetaD = 90;
        const static int samplingResolution_phiD = 360;

        BRDFReader();
        ~BRDFReader() = default;

        /**
        * @brief Read all the BRDFs stored in a given directory
        * @param fileDirectory the path of the directory where all the BRDFs are stored
        * @return Non-centered Z BRDFs data matrix where each row represents a BRDF
        *
        * Initializes the list of BRDFs filenames in the order in which they were read.
        */
        template<typename Scalar>
        Matrix<Scalar> createZ(const char *fileDirectory);

        /**
        * @brief Creates the centered ZZt matrix
        * @param fileDirectory the path of the directory where all the BRDFs are stored
        * @return the centered ZZt matrix
        *
        * Initializes the list of BRDFs filenames in the order in which they were read.
        *
        * As the set of BRDFs can be heavy, Z is not stored entirely inside the RAM.
        * Thus, a problem is not likely to occur if the RAM is too small compared to the set of BRDFs.
        */
        template <typename Scalar>
        Matrix<Scalar> createZZt_centered(const char *fileDirectory, RowVector<Scalar> &meanBRDF);


        /**
        * @brief Read a BRDF from a file
        * @param filePath the path of the BRDF's file
        * @return All the coefficients of a BRDF as a vector of scalars
        *
        * If the file is not found, returns an error
        */
        template<typename Scalar>
        RowVector<Scalar> read_brdf(const char *filePath);

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
        template <typename Scalar>
        static void lookup_brdf_val(const RowVector<Scalar>& brdf, double theta_in, double phi_in,
                             double theta_out, double phi_out, double& red_value, double& green_value, double& blue_value);

        /**
         * @return the list of BRDF filenames in the order in which they were read
         */
        inline const std::vector<std::string>& getBRDFFilenames () const {
            return brdf_filenames;
        }

        /**
        * @return the list of BRDF filepaths in the order in which they were read
        */
        inline const std::vector<std::string>& getBRDFFilePaths () const {
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

        constexpr static double red_scale = 1.0 / 1500.0;
        constexpr static double green_scale = 1.15 / 1500.0;
        constexpr static double blue_scale = 1.66 / 1500.0;

        /**
         * @brief Number of coefficients of each BRDF
         */
        const unsigned int num_coefficientsBRDF = 3 * samplingResolution_thetaH
                                                    * samplingResolution_thetaD
                                                    * samplingResolution_phiD / 2;

        /**
         * @brief the list of BRDF filenames in the order in which they were read
         */
        std::vector<std::string> brdf_filenames;

        /**
         * @brief the list of BRDF filePaths in the order in which they were read
         */
        std::vector<std::string> brdf_filePaths;

        /* ------------*/
        /* Functions */
        /* ------------*/

        /**
         * @brief Initializes the list of BRDFs filenames and the list of BRDF filePaths in the order in which they are read.
         * @param fileDirectory the path of the directory where all the BRDFs are stored
         */
        void extract_brdfFilePaths (const char *fileDirectory);


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
        static void std_coords_to_half_diff_coords(double theta_in, double phi_in, double theta_out, double phi_out,
                                            double& theta_half, double& phi_half, double& theta_diff, double& phi_diff);

        /**
         * @brief Rotates a vector along an axis
         * @param vector the vector to be rotated
         * @param axis the axis along which the vector is rotated
         * @param angle the angle of the rotation
         * @return the rotated vector
         */
        static Eigen::Vector3d rotate_vector(const Eigen::Vector3d &vector, const Eigen::Vector3d &axis, double angle);


        /**
         * @brief Computes a direction from angles
         * @param theta the theta angle
         * @param phi the phi angle
         * @return the direction
         */
        static Eigen::Vector3d compute_direction (double theta, double phi);


        /**
         * @brief Lookup theta_half index
         * @param theta_half the angle corresponding to the index
         * @return the index
         * @pre theta_half is between 0 and PI / 2
         * @post the result is between 0 and 89
         *
         * This is a non-linear mapping!
         */
        static unsigned int theta_half_index(double theta_half);

        /**
         * @brief Lookup theta_diff index
         * @param theta_diff the angle corresponding to the index
         * @return the index
         * @pre theta_diff is between 0 and PI / 2
         * @post the result is between 0 and 89
         */
        static unsigned int theta_diff_index(double theta_diff);

        /**
         * @brief Lookup phi_diff index
         * @param phi_diff the angle corresponding to the index
         * @return the index
         * @pre phi_diff is between 0 and pi
         * @post the result is between 0 and 179
         */
        static unsigned int phi_diff_index(double phi_diff);

        /* ------------*/
        /* Friends */
        /* ------------*/

        friend BRDFReaderTest;
    };
} // namespace ChefDevr

#include "BRDFReader.hpp"

#endif // BRDFREADER_H
