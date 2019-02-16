#ifndef MERL_READER_H_
#define MERL_READER_H_

#include "types.h"


namespace ChefDevr {

    class MERLReader {
    public:
        constexpr static int samplingResolution_thetaH = 90;
        constexpr static int samplingResolution_thetaD = 90;
        constexpr static int samplingResolution_phiD = 360;

        /**
         * @brief Number of coefficients of each BRDF
         */
        constexpr static unsigned int num_coefficientsBRDF = 3 * samplingResolution_thetaH
                                                             * samplingResolution_thetaD
                                                             * samplingResolution_phiD / 2;


        MERLReader() = delete;

        ~MERLReader() = delete;

        /**
         * @brief Read a BRDF from a file
         * @param filePaths Path of brdf file
         * @param brdf Brdf data (allocated in this function)
         * @return All the coefficients of a BRDF as a vector of scalars
         *
         * If the file is not found, returns an error
         */
        static void read_brdf(const char *filePaths, double *&brdf);

        /**
         * @brief Read a BRDF from a file
         * @param index_brdf the index of the path of the BRDF's file in brdf_filePaths
         * @return All the coefficients of a BRDF as a vector of scalars
         *
         * If the file is not found, returns an error
         */
        template<typename Scalar>
        static RowVector<Scalar> read_brdf(const char *filePath);

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
        template<typename Scalar>
        static void lookup_brdf_val(const RowVector<Scalar>& brdf, double theta_in, double phi_in,
                                    double theta_out, double phi_out, double &red_value, double &green_value,
                                    double &blue_value);

    class MERLReaderError : public std::runtime_error {
    public:
        explicit MERLReaderError(const std::string& msg) :
                    std::runtime_error(msg){}
    };

    private:
        /* ------------*/
        /* Attributes */
        /* ------------*/

        constexpr static double red_scale = 1.0 / 1500.0;
        constexpr static double green_scale = 1.15 / 1500.0;
        constexpr static double blue_scale = 1.66 / 1500.0;

        /* ------------*/
        /* Functions */
        /* ------------*/


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
                                                   double &theta_half, double &phi_half, double &theta_diff,
                                                   double &phi_diff);

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
        static Eigen::Vector3d compute_direction(double theta, double phi);


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

    };

}

#include "MERLReader.hpp"

#endif // MERL_READER_H_
