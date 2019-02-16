#include "MERLReader.h"

namespace ChefDevr {
    using namespace std;

    template<typename Scalar>
    void MERLReader::lookup_brdf_val(const RowVector<Scalar>& brdf, double theta_in, double phi_in,
                                     double theta_out, double phi_out, double &red_value, double &green_value,
                                     double &blue_value) {
        // Convert to halfangle / difference angle coordinates
        double theta_half, phi_half, theta_diff, phi_diff;
        std_coords_to_half_diff_coords(theta_in, phi_in, theta_out, phi_out, theta_half, phi_half, theta_diff,
                                       phi_diff);

        // Find index.
        // Note that phi_half is ignored, since isotropic BRDFs are assumed
        const unsigned int index = phi_diff_index(phi_diff) +
                                   theta_diff_index(theta_diff) * samplingResolution_phiD / 2 +
                                   theta_half_index(theta_half) * samplingResolution_phiD / 2 *
                                   samplingResolution_thetaD;

        const int stepBlue = samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD;

        red_value = (double) (brdf[index] * red_scale);
        green_value = (double) (brdf[index + stepBlue / 2] * green_scale);
        blue_value = (double) (brdf[index + stepBlue] * blue_scale);

        if (red_value < 0.0 || green_value < 0.0 || blue_value < 0.0) {
            throw MERLReaderError{"Below horizon."};
        }

    }

    template<typename Scalar>
    RowVector <Scalar> MERLReader::read_brdf(const char *filePath)
    {
        double* brdf_ptr;
        read_brdf(filePath, brdf_ptr);

        Eigen::Map<RowVector<double>> brdf_double{brdf_ptr, MERLReader::num_coefficientsBRDF};
        // clamp negative values to zero
        brdf_double = brdf_double.cwiseMax(0.0);

        const RowVector<Scalar> brdf = brdf_double.template cast<Scalar>();
        delete[] brdf_ptr;
        return brdf;
    }

} // ChefDevr
