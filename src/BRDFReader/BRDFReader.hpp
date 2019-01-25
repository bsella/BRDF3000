/**
 * @file BRDFReader.hpp
 */

namespace ChefDevr {

    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZ(const char *fileDirectory) {

    }

    template <typename Scalar>
    Vector<Scalar> BRDFReader::sampleBRDF(const Vector<Scalar> &brdf, unsigned int num_sampling) {
        const unsigned int num_phi = 4 * num_sampling;
        const double stepTheta = 0.5 * M_PI / num_sampling;
        const double stepPhi = 2.0 * M_PI / num_phi;
        const unsigned int num_coefficients = num_sampling * num_phi * num_sampling * num_phi;
        Vector<Scalar> sampled_brdf{num_coefficients};
        unsigned int index = 0;

        for (unsigned int i = 0; i < num_sampling; i++){
            const double theta_in = i * stepTheta;

            for (unsigned int j = 0; j < num_phi; j++){
                const double phi_in = j * stepPhi;
                for (unsigned int k = 0; k < num_sampling; k++){
                    const double theta_out = k * stepTheta;

                    for (unsigned int l = 0; l < num_phi; l++){
                        const double phi_out = l * stepPhi;
                        lookup_brdf_val(brdf, theta_in, phi_in, theta_out, phi_out, sampled_brdf[index], sampled_brdf[index + 1], sampled_brdf[index + 2]);
                        index += 3;
                    }
                }
            }
        }

        return sampled_brdf;
    }

    template <typename Scalar>
    Vector<Scalar> BRDFReader::read_brdf(const char *filePath) {
        const FILE *file = fopen(filePath, "rb");
        if (!file) {
            throw BRDFReaderError{string{"The file "} + filePath + " could not have been opened"};
        }

        int dims[3];
        const unsigned int num_valuesNeeded = samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD / 2;
        fread(dims, sizeof(int), 3, file);
        const unsigned int num_values = dims[0] * dims[1] * dims[2];
        if (num_values != num_valuesNeeded){
            fclose(file);
            throw BRDFReaderError{string{"Dimensions don't match : "} + to_string(num_values) + " is not equal to " + to_string(num_valuesNeeded)};
        }

        const unsigned int num_coefficents = 3 * num_values;
        Vector<Scalar> brdf{num_coefficents};
        fread(brdf.data(), sizeof(Scalar), num_coefficents, file);

        fclose(file);

        return brdf;
    }

    template <typename Scalar>
    void BRDFReader::lookup_brdf_val(const Vector<Scalar> &brdf, double theta_in, double phi_in, double theta_out, double phi_out, double& red_value, double& green_value, double& blue_value) {
        // Convert to halfangle / difference angle coordinates
        double theta_half, phi_half, theta_diff, phi_diff;
        std_coords_to_half_diff_coords(theta_in, phi_in, theta_out, phi_out, theta_half, phi_half, theta_diff, phi_diff);

        // Find index.
        // Note that phi_half is ignored, since isotropic BRDFs are assumed
        const unsigned int index =  phi_diff_index(phi_diff) +
                                    theta_diff_index(theta_diff) * samplingResolution_phiD * 0.5 +
                                    theta_half_index(theta_half) * samplingResolution_phiD * 0.5 * samplingResolution_thetaD;

        const double stepBlue = samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD;

        red_value = brdf[index] * red_scale;
        green_value = brdf[index + stepBlue * 0.5] * green_scale;
        blue_value = brdf[index + stepBlue] * blue_scale;

        if (red_value < 0.0 || green_value < 0.0 || blue_value < 0.0) {
            throw BRDFReaderError{"Below horizon."};
        }
    }

} // namespace ChefDevr
