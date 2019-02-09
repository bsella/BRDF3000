/**
 * @file BRDFReader.hpp
 */

#include <experimental/filesystem>


namespace ChefDevr {
    using namespace Eigen;
    using namespace std;


    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZ(const char *fileDirectory) {
        using namespace std::experimental::filesystem;

        if (!is_directory(fileDirectory)) {
            throw BRDFReaderError{"The directory " + std::string{fileDirectory} + " does not exist"};
        }

        std::vector<std::string> list_filePaths;
        for(const path &filePath : directory_iterator(fileDirectory)) {
            list_filePaths.push_back(filePath);
            const auto filename = filePath.filename();
            brdf_filenames.push_back(filename.string());
        }

        const auto num_brdfs = list_filePaths.size();
        const unsigned int num_coefficientsBRDF = 3 * samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD / 2;
        Matrix<Scalar> Z{num_brdfs, num_coefficientsBRDF};

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            Z.row(i) = read_brdf<Scalar>(num_coefficientsBRDF, list_filePaths[i].c_str());
            // clamp negative values to zero
            Z.row(i) = Z.row(i).cwiseMax(Scalar(0));
        }

        return Z;
    }

    template <typename Scalar>
    void BRDFReader::lookup_brdf_val(const RowVector<Scalar> &brdf, double theta_in, double phi_in,
                                     double theta_out, double phi_out, double& red_value, double& green_value, double& blue_value) {
        // Convert to halfangle / difference angle coordinates
        double theta_half, phi_half, theta_diff, phi_diff;
        std_coords_to_half_diff_coords(theta_in, phi_in, theta_out, phi_out, theta_half, phi_half, theta_diff, phi_diff);

        // Find index.
        // Note that phi_half is ignored, since isotropic BRDFs are assumed
        const unsigned int index =  phi_diff_index(phi_diff) +
                                    theta_diff_index(theta_diff) * samplingResolution_phiD / 2 +
                                    theta_half_index(theta_half) * samplingResolution_phiD / 2 * samplingResolution_thetaD;

        const double stepBlue = samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD;

        red_value = brdf[index] * red_scale;
        green_value = brdf[index + stepBlue * 0.5] * green_scale;
        blue_value = brdf[index + stepBlue] * blue_scale;

        if (red_value < 0.0 || green_value < 0.0 || blue_value < 0.0) {
            throw BRDFReaderError{"Below horizon."};
        }
        
    }

    template <typename Scalar>
    RowVector<Scalar> BRDFReader::read_brdf(unsigned int num_coefficientsNeeded, const char *filePath) {
        FILE *file = fopen(filePath, "rb");
        if (!file) {
            throw BRDFReaderError{string{"The file "} + filePath + " could not have been opened"};
        }

        unsigned int dims[3];
        auto numElementsRead = fread(dims, sizeof(unsigned int), 3, file);
        if (numElementsRead < 3) {
            std::fclose(file);
            throw BRDFReaderError{"The dimensions of the brdf has not been successfully read"};
        }

        const unsigned int num_coefficients = dims[0] * dims[1] * dims[2] * 3;
        if (num_coefficients != num_coefficientsNeeded){
            std::fclose(file);
            throw BRDFReaderError{string{"Dimensions don't match : "} + to_string(num_coefficients) +
                                  " is not equal to " + to_string(num_coefficientsNeeded)};
        }

        RowVector<double> brdf_double{num_coefficients};
        numElementsRead = fread(brdf_double.data(), sizeof(double), num_coefficients, file);
        if (numElementsRead < num_coefficients) {
            std::fclose(file);
            throw BRDFReaderError{"The coefficients of the brdf has not been successfully read"};
        }
        
        std::fclose(file);

        const RowVector<Scalar> brdf = brdf_double.template cast<Scalar>();
        return brdf;
    }

} // namespace ChefDevr
