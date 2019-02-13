/**
 * @file BRDFReader.hpp
 */


namespace ChefDevr {
    using namespace Eigen;
    using namespace std;


    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZ(const char *fileDirectory) {
        extract_brdfFilePaths(fileDirectory);

        const auto num_brdfs = brdf_filePaths.size();
        Matrix<Scalar> Z{num_brdfs, num_coefficientsBRDF};

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            Z.row(i) = read_brdf<Scalar>(brdf_filePaths[i].c_str());
        }

        return Z;
    }

    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZZt_centered(const char *fileDirectory, RowVector<Scalar> &meanBRDF) {
        extract_brdfFilePaths(fileDirectory);

        const auto num_brdfs = brdf_filePaths.size();
        Matrix<Scalar> ZZt_centered{num_brdfs, num_brdfs};
        meanBRDF.resize(num_coefficientsBRDF);

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector<Scalar> brdf_first = read_brdf<Scalar>(brdf_filePaths[i].c_str());

            ZZt_centered(i, i) = brdf_first.dot(brdf_first);
            meanBRDF += brdf_first;

            for (unsigned int j = i + 1; j < num_brdfs; ++j) {
                const RowVector<Scalar> brdf_second = read_brdf<Scalar>(brdf_filePaths[j].c_str());
                const Scalar coefficient = brdf_first.dot(brdf_second);
                ZZt_centered(i, j) = ZZt_centered(j, i) = coefficient;
            }
        }

        meanBRDF /= num_brdfs;

        RowVector<Scalar> brdf_brdfMean{num_brdfs};
        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector<Scalar> brdf = read_brdf<Scalar>(brdf_filePaths[i].c_str());
            brdf_brdfMean(i) = brdf.dot(meanBRDF);
        }

        ZZt_centered.colwise() -= brdf_brdfMean.transpose();
        ZZt_centered.rowwise() -= brdf_brdfMean;
        ZZt_centered = ZZt_centered.array() + meanBRDF.dot(meanBRDF);

        return ZZt_centered;
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

        red_value = (double)(brdf[index] * red_scale);
        green_value = (double)(brdf[index + stepBlue * 0.5] * green_scale);
        blue_value = (double)(brdf[index + stepBlue] * blue_scale);

        if (red_value < 0.0 || green_value < 0.0 || blue_value < 0.0) {
            throw BRDFReaderError{"Below horizon."};
        }
        
    }

    template <typename Scalar>
    RowVector<Scalar> BRDFReader::read_brdf(const char *filePath) {
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
        if (num_coefficients != num_coefficientsBRDF){
            std::fclose(file);
            throw BRDFReaderError{string{"Dimensions don't match : "} + to_string(num_coefficients) +
                                  " is not equal to " + to_string(num_coefficientsBRDF)};
        }

        RowVector<double> brdf_double{num_coefficients};
        numElementsRead = fread(brdf_double.data(), sizeof(double), num_coefficients, file);
        if (numElementsRead < num_coefficients) {
            std::fclose(file);
            throw BRDFReaderError{"The coefficients of the brdf has not been successfully read"};
        }
        
        std::fclose(file);

        // clamp negative values to zero
        brdf_double = brdf_double.cwiseMax(0.0);

        const RowVector<Scalar> brdf = brdf_double.template cast<Scalar>();
        return brdf;
    }

} // namespace ChefDevr
