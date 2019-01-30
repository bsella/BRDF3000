/**
 * @file BRDFReader.hpp
 */

//#include <stxxl/bits/io/syscall_file.h>
#include <experimental/filesystem>
//#include <stxxl/bits/stream/stream.h>


namespace ChefDevr {
    using namespace Eigen;
    using namespace std;


    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZ(const char *fileDirectory) {
        using namespace std::experimental::filesystem;
        //using StreamType = stxxl::stream::iterator2stream<typename std::vector<Scalar>::const_iterator>;

        if (!is_directory(fileDirectory)) {
            throw BRDFReaderError{"The directory " + std::string{fileDirectory} + " does not exist"};
        }

        std::vector<std::string> list_filePaths;
        for(const path &filePath : directory_iterator(fileDirectory)) {
            list_filePaths.push_back(filePath);
        }

        const auto num_brdfs = list_filePaths.size();
        const unsigned int num_coefficientsBRDF = 3 * samplingResolution_thetaH * samplingResolution_thetaD * samplingResolution_phiD / 2;
        Matrix<Scalar> Z{num_coefficientsBRDF, num_brdfs};

        //const auto num_coefficients = num_brdfs * num_coefficientsBRDF;
        //stxxl::vector<Scalar, 1> Z_stxxl;
        //Z_stxxl.reserve(num_coefficients);

        //auto Z_iterator = Z_stxxl.begin();
        for (unsigned int i = 0; i < num_brdfs; ++i) {
            Z.col(i) = read_brdf<Scalar>(num_coefficientsBRDF, list_filePaths[i].c_str());

            //StreamType input{brdf.begin(), brdf.end()};
            //Z_iterator = stxxl::stream::materialize(input, Z_iterator);
        }

        //Map< Matrix<Scalar> > Z(&Z_stxxl[0], num_coefficientsBRDF, num_brdfs);

        return Z;
    }

    template <typename Scalar>
    void BRDFReader::lookup_brdf_val(const Vector<Scalar> &brdf, double theta_in, double phi_in,
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
    Vector<Scalar> BRDFReader::read_brdf(unsigned int num_coefficientsNeeded, const char *filePath) {
        FILE *file = fopen(filePath, "rb");
        if (!file) {
            throw BRDFReaderError{string{"The file "} + filePath + " could not have been opened"};
        }

        unsigned int dims[3];
        std::fread(dims, sizeof(unsigned int), 3, file);
        const unsigned int num_coefficients = dims[0] * dims[1] * dims[2] * 3;
        if (num_coefficients != num_coefficientsNeeded){
            std::fclose(file);
            throw BRDFReaderError{string{"Dimensions don't match : "} + to_string(num_coefficients) + " is not equal to " + to_string(num_coefficientsNeeded)};
        }

        Vector<Scalar> brdf{num_coefficients};
        std::fread(brdf.data(), sizeof(Scalar), num_coefficients, file);

        std::fclose(file);

        return brdf;
    }

} // namespace ChefDevr
