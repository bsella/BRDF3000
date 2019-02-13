#include <experimental/filesystem>


namespace ChefDevr
{

    template<typename Scalar>
    void BRDFReconstructorSmallStorage<Scalar>::reconstruct(RowVector <Scalar> &brdf_reconstructed,
                                                            const Vector <Scalar> &coord) const {
        using namespace std::experimental::filesystem;

        RowVector <Scalar> cov_vector(BRDFReconstructor<Scalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar>::X, coord,
                                 BRDFReconstructor<Scalar>::latentDim, BRDFReconstructor<Scalar>::nb_data);
        const RowVector <Scalar> cov_Kminus1 = cov_vector * _K_minus1;

        if (!is_directory(fileDirectory)) {
            throw BRDFReader::BRDFReaderError{"The directory " + std::string{fileDirectory} + " does not exist"};
        }

        std::vector<std::string> list_filePaths;
        for (const path &filePath : directory_iterator(fileDirectory)) {
            list_filePaths.push_back(filePath);
        }

        const auto num_brdfs = list_filePaths.size();
        brdf_reconstructed = BRDFReconstructor<Scalar>::meanBRDF;

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector <Scalar> brdf = reader.read_brdf<Scalar>(list_filePaths[i].c_str());
            brdf_reconstructed += cov_Kminus1(i) * (brdf - BRDFReconstructor<Scalar>::meanBRDF);
        }
    }

    template<typename Scalar>
    Scalar BRDFReconstructorSmallStorage<Scalar>::reconstructionError(const unsigned int brdfindex) const {
        using namespace std::experimental::filesystem;

        if (brdfindex < 0 || brdfindex >= BRDFReconstructor<Scalar>::nb_data) {
            std::cerr << "Given index for BRDF reconstruction is out of bounds !" << std::endl;
            return Scalar(-1);
        }

        if (!is_directory(fileDirectory)) {
            throw BRDFReader::BRDFReaderError{"The directory " + std::string{fileDirectory} + " does not exist"};
        }

        std::vector<std::string> list_filePaths;
        for (const path &filePath : directory_iterator(fileDirectory)) {
            list_filePaths.push_back(filePath);
        }

        const long num_BRDFCoefficients =  BRDFReconstructor<Scalar>::meanBRDF.cols();
        RowVector <Scalar> reconstructed(num_BRDFCoefficients);
        const Vector <Scalar> coord = BRDFReconstructor<Scalar>::X.segment(brdfindex * BRDFReconstructor<Scalar>::latentDim,
                                                                           BRDFReconstructor<Scalar>::latentDim);
        reconstruct(reconstructed, coord);
        const RowVector <Scalar> brdf_gt = reader.read_brdf<Scalar>(list_filePaths[brdfindex].c_str());

        const RowVector <Scalar> diff((reconstructed - brdf_gt));

        return diff.dot(diff) / num_BRDFCoefficients;
    }

}