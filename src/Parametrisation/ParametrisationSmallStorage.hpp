#include <experimental/filesystem>
#include "MERLReader.h"


namespace ChefDevr
{

    template<typename Scalar, typename RScalar>
    void BRDFReconstructorSmallStorage<Scalar, RScalar>::reconstruct(RowVector<RScalar> &brdf_reconstructed,
                                                            const Vector <Scalar> &coord) const {
        using namespace std::experimental::filesystem;

        RowVector <Scalar> cov_vector(BRDFReconstructor<Scalar, RScalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar, RScalar>::X, coord, BRDFReconstructor<Scalar, RScalar>::latentDim, BRDFReconstructor<Scalar, RScalar>::nb_data);
        
        const RowVector <Scalar> cov_Kminus1 = cov_vector * _K_minus1;

        const auto num_brdfs = _K_minus1.rows();
        brdf_reconstructed = BRDFReconstructor<Scalar, RScalar>::meanBRDF;

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector <Scalar> brdf = MERLReader::read_brdf<Scalar>(brdf_filePaths[i].c_str());
            brdf_reconstructed += cov_Kminus1(i) * (brdf - BRDFReconstructor<Scalar, RScalar>::meanBRDF);
        }
    }


    template<typename Scalar, typename RScalar>
    Scalar BRDFReconstructorSmallStorage<Scalar, RScalar>::reconstructionError(const unsigned int brdfindex) const {
        using namespace std::experimental::filesystem;

        if (brdfindex < 0 || brdfindex >= BRDFReconstructor<Scalar, RScalar>::nb_data) {
            std::cerr << "Given index for BRDF reconstruction is out of bounds !" << std::endl;
            return Scalar(-1);
        }

        const long num_BRDFCoefficients =  BRDFReconstructor<Scalar, RScalar>::meanBRDF.cols();
        RowVector <RScalar> reconstructed(num_BRDFCoefficients);
        const Vector <Scalar> coord = BRDFReconstructor<Scalar, RScalar>::X.segment(brdfindex * BRDFReconstructor<Scalar, RScalar>::latentDim, BRDFReconstructor<Scalar, RScalar>::latentDim);
        
        reconstruct(reconstructed, coord);
        const RowVector <Scalar> brdf_groundTruth = MERLReader::read_brdf<Scalar>(brdf_filePaths[brdfindex].c_str());

        const RowVector <Scalar> diff(reconstructed - brdf_groundTruth);

        return diff.dot(diff) / num_BRDFCoefficients;
    }

}
