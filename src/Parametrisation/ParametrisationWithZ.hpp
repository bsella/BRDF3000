

namespace ChefDevr {

    template<typename Scalar>
    void BRDFReconstructorWithZ<Scalar>::reconstruct(RowVector <Scalar> &brdf,
                                                     const Vector <Scalar> &coord) const {
        RowVector<Scalar> cov_vector(BRDFReconstructor<Scalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar>::X, coord,
                                 BRDFReconstructor<Scalar>::latentDim, BRDFReconstructor<Scalar>::nb_data);
        brdf.noalias() = cov_vector * Km1Zc + BRDFReconstructor<Scalar>::meanBRDF;
    }

    template<typename Scalar>
    void BRDFReconstructorWithZ<Scalar>::reconstructWithoutMean(RowVector <Scalar> &brdf,
                                                                const Vector <Scalar> &coord) const {
        RowVector<Scalar> cov_vector(BRDFReconstructor<Scalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar>::X, coord,
                                 BRDFReconstructor<Scalar>::latentDim, BRDFReconstructor<Scalar>::nb_data);
        brdf.noalias() = cov_vector * Km1Zc;
    }

    template<typename Scalar>
    Scalar BRDFReconstructorWithZ<Scalar>::reconstructionError(unsigned int brdfindex) const {
        if (brdfindex < 0 || brdfindex >= BRDFReconstructor<Scalar>::nb_data) {
            std::cerr << "Given index for BRDF reconstruction is out of bounds !" << std::endl;
            return Scalar(-1);
        }

        RowVector<Scalar> reconstructed(Zcentered.cols());
        const Vector <Scalar> coord = BRDFReconstructor<Scalar>::X.segment(brdfindex * BRDFReconstructor<Scalar>::latentDim,
                                                                  BRDFReconstructor<Scalar>::latentDim);
        reconstructWithoutMean(reconstructed, coord);
        const RowVector<Scalar> diff = reconstructed - Zcentered.row(brdfindex);

        return diff.dot(diff) / BRDFReconstructor<Scalar>::meanBRDF.cols();
    }

}
