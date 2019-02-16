

namespace ChefDevr {

    template<typename Scalar, typename RScalar>
    void BRDFReconstructorWithZ<Scalar, RScalar>::reconstruct(RowVector<RScalar> &brdf,
                                                     const Vector <Scalar> &coord) const {
        RowVector<Scalar> cov_vector(BRDFReconstructor<Scalar, RScalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar, RScalar>::X, coord,
                                 BRDFReconstructor<Scalar, RScalar>::latentDim, BRDFReconstructor<Scalar, RScalar>::nb_data);
        brdf.noalias() = cov_vector * Km1Zc + BRDFReconstructor<Scalar, RScalar>::meanBRDF;
    }
    
    template<typename Scalar, typename RScalar>
    void BRDFReconstructorWithZ<Scalar, RScalar>::reconstruct(Eigen::Map<RowVector<RScalar>> &brdf,
                                                     const Vector <Scalar> &coord) const {
        RowVector<Scalar> cov_vector(BRDFReconstructor<Scalar, RScalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar, RScalar>::X, coord,
                                 BRDFReconstructor<Scalar, RScalar>::latentDim, BRDFReconstructor<Scalar, RScalar>::nb_data);
        brdf.noalias() = cov_vector * Km1Zc + BRDFReconstructor<Scalar, RScalar>::meanBRDF;
    }

    template<typename Scalar, typename RScalar>
    void BRDFReconstructorWithZ<Scalar, RScalar>::reconstructWithoutMean(RowVector <RScalar> &brdf,
                                                                const Vector <Scalar> &coord) const {
        RowVector<Scalar> cov_vector(BRDFReconstructor<Scalar, RScalar>::nb_data);
        computeCovVector<Scalar>(cov_vector.data(), BRDFReconstructor<Scalar, RScalar>::X, coord,
                                 BRDFReconstructor<Scalar, RScalar>::latentDim, BRDFReconstructor<Scalar, RScalar>::nb_data);
        brdf.noalias() = cov_vector * Km1Zc;
    }

    template<typename Scalar, typename RScalar>
    Scalar BRDFReconstructorWithZ<Scalar, RScalar>::reconstructionError(unsigned int brdfindex) const {
        if (brdfindex < 0 || brdfindex >= BRDFReconstructor<Scalar, RScalar>::nb_data) {
            std::cerr << "Given index for BRDF reconstruction is out of bounds !" << std::endl;
            return Scalar(-1);
        }

        RowVector<RScalar> reconstructed(Zcentered.cols());
        const Vector <Scalar> coord = BRDFReconstructor<Scalar, RScalar>::X.segment(brdfindex * BRDFReconstructor<Scalar, RScalar>::latentDim,BRDFReconstructor<Scalar, RScalar>::latentDim);
        
        reconstructWithoutMean(reconstructed, coord);
        const RowVector<Scalar> diff = reconstructed - Zcentered.row(brdfindex);

        return diff.dot(diff) / BRDFReconstructor<Scalar, RScalar>::meanBRDF.cols();
    }

}
