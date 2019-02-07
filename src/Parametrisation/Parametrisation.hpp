#include <Eigen/Eigenvalues>
#include <vector>
#include <numeric>
/**
 * @file Parametrisation.hpp
 */
namespace ChefDevr
{
    
template <typename Scalar>
void centerMat(Matrix<Scalar>& Z, RowVector<Scalar>& meanBRDF)
{
    meanBRDF.noalias() = Z.colwise().mean();
    # pragma omp parallel for
    for(unsigned int i = 0; i<Z.rows(); ++i)
        Z.row(i) -= meanBRDF;
}

template <typename Scalar>
void computeCovVector (
    Scalar* cov_vector,
    const Vector<Scalar>& X,
    const Vector<Scalar>& coordRef,
    const unsigned int dim,
    const unsigned int nb_data)
{
    # pragma omp parallel for 
    for (unsigned int i = 0; i < nb_data; ++i){
        cov_vector[i] = covariance<Scalar>(coordRef, X.segment(i*dim,dim));
    }
}

template <typename Scalar>
void BRDFReconstructor<Scalar>::reconstruct (RowVector<Scalar>& brdf,
                                     const Vector<Scalar>& coord) const
{
    RowVector<Scalar> cov_vector(nb_data);
    computeCovVector<Scalar>(cov_vector.data(), X, coord, dim, nb_data);
    brdf.noalias() = cov_vector * Km1Zc + meanBRDF;
}

template <typename Scalar>
void BRDFReconstructor<Scalar>::reconstructWithoutMean (RowVector<Scalar>& brdf,
                                                        const Vector<Scalar>& coord) const
{
    RowVector<Scalar> cov_vector(nb_data);
    computeCovVector<Scalar>(cov_vector.data(), X, coord, dim, nb_data);
    brdf.noalias() = cov_vector * Km1Zc;
}

template <typename Scalar>
Scalar BRDFReconstructor<Scalar>::reconstructionError (const unsigned int brdfindex) const
{
    if (brdfindex < 0 || brdfindex >= nb_data){
        std::cerr << "Given index for BRDF reconstruction is out of bounds !" << std::endl;
        return Scalar(-1);
    }
    RowVector<Scalar> reconstructed(Zcentered.cols());
    reconstructWithoutMean(reconstructed, X.segment(brdfindex*dim,dim));
    //const Vector<Scalar> reconstructed(Zcentered.row(brdfindex).transpose()+meanBRDF);
    const RowVector<Scalar> diff((reconstructed - Zcentered.row(brdfindex)) );
    return (diff.dot(diff)/Zcentered.cols());
}

} // namespace ChevDevr
