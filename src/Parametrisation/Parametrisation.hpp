#include <Eigen/Eigenvalues>
#include <vector>
#include <numeric>

/**
 * @file Parametrisation.hpp
 */
#include <iostream>
namespace ChefDevr
{
    
template <typename Scalar>
void centerMat(Matrix<Scalar>& Z)
{
    Vector<Scalar> colMean(Z.rowwise().mean());
    unsigned int i(0);
    # pragma omp parallel for
    for(i=0; i<Z.cols(); ++i)
        Z.col(i) -= colMean;
}

template <typename Scalar>
void computeCovVector (
    Vector<Scalar>& cov_vector,
    const Matrix<Scalar>& X_reshaped,
    const Vector<Scalar>& coordRef,
    const unsigned int dim,
    const unsigned int nb_data)
{
    # pragma omp parallel for 
    for (unsigned int i = 0; i < nb_data; ++i){
        cov_vector[i] = covariance(coordRef, X_reshaped.col[i]);
    }
}

template <typename Scalar>
void BRDFReconstructor<Scalar>::reconstruct (Vector<Scalar>& brdf,
                                     const Vector<Scalar>& coord,
                                     const Scalar& mu)
{
    Vector<Scalar> cov_vector(nb_data);
    computeCovVector(cov_vector, X_reshaped, coord, dim, nb_data);
    brdf = cov_vector * K_minus1 * Z + mu;
}

} // namespace ChevDevr
