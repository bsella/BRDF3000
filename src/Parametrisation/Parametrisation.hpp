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
    const Vector<Scalar>&X,
    const unsigned int lv_num,
    const unsigned int dim,
    const unsigned int nb_data)
{
    auto X_reshaped(X.reshaped(dim, nb_data));
    auto latentRef(X_reshaped.col(lv_num)); // Latent variable that corresponds to the column of K 
    
    # pragma omp parallel for 
    for (unsigned int i = 0; i < nb_data; ++i){
        cov_vector[i] = covariance(latentRef, X_reshaped.col[i]);
    }
}

} // namespace ChevDevr
