#include <Eigen/Eigenvalues>
#include <vector>
#include <numeric>

/**
 * @file Parametrisation.hpp
 */

namespace ChefDevr
{
    
template <typename Scalar>
void centerMat(Matrix<Scalar>& Z)
{
    // TODO
}

template <typename Scalar>
void computeCovVector (
    Vector<Scalar>& cov_vector,
    const Vector<Scalar>&X,
    const unsigned int lv_num,
    const unsigned int dim)
{
    const unsigned int nb_data = X.rows()/dim;
    auto X_reshaped(X.reshaped(dim, nb_data));
    auto latentRef(X_reshaped.col(lv_num)); // Latent variable that corresponds to the column of K 
    
    # pragma omp parallel for 
    for (unsigned int i = 0; i < nb_data; ++i){
        cov_vector[i] = covariance(latentRef, X_reshaped.col[i]);
    }
}

} // namespace ChevDevr
