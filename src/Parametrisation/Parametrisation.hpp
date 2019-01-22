#include <Eigen/Eigenvalues>
#include <vector>
#include <numeric>

namespace ChefDevr
{
    
template <typename Scalar>
void centerMat(Matrix<Scalar>& Z)
{
    // TODO
}

template <typename Scalar>
Vector<Scalar> computeCovVector (
    const Vector<Scalar>&X,
    const unsigned int lv_num,
    const unsigned int dim)
{
    const unsigned int nb_data = X.rows()/dim;
    auto X_reshaped(X.reshaped(dim, nb_data));
    auto latentRef(X_reshaped.col(lv_num)); // Latent variable that corresponds to the column of K 
    Vector<Scalar> cov_vector(dim);
    
    # pragma omp parallel for 
    for (unsigned int i(0); i < nb_data; ++i){
        cov_vector[i] = covariance(latentRef, X_reshaped.col[i]);
    }
    
    return cov_vector;
}

} // namespace ChevDevr
