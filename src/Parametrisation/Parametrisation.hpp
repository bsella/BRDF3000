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
    for (unsigned int i = 0; i < nb_data; ++i){
        cov_vector[i] = covariance(latentRef, X_reshaped.col[i]);
    }
    
    return cov_vector;
}

template <typename Scalar>
Vector<Scalar> computePCA (const Matrix<Scalar>& Z, const unsigned int latentDim)
{
    unsigned char i;
    auto B(Z.transpose()*Z);
    
    // Use EigenSolver to compute eigen values and vectors
    auto solver(Eigen::EigenSolver<Matrix<Scalar>>(B, true));
    solver.compute();
    auto eigenValues(solver.eigenvalues());
    auto eigenVectors(solver.eigenvectors());
    
    // Use std::sort to sort an indices vector in the same way eigen values should be ordered
    // initialize original index locations
    std::vector<size_t> idx(eigenValues.size());
    std::iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in eigenValues
    std::sort(idx.begin(), idx.end(),
            [&eigenValues](size_t i1, size_t i2) {return eigenValues[i1] > eigenValues[i2];});
    
    // Build D the diagonal matrix with ordered eigen values in diagonal
    auto D(Matrix<Scalar>::Zero(latentDim, latentDim));
    for (i=0; i<latentDim; ++i)
    {
        // use sorted indices to retrieve eigen values in descending order
        D(i,i) = std::sqrt(eigenValues[idx[i]]);
    }
    
    // Build V the transposed matrix of ordered eigen vectors 
    Matrix<Scalar> V(latentDim, eigenVectors.cols());
    for (i=0; i<latentDim; ++i)
    {
        // use sorted indices to retrieve eigen vectors in descending order
        V.row(i) = eigenVectors.col(idx[i]);
    }
    
    // Latent variables for each BRDFs are now in column
    auto X(D*V);
    // Output X as column vector
    return X.reshaped();
}

} // namespace ChevDevr
