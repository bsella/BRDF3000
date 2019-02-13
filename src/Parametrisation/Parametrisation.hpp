//#include <vector>
//#include <numeric>

/**
 * @file Parametrisation.hpp
 */
namespace ChefDevr
{
    
template <typename Scalar>
void centerMat(Matrix<Scalar>& Z, RowVector<Scalar>& meanBRDF)
{
    meanBRDF.noalias() = Z.colwise().mean();
    Z.rowwise() -= meanBRDF;
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


} // namespace ChevDevr
