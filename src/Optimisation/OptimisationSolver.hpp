#include "OptimisationSolver.h"
#include <cmath>

namespace ChefDevr
{
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::computeOptimisation (
            const Matrix<Scalar>& Z,
            const int dim,
            float minStep)
    {
        return OptiResult();
    }
    
    template <typename Scalar>
    Scalar OptimisationSolver<Scalar>::cost(
        const Scalar& detK,
        const Matrix<Scalar>& K_minus1,
        const Matrix<Scalar>& ZZt,
        const unsigned int d
                                           )
    {
        Scalar trace(0.0);
        
        // Compute trace of K_minus1 * ZZt
        # pragma omp parallel for reduction(+:trace)
        for (unsigned int i(0); i < ZZt.cols(); ++i){
            trace += K_minus1.row(i).dot(ZZt.col(i));
        }
        
        return Scalar(0.5) * d * std::log(detK) + Scalar(0.5) * trace;
    }
    
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::computeExplorationDisplacement (
        const OptiResult& optiRes)
    {
        return OptiResult();
    }
    
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::computeLearnDisplacement (
        const OptiResult& optiRes)
    {
        return OptiResult();
    }
} // namespace ChefDevr
