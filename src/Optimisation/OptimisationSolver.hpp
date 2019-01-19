#include "OptimisationSolver.h"

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
        const Matrix<Scalar>& Z,
        const Matrix<Scalar>& Zt)
    {
        return 0;
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
