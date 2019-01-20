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
        const unsigned int nb_data
                                           )
    {
        Scalar trace(0);
        // Compute trace of K_minus1 * ZZt
        # pragma omp parallel for reduction(+:trace)
        for (unsigned int i(0); i < ZZt.cols(); ++i){
            trace += K_minus1.row(i).dot(ZZt.col(i));
        }
        return Scalar(0.5) * nb_data * std::log(detK) + Scalar(0.5) * trace;
    }
    
    template <typename Scalar>
    typename Vector<Scalar> OptimisationSolver<Scalar>::computeExplorationDisplacement (
        Vector<Scalar>& X,
        Matrix<Scalar>& K_minus1
        Scalar& detK,
        Scalar& cost,
        const Matrix<Scalar>& ZZt,
        const Scalar step,
        const unsigned char dim,
        const unsigned int nb_data)
    {
        const auto& nbcoefs(optiRes.X.rows());
        Scalar new_cost;
        Vector<Scalar> new_X(X);
        Vector<Scalar> cov_vector(X.rows()*0.5), new_cov_vector(X.rows()*0.5);
        Matrix<Scalar> new_K_minus1(K_minus1.rows(), K_minus1.cols());
        Vector<Scalar> X_move(nbcoefs);
        unsigned int lv_num;
        
        for (unsigned int i(0); i < nbcoefs;++i)
        {
            lv_num = i/dim;
            X_move[i] += step;
            new_X[i] += step;
            new_cov_vector = computeCovarianceVector(newX, lv_num)
            
            // Update K_minus1 and detK with Sherman-Morisson formula
            new_K_minus1 = updateInverse(K_minus1, lv_num, new_cov_vector);
            new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
            // Update cost
            new_cost = cost(new_detK, new_K_minus1, ZZt, nb_data);
            
            if (new_cost > cost){
                X_move[i] -= Scalar(2)*step;
                new_X[i] -= Scalar(2)*step;
                new_cov_vector = computeCovVector(newX, lv_num)
                
                // Update K_minus1 and detK with Sherman-Morisson formula
                new_K_minus1 = updateInverse(K_minus1, lv_num, new_cov_vector);
                new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
                // Update cost
                new_cost = cost(new_detK, new_K_minus1, ZZt, nb_data);
                
                if (new_cost > cost){
                    X_move[i] += step;
                    new_X[i] += step;
                }
                else{
                    cost = new_cost;
                    K_minus1 = new_K_minus1;
                    detK = new_detK;
                }
            }
            else{
                cost = new_cost;
                K_minus1 = new_K_minus1;
                detK = new_detK;
            }
        }
        X = X+X_move;
        return X_move;
    }
    
    template <typename Scalar>
    typename Matrix<Scalar> updateInverse (
        const Matrix<Scalar>& K_minus1,
        unsigned int lv_num,
        Vector<Scalar>& cov_vector)
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    typename Matrix<Scalar> updateDeterminant (
        const Matrix<Scalar>& K_minus1,
        unsigned int lv_num,
        Vector<Scalar>& cov_vector)
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::computeLearnDisplacement (
        const OptiResult& optiRes)
    {
        return OptiResult();
    }
    
} // namespace ChefDevr
