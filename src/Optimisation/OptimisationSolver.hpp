#include "OptimisationSolver.h"
#include "../Parametrisation/Parametrisation.h"

#include <cmath>


namespace ChefDevr
{
    template <typename Scalar>
    OptimisationSolver<Scalar>::OptimisationSolver(
        Scalar _minStep,
        Matrix<Scalar>& _Z,
        const unsigned int _dim) :
        
        minStep(_minStep),
        nb_data(_Z.cols()),
        Z(_Z),
        dim(_dim)
    {
        step = reduceStep;
        ZZt = Z*Z.transpose();
    }
    
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::optimizeMapping ()
    {
        X = computePCA<Scalar>(Z, dim);
        // Compute K
        // Compute detK
        // Compute K_minus1
        // Compute cost
        centerMat<Scalar>(Z);
        costval = cost();
        
        // while ...
        do
        {
            X = exploratoryMove();
            while(false)
            {
                X = patternMove();
            }
        }while(false/*change this of course*/);
        
        return OptiResult();
    }
    
    template <typename Scalar>
    Scalar OptimisationSolver<Scalar>::cost()
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
    Vector<Scalar> OptimisationSolver<Scalar>::exploratoryMove ()
    {
        const auto& nbcoefs(X.rows());
        Scalar new_costval, new_detK;
        Vector<Scalar> new_X(X);
        Vector<Scalar> cov_vector(nb_data), new_cov_vector(nb_data);
        Matrix<Scalar> new_K_minus1(K_minus1.rows(), K_minus1.cols());
        Vector<Scalar> X_move(nbcoefs);
        unsigned int lv_num;
        
        for (unsigned int i(0); i < nbcoefs;++i)
        {
            lv_num = i/dim;
            X_move[i] = step;
            new_X[i] += step;
            new_cov_vector = computeCovarianceVector(new_X, lv_num);
            
            // Update K_minus1 and detK with Sherman-Morisson formula
            new_K_minus1 = updateInverse(K_minus1, lv_num, new_cov_vector);
            new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
            // Update costval
            new_costval = cost(new_detK, new_K_minus1, ZZt, nb_data);
            
            if (new_costval > costval){
                X_move[i] = -step;
                new_X[i] -= Scalar(2)*step;
                new_cov_vector = computeCovVector(new_X, lv_num);
                
                // Update K_minus1 and detK with Sherman-Morisson formula
                new_K_minus1 = updateInverse(K_minus1, lv_num, new_cov_vector);
                new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
                // Update cost
                new_costval = cost(new_detK, new_K_minus1, ZZt, nb_data);
                
                if (new_costval > costval){
                    X_move[i] = Scalar(0);
                    new_X[i] += step;
                }
                else{
                    costval = new_costval;
                    // cost has changed -> keep new_K_minus1 and new_detK
                    K_minus1 = new_K_minus1;
                    detK = new_detK;
                }
            }
            else{
                costval = new_costval;
                // cost has changed -> keep new_K_minus1 and new_detK
                K_minus1 = new_K_minus1;
                detK = new_detK;
            }
        }
        // X = X+X_move; <=>
        X = new_X;
        return X_move;
    }
    
    template <typename Scalar>
    Matrix<Scalar> OptimisationSolver<Scalar>::computeInverse (unsigned int lv_num, Vector<Scalar>& cov_vector)const
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    Scalar OptimisationSolver<Scalar>::computeDeterminant (unsigned int lv_num, Vector<Scalar>& cov_vector)const
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    Vector<Scalar> OptimisationSolver<Scalar>::patternMove ()
    {
        return OptiResult();
    }
} // namespace ChefDevr
