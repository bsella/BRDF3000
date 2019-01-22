#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include "../Parametrisation/types.h"
#include <cmath>

namespace ChefDevr
{
    /**
     * @brief Class that solves the optimisation problem defined in the research paper:
     * A Versatile Parametrisation for Measured Materials Manifold
     * Computes an optimised mapping from BRDFs space to a latent space
     */
    template <typename Scalar>
    class OptimisationSolver {
    public:
        struct OptiResult{
            OptiResult(const OptiResult& other):
                K(other.K),
                K_minus1(other.K_minus1),
                X(other.X),
                cost(other.cost){}
                
            /**
             * @brief Variance covariance matrix of latent variables : Mapping matrix
             */
            Matrix<Scalar> K;
            
            /**
             * @brief Inverse of K : Inverse mapping matrix
             */
            Matrix<Scalar> K_minus1;
            
            /**
             * @brief Latent variables vector
             */
            Vector<Scalar> X;
            
            /**
             * @brief Value of the cost function for this solution
             */
            Scalar cost;
        };
        
        OptimisationSolver(
            Scalar minStep,
            Matrix<Scalar>& Z,
            unsigned int dim);
        
        ~OptimisationSolver(){}
        
        /**
        * @brief Computes the optimized parametrisation of the BRDFs manifold
        * @return Optimisation result (cf OptiResult struct)
        * 
        * Uses Hook & Jeeves method to solve the optimisation
        */
        OptiResult optimizeMapping ();
        
    private:

        /**
         * @brief Factor of step reduction
         */
        static constexpr Scalar reduceStep = .5f;

        /** @brief Step below wich solution is considered optimal */
        const Scalar minStep;

        /**
         * @brief Value of the step for Hooke & Jeeves method
         */
        Scalar step;

        /**
         * @brief Number of BRDFs in the Z matrix
         */
        const unsigned int nb_data;

        /**
         * @brief BRDFs data matrix
         */
        const Matrix<Scalar>& Z;

        /**
         * @brief Z*Ztransposed
         */
        Matrix<Scalar> ZZt;
        

        /**
         * @brief Dimension of produced latent space
         */
        unsigned int dim;

        /**
         * @brief Latent variables vector
         */
        Vector<Scalar> X;

        /**
         * @brief Inverse of K : Inverse mapping matrix
         *
         * We do not store K because the algorithm doesnt require it directly :
         * We compute only columns of K when necessary instead -> covariance vectors
        */
        Matrix<Scalar> K_minus1;

        /**
         * @brief Determinant of K
         */
        Scalar detK;

        /**
         * @brief Value of the cost function for this solution
         */
        Scalar costval;
        
        /**
        * @brief Computes the current cost of the solution
        * @return Current cost of the solution
        */
        Scalar cost ();

        /**
        * @brief Finds a displacement vector of X that improves the solution
        * @return A displacement vector of X that improves the solution
        * 
        * Adds and substracts _reduceStep_ from each element
        * and reevaluate the cost function to check for better solutions
        */
        Vector<Scalar> exploratoryMove ();

        /**
        * @brief Apply a previously computed displacement to each element of X
        * @return New modified vector
        */
        Vector<Scalar> patternMove ();
        
        /**
         * @brief Computes the inverse matrix K_minus1 with Sherman-Morisson formula
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         * @return The new K_minus1
         */
        Matrix<Scalar> computeInverse (unsigned int lv_num, Vector<Scalar>& cov_vector) const;
        
        /**
         * @brief Computes the determinant of the matrix K with Sherman-Morisson formula
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         * @return The new determinant of K
         */
        Scalar computeDeterminant (unsigned int lv_num, Vector<Scalar>& cov_vector) const;
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
