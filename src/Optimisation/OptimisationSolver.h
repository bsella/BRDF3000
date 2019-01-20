#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include "../types.h"
#include <cmath>

namespace ChefDevr
{
    template <typename Scalar>
    class OptimisationSolver {
    public:
        struct OptiResult
        {
            OptiResult(const OptiResult& other):
                K(other.K),
                K_minus1(other.K_minus1),
                X(other.X),
                cost(other.cost){}
                
            Matrix<Scalar> K; // Mapping matrix
            Matrix<Scalar> K_minus1; // Inverse mapping matrix
            Vector<Scalar> X; // Latent variables vector
            Scalar cost; // Value of the cost function for this solution
        };
        
        OptimisationSolver() = delete;
        ~OptimisationSolver() = delete;
        
        /**
        * @brief  Compute an optimised mapping from BRDFs space to a latent space
        * whose dimension is specified by dim
        * @param  Z BRDFs data matrix
        * @param  dim Dimension of latent space
        * @param  minStep Step beyond which convergence is considered being reached
        * @return Structure containing optimized mapping, inverse mapping and latent variables
        * 
        * Uses Hook & Jeeves method to solve the optimisation
        */
        static OptiResult computeOptimisation (
            const Matrix<Scalar>& Z,
            const int dim,
            float minStep);
        
    private:
        static constexpr Scalar reduceStep = .5f;
        static constexpr Scalar mu = 0.0001f;
        static constexpr Scalar l = 1.0f;
        
        /**
        * @brief Computes the cost of the inverse mapping defined by K_minus1
        * @param detK Determinant of the forward mapping matrix K
        * @param K_minus1 Inverse mapping matrix
        * @param ZZt Z*Ztransposed
        * @param nb_data The number of BRDFs in Z (Z.cols())
        */
        static Scalar cost(
            const Scalar& detK,
            const Matrix<Scalar>& K_minus1,
            const Matrix<Scalar>& ZZt,
            const unsigned int nb_data);

        /**
        * @brief Finds a displacement vector of X that improves the solution
        * @return Displacement vector
        * @param X Latent variables vector
        * @param K_minus1 Inverse mapping matrix
        * @param detK Determinant of the mapping matrix K
        * @param cost Value of the cost function for this solution
        * @param ZZt Z*Ztransposed
        * @param step Step
        * @param dim Latent space dimension
        * @return A displacement vector of X that improves the solution
        * 
        * Adds and substracts _reduceStep_ from each element
        * and reevaluate the cost function to check for better solutions
        */
        static Vector<Scalar> computeExplorationDisplacement (
            Vector<Scalar>& X,
            Matrix<Scalar>& K_minus1
            Scalar& detK,
            Scalar& cost,
            const Matrix<Scalar>& ZZt,
            const Scalar step,
            const unsigned char dim,
            const unsigned int nb_data);

        /**
        * @brief Apply a previously computed displacement to each element of X
        * @return New modified vector
        */
        static OptiResult computeLearnDisplacement (const OptiResult& optiRes);
        
        /**
         * @brief Updates the inverse matrix K_minus1 with Sherman-Morisson formula
         * @param K_minus1 Inverse mapping matrix
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         */
        static Matrix<Scalar> updateInverse (
            const Matrix<Scalar>& K_minus1,
            unsigned int lv_num,
            Vector<Scalar>& cov_vector);
        
        /**
         * @brief Updates the determinant of the matrix K with Sherman-Morisson formula
         * @param K_minus1 Inverse mapping matrix
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         */
        static Matrix<Scalar> updateDeterminant (
            const Matrix<Scalar>& K_minus1,
            unsigned int lv_num,
            Vector<Scalar>& cov_vector);
        
        /**
         * @brief Computes the covariance vector (usually stored in K)
         * of the lv_num'th latent variable
         * @param X Latent variables vector
         * @param lv_num Number of the latent variable for the cov vector to be computed
         */
        static Vector<Scalar> computeCovVector (
            const Vector<Scalar>& X,
            unsigned int lv_num);
        
        /**
         * @brief Covariance function given in the research paper :
         * A Versatile Parametrization for Measured Materials Manifold
         * @param x1 first latent variable
         * @param x2 second latent variable
         * @return Covariance value
         */
        static inline Scalar covariance(Vector<Scalar>& x1, Vector<Scalar>& x2)
        {
            const Scalar norm_x1_x2((x1-x2).norm());
            const Scalar exp_part(std::exp(-norm_x1_x2*norm_x1_x2/(Scalar(2)*l*l)));
            // dirac(x1-x2) <=> norm(x1-x2) == 0
            return norm_x1_2 < 0.00001 ? mu + exp_part : exp_part;
        }
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
