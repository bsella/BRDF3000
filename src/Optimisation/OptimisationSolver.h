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
                
            Matrix<Scalar> K; /** Variance covariance matrix of latent variables : Mapping matrix */
            Matrix<Scalar> K_minus1; /** Inverse of K : Inverse mapping matrix */
            Vector<Scalar> X; /** Latent variables vector */
            Scalar cost; /** Value of the cost function for this solution */
        };
        
        OptimisationSolver(
            Scalar minStep,
            const Matrix<Scalar>& Z,
            const unsigned int dim);
        
        ~OptimisationSolver(){}
        
        /**
        * @brief Compute an optimised mapping from BRDFs space to a latent space
        * @return Optimisation result (cf OptiResult struct)
        * 
        * Uses Hook & Jeeves method to solve the optimisation
        */
        OptiResult computeOptimisation ();
        
    private:
        static constexpr Scalar reduceStep = .5f; /** Factor of step reduction */
        static constexpr Scalar mu = 0.0001f; /** mu constant that helps interpolating data while keeping good solution */
        static constexpr Scalar l = 1.0f; /** l constant defined in the research paper */
        const Scalar minStep; /** Step below wich solution is considered optimal */
        Scalar step; /** Value of the step for Hooke & Jeeves method */
        
        const unsigned int nb_data; /** Number of BRDFs in the Z matrix */
        const Matrix<Scalar>& Z; /** BRDFs data matrix */
        Matrix<Scalar> ZZt; /** Z*Ztransposed */
        
        unsigned int dim; /** Dimension of produced latent space */
        Vector<Scalar> X; /** Latent variables vector */
        Matrix<Scalar> K_minus1; /** Inverse of K : Inverse mapping matrix
        NB : We do not store K because the algorithm doesnt require it directly : we compute only columns of K when necessary instead -> covariance vectors */
        Scalar detK; /** Determinant of K */
        Scalar costval; /** Value of the cost function for this solution */
        
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
        Vector<Scalar> computeExplorationDisplacement ();

        /**
        * @brief Apply a previously computed displacement to each element of X
        * @return New modified vector
        */
        OptiResult computeLearnDisplacement (const OptiResult& optiRes);
        
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
        
        /**
         * @brief Computes the covariance vector (that should be stored in K)
         * of the lv_num'th latent variable
         * @param lv_num Number of the latent variable for the cov vector to be computed
         */
        Vector<Scalar> computeCovVector (unsigned int lv_num) const;
        
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
            return norm_x1_x2 < 0.00001 ? mu + exp_part : exp_part;
        }
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
