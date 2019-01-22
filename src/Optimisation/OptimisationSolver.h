#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include "../types.h"
#include <cmath>

namespace ChefDevr
{
    /**
     * @brief Class that solves the optimisation problem defined in the research paper:
     * A Versatile Parametrization for Measured Materials Manifold
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
                
            /** Variance covariance matrix of latent variables : Mapping matrix */
            Matrix<Scalar> K;
            
            /** Inverse of K : Inverse mapping matrix */
            Matrix<Scalar> K_minus1;
            
            /** Latent variables vector */
            Vector<Scalar> X;
            
            /** Value of the cost function for this solution */
            Scalar cost;
        };
        
        OptimisationSolver(
            Scalar minStep,
            Matrix<Scalar>& Z,
            const unsigned int dim);
        
        ~OptimisationSolver(){}
        
        /**
        * @brief Computes the optimized parametrization of the BRDFs manifold
        * @return Optimisation result (cf OptiResult struct)
        * 
        * Uses Hook & Jeeves method to solve the optimisation
        */
        OptiResult optimizeMapping ();
        
    private:

        /** Factor of step reduction */
        static constexpr Scalar reduceStep = .5f;

        /** mu constant that helps interpolating data while keeping good solution */
        static constexpr Scalar mu = 0.0001f;

        /** l constant defined in the research paper */
        static constexpr Scalar l = 1.0f;

        /** Step below wich solution is considered optimal */
        const Scalar minStep;

        /** Value of the step for Hooke & Jeeves method */
        Scalar step;

        /** Number of BRDFs in the Z matrix */
        const unsigned int nb_data;

        /** BRDFs data matrix */
        const Matrix<Scalar>& Z;

        /** Z*Ztransposed */
        Matrix<Scalar> ZZt;
        

        /** Dimension of produced latent space */
        unsigned int dim;

        /** Latent variables vector */
        Vector<Scalar> X;

        /** Inverse of K : Inverse mapping matrix
        NB : We do not store K because the algorithm doesnt require it directly : 
        we compute only columns of K when necessary instead -> covariance vectors 
        */
        Matrix<Scalar> K_minus1;

        /** Determinant of K */
        Scalar detK;

        /** Value of the cost function for this solution */
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
        
        /**
         * @brief Computes the covariance vector (that should be stored in K)
         * of the lv_num'th latent variable
         * @param lv_num Number of the latent variable for the cov vector to be computed
         */
        Vector<Scalar> computeCovVector (unsigned int lv_num) const;
        
        /**
         * @brief Centers the Z BRDFs data matrix
         * Modifies the Z member variable 
         */
        void centerZ();
        
        /**
         * @brief Covariance function given in the research paper :
         * A Versatile Parametrization for Measured Materials Manifold
         * @param x1 first latent variable
         * @param x2 second latent variable
         * @return Covariance value
         */
        static inline Scalar covariance (Vector<Scalar>& x1, Vector<Scalar>& x2){
            const Scalar norm_x1_x2((x1-x2).norm());
            const Scalar exp_part(std::exp(-norm_x1_x2*norm_x1_x2/(Scalar(2)*l*l)));
            // dirac(x1-x2) <=> norm(x1-x2) == 0
            return norm_x1_x2 < 0.00001 ? mu + exp_part : exp_part;
        }
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
