#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

/** 
 * @file OptimisationSolver.h
 */

#include "../Parametrisation/types.h"
#include <cmath>


namespace ChefDevr
{
    /**
     * @brief Class that solves the optimisation problem defined in the research paper:
     * A Versatile Parametrisation for Measured Materials Manifold
     * Computes an optimised mapping from BRDFs space to a latent space
     * @tparam Scalar The type of scalar number to do computations with.
     * The precision of this type is crucial to reconstruct an accurate BRDF from the latent space
     */
    template <typename Scalar>
    class OptimisationSolver {
    public:
        
        /**
         * @brief Constructor
         * @param minStep Step below wich solution is considered optimal
         */
        OptimisationSolver(
            Scalar minStep,
            Matrix<Scalar>& Z,
            unsigned int latentDim);
        
        ~OptimisationSolver(){}
        
        /**
        * @brief Computes the optimized parametrisation of the BRDFs manifold
        * Uses Hook & Jeeves method to solve the optimisation
        */
        void optimizeMapping ();
        
        /**
         * @return A reference of the inverse mapping matrix
         */
        inline const Matrix<Scalar>& getInverseMapping() const { return K_minus1; }
        
        /**
         * @return A reference of the latent variables vector
         */
        inline const Vector<Scalar>& getLatentVariables() const { return X; }
        
        /**
         * @return The value of the cost function for the solution
         */
        inline Scalar getCostValue() const { return costval; }
        
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
        unsigned int latentDim;

        /**
         * @brief Latent variables vector
         */
        Vector<Scalar> X;
        
        /**
         * @brief The displacement vector of X that should improve the solution
         */
        Vector<Scalar> X_move;

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
        * @brief Computes the cost of the solution defined by K_minus1
        * @param K_minus1 Inverse mapping
        * @param detK Determinant of the matrix K
        * @return Cost of the solution
        */
        Scalar cost (const Matrix<Scalar>& K_minus1, const Scalar& detK) const;

        /**
        * @brief Updates the displacement vector of X that improves the solution (X_move)
        * 
        * Adds and substracts step from each element
        * and reevaluate the cost function to check for better solutions
        */
        void exploratoryMove ();

        /**
        * @brief Apply X_move to the latent variable vector X
        * Updates new_X, new_K_minus1, new_detK accordingly
        * @param new_X Latent variables vector the apply the move on
        * @param new_K_minus1 New inverse mapping matrix K_minus1 to fill
        * @param new_detK New determinant of K to fill
        */
        void patternMove (Vector<Scalar>& new_X, Matrix<Scalar>& new_K_minus1, Scalar& new_detK) const;
        
        /**
         * @brief Computes the new inverse matrix K_minus1 with Sherman-Morisson formula
         * @param old_K_minus1 K_minus1 matrix before K had changed
         * @param new_K_minus1 K_minus1 matrix after K has changed
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         */
        void computeInverse (
            const Matrix<Scalar>& old_K_minus1,
            Matrix<Scalar>& new_K_minus1,
            unsigned int lv_num,
            Vector<Scalar>& cov_vector) const;
        
        /**
         * @brief Computes new the determinant of the matrix K with Sherman-Morisson formula
         * @param new_K_minus1 K_minus1 matrix after K has changed
         * @param lv_num Number of the latent variable that has changed
         * @param cov_vector The column of K that has changed
         * @return The new determinant of K
         */
        Scalar computeDeterminant (
            const Matrix<Scalar>& new_K_minus1,
            unsigned int lv_num,
            Vector<Scalar>& cov_vector) const;
        
        /**
        * @brief Initializes the latent coordinates vector X by applying the PCA method
        * on the Z matrix
        * 
        * Uses the Matusik method found in the paper
        * "A data-driven reflectance model"
        */
        void initX ();
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
