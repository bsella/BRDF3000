#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

/** 
 * @file OptimisationSolver.h
 */

#include "Parametrisation/types.h"
#include "../../tests/OptimisationTest.h"
#include <cmath>


namespace ChefDevr
{
    /**
     * @brief Class that solves the optimisation problem defined in the research paper:
     * A Versatile Parametrisation for Measured Materials Manifold.
     * Computes an optimised mapping from BRDFs space to a latent space
     * @tparam Scalar The type of scalar number to do computations with.
     * The precision of this type is crucial to reconstruct an accurate BRDF from the latent space
     */
    template <typename Scalar>
    class OptimisationSolver {
    public:
        
        /**
         * @brief Constructor
         * @param minStep Step below which solution is considered optimal
         * @param ZZt Z * (Z transposed) where Z is the centered BRDF data matrix (BRDFs stored in row major)
         * @param latentDim Dimension of optimised latent space
         */
        OptimisationSolver(
            long num_BRDFCoefficients,
            Scalar minStep,
            const Matrix<Scalar>& ZZt,
            unsigned int latentDim);
        
        ~OptimisationSolver() = default;

        /**
         * @brief Computes the optimized parametrisation of the BRDFs manifold.
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
         * @return A reference of the value of the cost function for the solution
         */
        inline const Scalar& getCostValue() const { return costval; }
        
    private:

        /**
         * @brief Initial value of the step
         */
        static constexpr Scalar step0 = .5f;
        
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
        const long nb_data;

        /**
         * @brief Number of coefficients of each BRDF stored inside the Z matrix
         */
        const long num_BRDFCoefficients;
        
        /**
         * @brief Z*Ztransposed
         */
        const Matrix<Scalar> ZZt;

        /**
         * @brief Dimension of produced latent space
         */
        unsigned int latentDim;

        /**
         * @brief Latent variables vector
         */
        Vector<Scalar> X;
        
        /**
         * @brief The movement vector of X that should improve the solution
         */
        Vector<Scalar> X_move;

        /**
         * @brief Inverse of K : Inverse mapping matrix
         *
         * We do not store K because the algorithm doesnt require it directly :
         * we compute only columns of K when necessary instead -> covariance vectors
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
         * @param cost value of the cost to fill
         * @param K_minus1 Inverse mapping
         * @param detK Determinant of the matrix K
         */
        void cost (Scalar& cost, const Matrix<Scalar>& K_minus1, const Scalar& detK) const;

        /**
         * @brief Updates the movement vector of X that improves the solution (X_move)
         * @return True if has moved, false otherwise
         * 
         * Adds and substracts step from each element
         * and reevaluate the cost function to check for better solutions
         */
        bool exploratoryMove ();

        /**
         * @brief Apply X_move to the latent variable vector X.
         * Updates new_X, new_K_minus1, new_detK accordingly
         * @param new_X Latent variables vector to apply the move on
         * @param new_K_minus1 New inverse mapping matrix K_minus1 to fill
         * @param new_detK New determinant of K to fill
         * @return true if effectively moved, false otherwise
         */
        bool patternMove (Vector<Scalar>& new_X, Matrix<Scalar>& new_K_minus1, Scalar& new_detK) const;
        
        /**
         * @brief Computes the new inverse matrix K_minus1 and the new determinant of K
         * using Sherman-Morisson formula
         * @param old_K_minus1 K_minus1 matrix before K had changed
         * @param new_K_minus1 K_minus1 matrix after K has changed (computed in this function)
         * @param old_detK Determinant of K before K had changed
         * @param new_detK Determinant of K after K has changed (computed in this function)
         * @param lv_num Number of the latent variable that has changed
         * @param diff_cov_vector The difference between
         * the column/row of K that has changed and the column/row before it had change.
         * Although this parameter is not const it remains unchanged when the function returns
         *
         * @pre old_detK != 0
         *
         * Computes new_K_minus_1 and new_detK from their previous values using Sherman-Morisson formula.
         * The formula for the inverse and determinant are applied twice : one for the row modification, and another one for the column modification.
         */
        void shermanMorissonUpdate (
            const Matrix<Scalar>& old_K_minus1,
            Matrix<Scalar>& new_K_minus1,
            const Scalar& old_detK,
            Scalar& new_detK,
            unsigned int lv_num,
            Vector<Scalar>& diff_cov_vector) const;
        
        /**
         * @brief Initializes the latent coordinates vector X by applying the PCA method
         * on the Z matrix and reducing its coordinates between [-1; 1]
         * @param K Variance-covariance matrix used to compute eigen vectors & values
         */
        void initX (const Matrix<Scalar>& K);

        /* ------------*/
        /* Friends */
        /* ------------*/

        friend OptimisationTest;
    };
} // namespace ChefDevr

#include "OptimisationSolver.hpp"

#endif // OPTIMISATIONSOLVER_H
