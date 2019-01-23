#ifndef PARAMETRISATION__H
#define PARAMETRISATION__H

/**
 * @file Parametrisation.h
 * @brief Functions and classes that are common to the Optimisation module and BRDF Explorer module
 */ 


#include "./types.h"

namespace ChefDevr
{
    /**
    * @brief Class used for BRDF reconstruction from latent space coordinates
    */
    template <typename Scalar>
    class BRDFReconstructor
    {
    public:
        /**
        * @brief Constructor of the class
        * @param Z BRDFs data matrix,
        * @param K_minus1 Inverse mapping matrix
        * @param mu Value of the mu constant that helps interpolation source data
        * @param dim Dimension of the latent space
        */
        BRDFReconstructor (
            const Matrix<Scalar>& _Z, 
            const Matrix<Scalar>& _K_minus1,
            const Scalar _mu,
            const unsigned int dim):
            
            Z(_Z),
            K_minus1(_K_minus1),
            mu(_mu){}
        
        ~BRDFReconstructor(){}
        
        /**
        * @brief Reconstructs a BRDF for latent space coordinates
        * @param coord Coordinates of the latent space point to recontruct as a BRDF
        * @return The BRDF data as a column vector
        */
        Vector<Scalar> reconstruct (Vector<Scalar> coord);
        
    private:
        /** BRDFs data matrix */
        const Matrix<Scalar>& Z;
        
        /** Inverse mapping matrix */
        const Matrix<Scalar>& K_minus1;
        
        /** Value of the mu constant that helps interpolation source data */
        const Scalar mu;
        
        /** Dimension of the latent space */
        const unsigned int dim;
        
    };
    
    /**
    * @brief Covariance function given in the research paper :
    * A Versatile Parametrization for Measured Materials Manifold
    * @param x1 First latent variable
    * @param x2 Second latent variable
    * @param mu The constant that helps interpolating data while keeping good solution
    * @param l Constant defined in the research paper
    * @return Covariance value
    */
    template <typename Scalar>
    inline Scalar covariance (
        Vector<Scalar>& x1,
        Vector<Scalar>& x2,
        const Scalar mu = 0.0001f,
        const Scalar l  = 1.0f)
    {
        const Scalar norm_x1_x2((x1-x2).norm());
        const Scalar exp_part(std::exp(-norm_x1_x2*norm_x1_x2/(Scalar(2)*l*l)));
        // dirac(x1-x2) <=> norm(x1-x2) == 0
        return norm_x1_x2 < Scalar(0.00001f) ? mu + exp_part : exp_part;
    }
    
    /**
     * @brief Centers matrix by sustracting mean to all columns
     * @param Z Matrix to center
     * Z should be in column major
     */
    template <typename Scalar>
    void centerMat(Matrix<Scalar>& Z);
    
    /**
    * @brief Computes the covariance column vector for the lv_num'th latent variable
    * @param cov_vector The covariance column vector to fill
    * @param lv_num Number of the latent variable for the cov vector to be computed
    * @param nb_data Number of data 
    * @return Covariance column vector
    * 
    * Covariance is computed with the "covariance" function
    */
    template <typename Scalar>
    void computeCovVector (
        Vector<Scalar>& cov_vector,
        const Vector<Scalar>&X,
        const unsigned int lv_num,
        const unsigned int dim);
    
} // ChefDevr

#include "Parametrisation.hpp"

#endif // PARAMETRISATION__H
