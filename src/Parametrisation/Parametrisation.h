#ifndef PARAMETRISATION__H
#define PARAMETRISATION__H

/**
 * @file Parametrisation.h
 * @brief Functions and classes related to the BRDF space parametrisation
 * that are common to the Optimisation module and BRDF Explorer module
 */ 

#include <iostream>
#include "types.h"

#define MU_DEFAULT Scalar(0.0001f)
#define L_DEFAULT Scalar(1.0f)

namespace ChefDevr
{
    /**
     * @brief Class that allows BRDF reconstruction from latent space coordinates
     * @tparam Scalar The type of the values used to reconstruct a BRDF.
     * The precision of this type is crucial to reconstruct an accurate BRDF from the latent space.
     */
    template <typename Scalar>
    class BRDFReconstructor
    {
    public:
        /**
         * @brief Constructor of the class
         * @param _Zcentered Centered BRDFs data matrix (BRDFs stored in row major),
         * @param _K_minus1 Inverse mapping matrix
         * @param _X Latent variables vector
         * @param _meanBRDF The mean BRDF (mean of the rows of Z before it was centered)
         * @param _dim Dimension of the latent space
         * @param _mu Value of the mu constant that helps interpolation source data
         * @param _l Constant defined in the research paper
         */
        BRDFReconstructor (
            const Matrix<Scalar>& _Zcentered, 
            const Matrix<Scalar>& _K_minus1,
            const Vector<Scalar>& _X,
            const RowVector<Scalar>& _meanBRDF,
            const unsigned int _dim,
            const Scalar _mu = MU_DEFAULT,
            const Scalar _l = L_DEFAULT):
            
            Zcentered(_Zcentered),
            K_minus1(_K_minus1),
            X(_X),
            meanBRDF(_meanBRDF),
            dim(_dim),
            nb_data(Zcentered.rows()),
            mu(_mu),
            l(_l){}
        
        ~BRDFReconstructor(){}
        
        /**
         * @brief Reconstructs a BRDF for latent space coordinates
         * @param brdf The brdf data vector to fill
         * @param coord Coordinates of the latent space point to recontruct as a BRDF
         * @return The BRDF data as a column vector
         */
        void reconstruct (RowVector<Scalar>& brdf,
                          const Vector<Scalar>& coord) const;
        
        Scalar reconstructionError (unsigned int brdfindex) const;
        
    private:
        /** 
         * @brief BRDFs data matrix
         * 
         * Each row represents a BRDF
         */
        const Matrix<Scalar>& Zcentered;
        
        /** 
         * @brief Inverse mapping matrix
         */
        const Matrix<Scalar>& K_minus1;
        
        /** 
         * @brief Latent variables vector
         */
        const Vector<Scalar>& X;
        
        /**
         * @brief Mean column of Z (before it was centered)
         */
        const RowVector<Scalar>& meanBRDF;
        
        /** 
         * @brief Dimension of the latent space
         */
        const unsigned int dim;
        
        /**
         * @brief Number of BRDFs in Z
         */
        const unsigned int nb_data;
        
        /** 
         * @brief Value of the mu constant that helps interpolation source data
         */
        const Scalar mu;
        
        /** 
         * @brief l constant used in covariance function
         * 
         * See paper "A Versatile Parametrisation for Measured Materials Manifold"
         */
        const Scalar l;
        
        /**
         * @brief Reconstructs a BRDF for latent space coordinates without adding the mean
         * @param brdf The brdf data vector to fill
         * @param coord Coordinates of the latent space point to recontruct as a BRDF
         * @return The BRDF data as a column vector
         */
        void reconstructWithoutMean (RowVector<Scalar>& brdf,
                                    const Vector<Scalar>& coord) const;
        
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
        const Vector<Scalar>& x1,
        const Vector<Scalar>& x2,
        const Scalar mu = MU_DEFAULT,
        const Scalar l  = L_DEFAULT)
    {
        const Scalar sqnorm_x1_x2((x1-x2).squaredNorm());
        const Scalar exp_part(std::exp(-sqnorm_x1_x2/(Scalar(2)*l*l)));
        // dirac(x1-x2) <=> norm(x1-x2) == 0
        return sqnorm_x1_x2 < std::numeric_limits<Scalar>::epsilon() ? mu + exp_part : exp_part;
    }
    
    /**
     * @brief Centers matrix by sustracting mean to all columns
     * @param Z Matrix to center
     * @param meanBRDF Mean column of Z (filled in the function)
     * Z should be in column major
     */
    template <typename Scalar>
    void centerMat(Matrix<Scalar>& Z, RowVector<Scalar>& meanBRDF);
    
    /**
     * @brief Computes the covariance column vector for the coordRef coordinates variable
     * @param cov_vector The covariance column vector to fill
     * @param X_reshaped Latent variables matrix in the form such that each column is a latent variable with a number of rows equal to dim.
     * @param coordRef Coordinates to compare with every latent variable 
     * @param dim Dimension of latent space
     * @param nb_data Number of data 
     * @return Covariance column vector
     * 
     * Covariance is computed with the "covariance" function
     */
    template <typename Scalar>
    void computeCovVector (
        Scalar* cov_vector,
        const Vector<Scalar>& X,
        const Vector<Scalar>& coordRef,
        unsigned int dim,
        unsigned int nb_data);
        
} // ChefDevr

#include "Parametrisation.hpp"

#endif // PARAMETRISATION__H
