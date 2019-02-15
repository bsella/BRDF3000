#ifndef PARAMETRISATION_WITH_Z__H
#define PARAMETRISATION_WITH_Z__H

#include "Parametrisation.h"


/**
 * @file ParametrisationWithZ.h
 * @brief Does the BRDF space parametrisation
 */
namespace ChefDevr {

    template <typename Scalar>
    class BRDFReconstructorWithZ : public BRDFReconstructor<Scalar>
    {
    public:
        /**
         * @brief Constructor of the class
         * @param _Zcentered Centered BRDFs data matrix (BRDFs stored in row major),
         * @param _K_minus1 Inverse mapping matrix
         * @param _X Latent variables vector
         * @param _meanBRDF The mean BRDF (mean of the rows of Z before it was centered)
         * @param _latentDim Dimension of the latent space
         * @param _mu Value of the mu constant that helps interpolation source data
         * @param _l Constant defined in the research paper
         */
        BRDFReconstructorWithZ (
                const Matrix<Scalar>& _Zcentered,
                const Matrix<Scalar>& _K_minus1,
                const Vector<Scalar>& _X,
                const RowVector<Scalar>& _meanBRDF,
                const unsigned int _latentDim,
                const Scalar _mu = MU_DEFAULT,
                const Scalar _l = L_DEFAULT):
                BRDFReconstructor<Scalar>(_K_minus1, _X, _meanBRDF, _latentDim, _mu, _l),
                Zcentered(_Zcentered),
                Km1Zc(_K_minus1*_Zcentered)
        {}

        ~BRDFReconstructorWithZ() = default;


        /**
         * @brief Reconstructs a BRDF from its latent space coordinates
         * @param brdf The brdf data vector to fill
         * @param coord Coordinates of the latent space point to recontruct as a BRDF
         * @return The BRDF data as a row vector
         */
        void reconstruct (RowVector<Scalar>& brdf, const Vector<Scalar>& coord) const override;

        /**
         * @brief Computes the error between a reference brdf and this brdf reconstructed from its latent coordinates
         * @param brdfindex : The index of the brdf in the list of brdfs read to construct Z
         * @return the mean square error between a reference brdf and its reconstruction
         */
        Scalar reconstructionError (unsigned int brdfindex) const override;

    private:

        /**
         * @brief Centered BRDFs data matrix (BRDFs stored in row major)
         */
        const Matrix<Scalar>& Zcentered;

        /**
         * @brief K_minus1 times Z centered
         */
        const Matrix<Scalar> Km1Zc;

        /**
         * @brief Reconstructs a BRDF for latent space coordinates without adding the mean
         * @param brdf The brdf data vector to fill
         * @param coord Coordinates of the latent space point to recontruct as a BRDF
         * @return The BRDF data as a column vector
         */
        void reconstructWithoutMean (RowVector<Scalar>& brdf,
                                     const Vector<Scalar>& coord) const;

    };

} // ChefDevr

#include "ParametrisationWithZ.hpp"

#endif // PARAMETRISATION_WITH_Z__H
