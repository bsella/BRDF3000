#ifndef PARAMETRISATION_WITH_SMALL_STORAGE__H
#define PARAMETRISATION_WITH_SMALL_STORAGE__H

#include "Parametrisation.h"

#include <BRDFReader/BRDFReader.h>


/**
 * @file ParametrisationSmallStorage.h
 * @brief Does the BRDF space parametrisation using a small amount of Ram
 */
namespace ChefDevr {

    template <typename Scalar>
    class BRDFReconstructorSmallStorage : public BRDFReconstructor<Scalar>
    {
    public:
        /**
         * @brief Constructor of the class
         * @param _K_minus1 Inverse mapping matrix
         * @param _X Latent variables vector
         * @param _meanBRDF The mean BRDF (mean of the rows of Z before it was centered)
         * @param _latentDim Dimension of the latent space
         * @param reader The reader of BRDFs
         * @param _mu Value of the mu constant that helps interpolation source data
         * @param _l Constant defined in the research paper
         */
        BRDFReconstructorSmallStorage (
                const Matrix<Scalar>& _K_minus1,
                const Vector<Scalar>& _X,
                const RowVector<Scalar>& _meanBRDF,
                const unsigned int _latentDim,
                BRDFReader &reader,
                const Scalar _mu = MU_DEFAULT,
                const Scalar _l = L_DEFAULT):
                BRDFReconstructor<Scalar>(_K_minus1, _X, _meanBRDF, _latentDim, _mu, _l),
                _K_minus1{_K_minus1},
                reader{reader}
        {}

        ~BRDFReconstructorSmallStorage() = default;


        /**
         * @brief Reconstructs a BRDF from its latent space coordinates
         * @param[out] brdf The brdf data vector to fill
         * @param[in] coord Coordinates of the latent space point from which a BRDF is reconstructed
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
         * @brief Inverse of K : Inverse mapping matrix
        */
        const Matrix<Scalar>& _K_minus1;

        /**
         * @brief The reader of brdfs
         */
        BRDFReader &reader;

    };

} // ChefDevr

#include "ParametrisationSmallStorage.hpp"

#endif // PARAMETRISATION_WITH_SMALL_STORAGE__H