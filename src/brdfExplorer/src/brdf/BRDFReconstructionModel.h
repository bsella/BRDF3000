#ifndef BRDF_RECONSTRUCTION_MODEL_H
#define BRDF_RECONSTRUCTION_MODEL_H

#include "BRDFReconstructed.h

/**
 * @file BRDFReconstructionModel.h
 * @tparam Scalar The type of the values used to reconstruct a BRDF.
 * The precision of this type is crucial to reconstruct an accurate BRDF from the latent space.
 */


template <typename Scalar>
class BRDFReconstructionModel {
public:

    /**
     * @brief Constructs a BRDF from its latent coordinates
     * @param x the first coordinate of the BRDF in the latent space
     * @param y the second coordinate of the BRDF in the latent space
     * @return the reconstructed BRDF
     */
    BRDFReconstructed createBRDFFromLSCoord (float x, float y);

private:
    /**
     * @brief the matrix allowing to reconstruct a BRDF from its latent coordinates
     *
     * This is the inverse covariance matrix of the latent variables.
     * The latent variables are the BRDFs taken as reference in the latent space.
     */
    Matrix<Scalar> K_minus1;

    /**
     * @brief The BRDF matrix
     *
     * Each column contains the values of one BRDF
     */
    Matrix<Scalar> Z;
};

#endif
