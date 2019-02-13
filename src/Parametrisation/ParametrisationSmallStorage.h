#ifndef PARAMETRISATION_WITH_SMALL_STORAGE__H
#define PARAMETRISATION_WITH_SMALL_STORAGE__H

#include "Parametrisation.h"

#include <BRDFReader/BRDFReader.h>


namespace ChefDevr {

    template <typename Scalar>
    class BRDFReconstructorSmallStorage : public BRDFReconstructor<Scalar>
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
        BRDFReconstructorSmallStorage (
                const Matrix<Scalar>& _K_minus1,
                const Vector<Scalar>& _X,
                const RowVector<Scalar>& _meanBRDF,
                const unsigned int _latentDim,
                BRDFReader &reader,
                const char *fileDirectory,
                const Scalar _mu = MU_DEFAULT,
                const Scalar _l = L_DEFAULT):
                BRDFReconstructor<Scalar>(_K_minus1, _X, _meanBRDF, _latentDim, _mu, _l),
                _K_minus1{_K_minus1},
                reader{reader},
                fileDirectory(fileDirectory)
        {}

        ~BRDFReconstructorSmallStorage() = default;


        /**
         * @brief Reconstructs a BRDF for latent space coordinates
         * @param brdf The brdf data vector to fill
         * @param coord Coordinates of the latent space point to recontruct as a BRDF
         * @return The BRDF data as a column vector
         */
        void reconstruct (RowVector<Scalar>& brdf, const Vector<Scalar>& coord) const override;

        Scalar reconstructionError (unsigned int brdfindex) const override;

    private:

        const Matrix<Scalar>& _K_minus1;

        BRDFReader &reader;

        const char *fileDirectory;

    };

} // ChefDevr

#include "ParametrisationSmallStorage.hpp"

#endif // PARAMETRISATION_WITH_SMALL_STORAGE__H