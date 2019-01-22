#ifndef PCASOLVER_H
#define PCASOLVER_H

#include "../types.h"

namespace ChefDevr
{
    /**
     * @brief This class is used to compute the PCA of the Z BRDFs data matrix
     * and initialize the X latent variable column vector
     */
    template <typename Scalar>
    class PCASolver{
    public:
        
        PCASolver () = delete;
        ~PCASolver () = delete;
        
        /**
        * @brief Generate a column vector of latent coordinates by applying the PCA method
        * on the Z matrix
        * @param Z The matrix of the BRDFs data 
        * @param latentDim Latent dimension (2 by default)
        * @return A matrix "reduced" in latent dimension
        * 
        * Uses the Matusik method found in the paper
        * "A data-driven reflectance model"
        */
        static Vector<Scalar> computePCA (const Matrix<Scalar>& Z, unsigned char latentDim=2);
    };
} // namespace ChefDevr

#endif // PCASOLVER_H
