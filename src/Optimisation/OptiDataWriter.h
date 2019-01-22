#ifndef OPTI_DATA_WRITER__H
#define OPTI_DATA_WRITER__H

#include "../Parametrisation/types.h"
#include "../Optimisation/OptimisationSolver.h"

#include <string>
#include <vector>

namespace ChefDevr
{

    /**
    * @brief The purpose of this class is to write the data needed 
    * for a BRDFs reconstruction from latent space
    * and for a latent space representation
    */
    template <typename Scalar>
    class OptiDataWriter
    {
    public:
        OptiDataWriter (
            const Matrix<Scalar>& _Z,
            const typename OptimisationSolver<Scalar>::OptiResult& _optiRes):
            
            Z(_Z),
            optiRes(_optiRes){}
            
        ~OptiDataWriter (){}
        
        /**
         * @brief Writes the image of the latent space on the disk
         * @param path Path of the folder to write into
         * @param width Width of the image in number of pixels
         * @param height Height of the image in number of pixels
         */
        void writeAlbedoMap (
            const std::string& path,
            const unsigned width = 1024,
            const unsigned height = 1024);
        
        /**
         * @brief Writes the data that defines the parametrization (without Z)
         * @param path Path of the folder to write in=
         * @param brdfsFiles List of the BRDF files used for the parametrization
         * in the right order (consistent with Z matrix)
         */
        void writeParametrizationData (
            const std::string& path,
            const std::vector<std::string>& brdfsFiles);
        
    private:
        /** BRDFs data matrix */
        const Matrix<Scalar>& Z;
        
        /** Optimisation result that defines latent space */
        const typename OptimisationSolver<Scalar>::OptiResult& optiRes;
    };
} // ChefDevr

#include "OptiDataWriter.hpp"

#endif // OPTI_DATA_WRITER__H
