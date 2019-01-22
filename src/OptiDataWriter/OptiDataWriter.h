#ifndef OPTI_DATA_WRITER__H
#define OPTI_DATA_WRITER__H

#include "../types.h"
#include "../Optimisation/OptimisationSolver.h"

#include <string>
#include <vector>

namespace ChefDevr
{

    /**
    * @brief This class' purpose is to write the data needed for a BRDFs reconstruction from latent space
    * and latent space representation
    */
    template <typename Scalar>
    class OptiDataWriter
    {
    public:
        OptiDataWriter ()  = delete;
        ~OptiDataWriter () = delete;
        
        /**
         * @brief Writes the image of the latent space on the disk
         * @param path Path of the folder to write in
         * @param optiRes Optimisation result that defines latent space
         * @param Z BRDFs data matrix
         * @param width Width of the image in number of pixels
         * @param height Height of the image in number of pixels
         */
        static void writeAlbedoMap (
            const std::string& path,
            const OptiResult& optiRes,
            const Matrix<Scalar>& Z,
            const unsigned width = 1024,
            const unsigned height = 1024);
        
        /**
         * @brief Writes the data that defines the parametrization (without Z)
         * @param path Path of the folder to write in
         * @param optiRes Optimisation result that defines latent space
         * @param brdfsFiles List of the BRDF files used for the parametrization
         * in the right order (consistent with Z matrix)
         */
        static void writeParametrizationData (
            const std::string& path,
            const OptiResult& optiRes,
            const std::vector<std::string>& brdfsFiles);
        
    private:
        
    };
} // ChefDevr

#endif // OPTI_DATA_WRITER__H
