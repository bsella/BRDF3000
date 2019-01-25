#ifndef OPTI_DATA_WRITER__H
#define OPTI_DATA_WRITER__H

/**
 * @file OptiDataWriter.h
 */

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
         * @brief Writes the data that defines the parametrisation (without Z)
         * @param path Path of the folder to write in
         * @param brdfsFiles List of the BRDF files used for the parametrization
         * in the right order (consistent with Z matrix)
         * 
         * <table>
         * <caption id="multi_row">File format</caption>
         * <tr><th>filename = paramtrzdata
         * <tr><td>propriete intellectuelle et commentaires
         * <tr><td>scalar_type
         * <tr><td>"inverse mapping matrix"
         * <tr><td>nb_rows nb_cols
         * <tr><td>K_minus1
         * <tr><td>"latent variables"
         * <tr><td>number_of_latent_variables
         * <tr><td>latent_space_dimension
         * <tr><td>brdf_filename x11 x12 ... x1d
         * <tr><td>brdf_filename x21 x22 ... x1d
         * <tr><td>...
         * <tr><td>brdf_filename xn1 xn2 ... xnd
         * </table>
         * 
         */
        void writeParametrisationData (
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
