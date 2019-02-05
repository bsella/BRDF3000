#ifndef OPTI_DATA_WRITER__H
#define OPTI_DATA_WRITER__H

/**
* @file OptiDataWriter.h
* @brief Provides the functions to write latent space data needed for BRDFs reconstruction
* and latent space visualization
*/

#include "Parametrisation/types.h"
#include "Optimisation/OptimisationSolver.h"

#include <string>
#include <vector>


namespace ChefDevr
{
    /**
    * @brief Writes the image of the latent space on the disk
    * @param path Path of the folder to write into
    * @param Z Centered BRDF data matrix (BRDFs stored in row major)
    * @param X Latent variables vector
    * @param K_minus1 Inverse mapping matrix
    * @param meanBRDF Mean column of Z (before it was centered)
    * @param albedoSampling Sampling resolution in phi for albedo computation
    * @param width Width of the image in number of pixels
    * @param height Height of the image in number of pixels
    */
    template <typename Scalar>
    void writeAlbedoMap (
        const std::string& path,
        const Matrix<Scalar>& Z,
        const Vector<Scalar>& X,
        const Matrix<Scalar>& K_minus1,
        const Vector<Scalar>& meanBRDF,
        const unsigned int latentDim,
        const unsigned int albedoSampling = 4,
        const unsigned int width = 16,
        const unsigned int height = 16);
    
    /**
    * @brief Writes the data that defines the parametrisation (without Z)
    * @param path Path of the folder to write in
    * @param brdfsFilenames List of the BRDF files used for the parametrization
    * in the right order (consistent with Z matrix)
    * @param X Latent variables vector
    * @param K_minus1 Inverse mapping matrix
    * @param latentDim Dimension of latent space
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
    template <typename Scalar>
    void writeParametrisationData (
        const std::string& path,
        const std::vector<std::string>& brdfsFilenames,
        const Vector<Scalar>& X,
        const Matrix<Scalar>& K_minus1,
        const unsigned int latentDim);
} // ChefDevr

#include "OptiDataWriter.hpp"

#endif // OPTI_DATA_WRITER__H
