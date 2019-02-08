#ifndef ALBEDO__H 
#define ALBEDO__H

#include "Parametrisation/Parametrisation.h"
#include "BRDFReader/BRDFReader.h"

/**
 * @file Albedo.h
 */

namespace ChefDevr
{
    /**
     * @brief Provides albedo computation using the most efficient parallelization
     * solution made available by the material
     */
    class Albedo
    {
    public:
        /**
        * @brief Computes the albedo of a BRDF
        * @param brdf BRDF in the format defined in Methods & Algorithm report
        * @param[out] r red value for the albedo
        * @param[out] g green value for the albedo
        * @param[out] b blue value for the albedo
        * @param num_sampling the number of phi angles to sample
        */
        template <typename Scalar>
        static void computeAlbedo (
            const RowVector<Scalar>& brdf,
            double& r, double& g, double& b,
            const unsigned int num_sampling);
    
    private:
        /**
        * @brief Computes the albedo of a BRDF in parallel with OpenMP
        * @param brdf BRDF in the format defined in Methods & Algorithm report
        * @param[out] r red value for the albedo
        * @param[out] g green value for the albedo
        * @param[out] b blue value for the albedo
        * @param num_sampling the number of phi angles to sample
        */
        template <typename Scalar>
        static void computeAlbedoOpenMP (
            const RowVector<Scalar>& brdf,
            double& r, double& g, double& b,
            const unsigned int num_sampling);
        
        /**
        * @brief Computes the albedo of a BRDF in parallel with Nvidia Cuda
        * @param brdf Resampled BRDF in the format defined in Methods & Algorithm report
        * @param[out] r red value for the albedo
        * @param[out] g green value for the albedo
        * @param[out] b blue value for the albedo
        * @param num_sampling the number of phi angles to sample
        */
        template <typename Scalar>
        static void computeAlbedoCuda (
            const RowVector<Scalar>& brdf,
            double& r, double& g, double& b,
            const unsigned int num_sampling);
    };

} // ChefDevr

#include "Albedo.hpp"

#endif
