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
    * @brief Simple (red, green, blue) color structure
    */
    struct Color
    {
        double r, g, b;
    };
    
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
        * @param num_sampling the number of phi angles to sample
        * @param albedo Albedo colour (r g b) to fill with result
        */
        template <typename Scalar>
        static void computeAlbedo (const Vector<Scalar>& brdf, Color& albedo, unsigned int num_sampling);
    
    private:
        /**
        * @brief Computes the albedo of a BRDF in parallel with OpenMP
        * @param brdf BRDF in the format defined in Methods & Algorithm report
        * @param albedo Albedo colour (r g b) to fill with result
        * @param num_sampling the number of phi angles to sample
        */
        template <typename Scalar>
        static void computeAlbedoOpenMP (const Vector<Scalar>& brdf, Color& albedo, unsigned int num_sampling);
        
        /**
        * @brief Computes the albedo of a BRDF in parallel with Nvidia Cuda
        * @param brdf Resampled BRDF in the format defined in Methods & Algorithm report
        * @param albedo Albedo colour (r g b) to fill with result
        */
        template <typename Scalar>
        static void computeAlbedoCuda (const Vector<Scalar>& brdf, Color& albedo);
    };

} // ChefDevr

#include "Albedo.hpp"

#endif
