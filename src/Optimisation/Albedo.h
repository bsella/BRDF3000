#ifndef ALBEDO__H 
#define ALBEDO__H

#include "../Parametrisation/Parametrisation.h"
#include "../BRDFReader/BRDFReader.h"

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
        * @param brdf Resampled BRDF in the format defined in Methods & Algorithm report
        * @param albedo Albedo colour (r g b) to fill with result
        */
        static void computeAlbedo (const ResampledBRDF& brdf, Color& albedo);
    
    private:
        /**
        * @brief Computes the albedo of a BRDF in parallel with OpenMP
        * @param brdf Resampled BRDF in the format defined in Methods & Algorithm report
        * @param albedo Albedo colour (r g b) to fill with result
        */
        static void computeAlbedoOpenMP (const ResampledBRDF& brdf, Color& albedo);
        
        /**
        * @brief Computes the albedo of a BRDF in parallel with Nvidia Cuda
        * @param brdf Resampled BRDF in the format defined in Methods & Algorithm report
        * @param albedo Albedo colour (r g b) to fill with result
        */
        static void computeAlbedoCuda (const ResampledBRDF& brdf, Color& albedo);
    };

} // ChefDevr

#endif
