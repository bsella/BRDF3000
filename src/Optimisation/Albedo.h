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
        * @param resampled_brdf 
        * @return Albedo colour (r g b)
        */
        static Color computeAlbedo (const ResampledBRDF& brdf);
    
    private:
        /**
        * @brief Computes the albedo of a BRDF in parallel with OpenMP
        * @return Albedo colour (r g b)
        */
        static Color computeAlbedoOpenMP (const ResampledBRDF& brdf);
        
        /**
        * @brief Computes the albedo of a BRDF in parallel with Nvidia Cuda
        * @return Albedo colour (r g b)
        */
        static Color computeAlbedoCuda (const ResampledBRDF& brdf);
    };

} // ChefDevr

#endif
