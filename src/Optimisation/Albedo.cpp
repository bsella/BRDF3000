#include "Albedo.h" 
#include <cmath>


/**
 * @brief Albedo.cpp
 */

namespace ChefDevr
{
    Color Albedo::computeAlbedo (const ResampledBRDF& brdf)
    {
        return computeAlbedoOpenMP(brdf);
    }
    
    Color Albedo::computeAlbedoOpenMP (const ResampledBRDF& brdf)
    {
        return Color(); 
    }
    
    Color Albedo::computeAlbedoCuda (const ResampledBRDF& brdf)
    {
        return Color();
    }
} // namespace ChefDevr
