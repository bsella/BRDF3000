#ifndef ALBEDO__H 
#define ALBEDO__H

#include "../Parametrisation/Parametrisation.h"

namespace ChefDevr
{
    struct Color
    {
        double r, g, b;
    }
    
    /**
     * @brief Computes the albedo of a BRDF
     * @return Albedo colour (r g b)
     */
    template <typename Scalar>
    Color computeAlbedo (Vector<Scalar>& brdf);
    
    /**
     * @brief Computes the albedo of a BRDF on a NVIDIA Cuda enabled graphics card
     * @return Albedo colour (r g b)
     */
    template <typename Scalar>
    Color computeAlbedoCuda (Vector<Scalar>& brdf);

} // ChefDevr

#include "Albedo.hpp"

#endif
