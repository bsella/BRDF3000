#include "OptiDataWriter.h"
#include "Albedo.h"
#include "../Parametrisation/Parametrisation.h"
#include "../../lib/bitmap/bitmap_image.hpp"
#include <iostream>
#include <cmath>

/**
 * @file OptiDataWriter.hpp
 */

namespace ChefDevr
{
    template <typename Scalar>
    void writeAlbedoMap (
        const std::string& path,
        const Matrix<Scalar>& Z,
        const Vector<Scalar>& X,
        const Matrix<Scalar>& K_minus1,
        const Vector<Scalar>& meanBRDF,
        const unsigned int latentDim,
        const unsigned width,
        const unsigned height)
    {
        if (latentDim != 2){
            std::cerr << "Cannot produce a map of a latent space whose dim is different than 2" << std::endl;
            return;
        }
        bitmap_image map(width, height);
        int color_res(std::pow(2, map.bytes_per_pixel()));
        Color albedo;
        
        BRDFReconstructor<Scalar> reconstructor(
            Z, K_minus1, X.reshaped(latentDim, Z.cols()), meanBRDF, latentDim);
        Scalar xstep(2/width), ystep(2/height);
        Vector<Scalar> brdf(Z.rows());
        Vector<Scalar> coord(2);
        coord << -1, -1;
        
        for (unsigned int pixx(0); pixx < width; ++pixx)
        {
            coord[0] += xstep;
            for (unsigned int pixy(0); pixy < height; ++pixy)
            {
                coord[1] += ystep;
                reconstructor.reconstruct(brdf, coord);
                computeAlbedo(brdf, albedo);
                map.set_pixel(pixx, pixy,
                             albedo.r*color_res, albedo.g*color_res, albedo.b*color_res);
            }
            coord[1] = -1;
        }
        map.save_image(path.c_str());
    }
} // namespace ChevDevr  
