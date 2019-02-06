#include "OptiDataWriter.h"
#include "Albedo.h"
#include "Parametrisation/Parametrisation.h"
#include "bitmap_image.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <experimental/filesystem>
#include <string>

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
        const RowVector<Scalar>& meanBRDF,
        const unsigned int latentDim,
        const unsigned int albedoSampling,
        const unsigned int width,
        const unsigned int height)
    {
        if (latentDim != 2){
            std::cerr << "Cannot produce a map of a latent space whose dim is different than 2" << std::endl;
            return;
        }
        bitmap_image map(width, height);
        int color_res(std::pow(2, map.bytes_per_pixel()));
        Color albedo;
        
        BRDFReconstructor<Scalar> reconstructor(
            Z, K_minus1, X, meanBRDF, latentDim);
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
                Albedo::computeAlbedo<Scalar>(brdf, albedo, albedoSampling);
                map.set_pixel(pixx, pixy,
                             albedo.r*color_res, albedo.g*color_res, albedo.b*color_res);
            }
            coord[1] = -1;
        }
        map.save_image(path.c_str());
    }
    
    template <typename Scalar>
    void writeParametrisationData (
        const std::string& path,
        const std::vector<std::string>& brdfsFilenames,
        const Vector<Scalar>& X,
        const Matrix<Scalar>& K_minus1,
        const unsigned int latentDim)
    {
        std::ofstream file(path);
        
        if (!file.is_open()){
            std::cerr << "Could not create file \"" << path<< "\"" << std::endl;
        }
        unsigned int nbLVars(X.rows()/latentDim);
        
        file << "# This file contains the results of an optimised parametrisation of the measured materials manifold" << std::endl;
        file << "# " << typeid(Scalar).name() << std::endl;
        file << "# inverse mapping matrix" << std::endl;
        file << K_minus1.rows() << " " << K_minus1.cols() << std::endl;
        file << K_minus1 << std::endl;
        file << "# latent variables" << std::endl;
        file << nbLVars << std::endl;
        file << latentDim << std::endl;

        for (unsigned int ivar(0); ivar < nbLVars; ++ivar)
        {
            file << brdfsFilenames[ivar] << " " <<  X.segment(ivar*latentDim, latentDim).transpose() << std::endl;
        }
    }
} // namespace ChevDevr  
