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
    void progressBar(double progress)
    {
        const int barWidth = 70;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
    }

    /*
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
        const unsigned int height,
        const double latentWidth,
        const double latentHeight)
    {
        BRDFReconstructor<Scalar> reconstructor(
            Z, K_minus1, X, meanBRDF, latentDim);
        writeAlbedoMap(path, reconstructor, albedoSampling, width, height);
    }
    */
    
    template <typename Scalar>
    void writeAlbedoMap (
        const std::string& path,
        const BRDFReconstructor<Scalar>* reconstructor,
        const unsigned int albedoSampling,
        const unsigned int width,
        const unsigned int height,
        const double latentWidth,
        const double latentHeight)
    {
         if (reconstructor->getLatentDim() != 2){
            std::cerr << "Cannot produce a map of a latent space whose dim is different than 2" << std::endl;
            return;
        }
        bitmap_image map(width, height);
        const double color_max(255);
        double r, g, b;
        Scalar xstep(latentHeight/width), ystep(latentHeight/height);
        RowVector<Scalar> brdf(reconstructor->getBRDFCoeffNb());
        Vector<Scalar> coord(2);
        coord << (xstep-latentWidth)*0.5, (ystep-latentHeight)*0.5;
        std::cout << "Compute albedo map" << std::endl;
        for (unsigned int pixx(0); pixx < width; ++pixx)
        {
            coord[0] += xstep;
            for (unsigned int pixy(0); pixy < height; ++pixy)
            {
                progressBar(double(pixx*height+pixy) / (width*height));
                coord[1] += ystep;
                reconstructor->reconstruct(brdf, coord);
                // clamp BRDF values in [0; +inf)
                brdf = brdf.cwiseMax(Scalar(0));
                Albedo::computeAlbedo<Scalar>(brdf, r, g, b, albedoSampling);
                map.set_pixel(pixx, pixy,
                              r*color_max,
                              g*color_max,
                              b*color_max);
            }
            coord[1] = (ystep-latentHeight)*0.5;
        }
        std::cout << std::endl;
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

        file.precision(100);
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
