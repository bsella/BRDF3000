#include <chrono>
#include <cstdio>

#include "Parametrisation/types.h"
#include "Parametrisation/Parametrisation.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"
#include "Optimisation/OptiDataWriter.h"
#include "Optimisation/Albedo.h"

using namespace ChefDevr;
using Scalar = long double;

template <typename Scalar>
void writeBRDF(const std::string& path, const RowVector<Scalar>& brdf)
{
    std::ofstream file(path, ios::binary);
    
    if (!file.is_open()){
        std::cerr << "Could not create file \"" << path<< "\"" << std::endl;
    }
    
    int dims[3]{
        BRDFReader::samplingResolution_thetaH,
        BRDFReader::samplingResolution_thetaD,
        BRDFReader::samplingResolution_phiD/2};
        
    file.write(reinterpret_cast<char*>(&dims[0]), sizeof(int));
    file.write(reinterpret_cast<char*>(&dims[1]), sizeof(int));
    file.write(reinterpret_cast<char*>(&dims[2]), sizeof(int));
    
    double conv;
    for (unsigned int i(0); i < brdf.cols(); ++i)
    {    
        conv = brdf[i];
        file.write(reinterpret_cast<char*>(&conv), sizeof(double));
    }
}

int main(int numArguments, const char *argv[]) {
    std::cout << std::endl;
    // Récupérer les arguments
    // Afficher le message d'aide si faut
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double, std::milli> duration;
    BRDFReader reader;
    const unsigned int dim = 2;
    const Scalar minStep = 0.005;
    const char *brdfsDir = "../data/";
    const std::string mapPath("../map.bmp"), optiDataPath("../paramtrzData");
    const unsigned int mapWidth(32), mapHeight(32), albedoSampling(16);
    const unsigned int reconstBRDFindex(0);
    double r, g, b;

    auto Z = reader.createZ<Scalar>(brdfsDir);
    RowVector<Scalar> meanBRDF(Z.cols());
    centerMat(Z, meanBRDF);
    RowVector<Scalar> brdf_r(Z.cols());
    
    
    OptimisationSolver<Scalar> optimizer{minStep, Z, dim};
    
    std::cout << "==== optimize mapping ===" << std::endl;
    start = std::chrono::system_clock::now();
    optimizer.optimizeMapping();
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Optimisation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;
    
    BRDFReconstructor<Scalar> reconstructor(
                Z,
                optimizer.getInverseMapping(),
                optimizer.getLatentVariables(),
                meanBRDF,
                dim);
    
    start = std::chrono::system_clock::now();
    reconstructor.reconstruct(brdf_r, optimizer.getLatentVariables().segment(reconstBRDFindex*dim,dim));
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout <<"Reconstruction of " << reader.getBRDFFilenames()[reconstBRDFindex] << " took " << duration.count()*0.001 << " seconds" << std::endl << std::endl;
    writeBRDF<Scalar>("../r_" + reader.getBRDFFilenames()[reconstBRDFindex], brdf_r);
    
    for (unsigned int i(0); i < std::min(static_cast<int>(Z.rows()), 5); ++i)
    {
        std::cout << "Reconstruction error for " << reader.getBRDFFilenames()[i] <<  " : " << reconstructor.reconstructionError(0) << std::endl;
    }
    std::cout << std::endl;
        
    start = std::chrono::system_clock::now();
    Albedo::computeAlbedo<Scalar>(brdf_r, r, g, b, albedoSampling);
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Albedo computing took " << duration.count()*0.001 << " seconds" << std::endl << std::endl;

    start = std::chrono::system_clock::now();
    writeAlbedoMap<Scalar>(
        mapPath,
        reconstructor,
        albedoSampling,
        mapWidth,
        mapHeight);
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Map computing took " << duration.count()*0.001 << " seconds" << std::endl;

    writeParametrisationData<Scalar>(
        optiDataPath,
        reader.getBRDFFilenames(),
        optimizer.getLatentVariables(),
        optimizer.getInverseMapping(),
        dim);
    
    exit(EXIT_SUCCESS);
}
