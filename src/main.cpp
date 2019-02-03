#include <chrono>

#include "Parametrisation/types.h"
#include "Parametrisation/Parametrisation.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"
#include "Optimisation/OptiDataWriter.h"

using namespace ChefDevr;
using Scalar = double;

void writeBRDF(const std::string& path, const Vector<Scalar>& brdf)
{
    std::ofstream file(path);
    if (!file.is_open()){
        std::cerr << "Could not create file \"" << path<< "\"" << std::endl;
    }
    file << brdf;
}

int main(int numArguments, const char *argv[]) {
    std::cout << std::endl;
    // Récupérer les arguments
    // Afficher le message d'aide si faut
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double, std::micro> duration;
    BRDFReader reader;
    const unsigned int dim = 2;
    const Scalar minStep = 0.01;
    const char *brdfsDir = "../data/";
    const std::string mapPath("../map.bmp");
    const unsigned int mapWidth(8), mapHeight(8), albedoSampling(4);

    auto Z = reader.createZ<Scalar>(brdfsDir);
    RowVector<Scalar> meanBRDF(Z.cols());
    centerMat(Z, meanBRDF);
    Vector<Scalar> brdf0_r(Z.cols());
    
    
    OptimisationSolver<Scalar> optimizer{minStep, Z, dim};
    
    start = std::chrono::system_clock::now();
    optimizer.optimizeMapping();
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Optimisation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;
    
    BRDFReconstructor<Scalar> reconstructor(
                Z,
                optimizer.getInverseMapping(),
                optimizer.getLatentVariables(),
                meanBRDF.transpose(),
                dim);
    
    start = std::chrono::system_clock::now();
    reconstructor.reconstruct(brdf0_r, optimizer.getLatentVariables().segment(0,dim));
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Reconstruction took " << duration.count()*0.000001 << " seconds" << std::endl << std::endl;
    writeBRDF("../brdf0.binary", brdf0_r);
    
    
    std::cout << "Reconstruction error for BRDF n°0 : " << reconstructor.reconstructionError(0) << std::endl << std::endl;
    
    /*
    start = std::chrono::system_clock::now();
    writeAlbedoMap<Scalar>(mapPath,
                   Z,
                   optimizer.getLatentVariables(),
                   optimizer.getInverseMapping(),
                   meanBRDF,
                   dim,
                   albedoSampling,
                   mapWidth,
                   mapHeight);
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Map computing took " << duration.count()*0.000001 << " seconds" << std::endl;
    */
    
    exit(0);
}
