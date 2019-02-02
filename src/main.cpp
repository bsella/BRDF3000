#include <chrono>

#include "Parametrisation/types.h"
#include "Parametrisation/Parametrisation.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"
#include "Optimisation/OptiDataWriter.h"

using namespace ChefDevr;
using Scalar = double;

int main(int numArguments, const char *argv[]) {
    // Récupérer les arguments
    // Afficher le message d'aide si faut
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double, std::micro> duration;
        
    BRDFReader reader;
    const unsigned int dim = 2;
    const Scalar minStep = 0.5;
    const char *brdfsDir = "../data/";
    const std::string mapPath("../map.bmp");
    const unsigned int mapWidth(8), mapHeight(8), albedoSampling(4);

    auto Z = reader.createZ<Scalar>(brdfsDir);
    std::cout << Z.row(0).transpose() << std::endl;
    RowVector<Scalar> meanBRDF(Z.cols());
    centerMat(Z, meanBRDF);
    
    
    OptimisationSolver<Scalar> optimisation{minStep, Z, dim};
    
    start = std::chrono::system_clock::now();
    optimisation.optimizeMapping();
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Optimisation took " << duration.count() << " micro seconds" << std::endl;
    
    BRDFReconstructor<Scalar> reconstructor(
                Z,
                optimisation.getInverseMapping(),
                optimisation.getLatentVariables(),
                meanBRDF.transpose(),
                dim);
    
    std::cout << "Reconstruction error for BRDF n°0 : " << reconstructor.reconstructionError(0) << std::endl;
    
    /*
    start = std::chrono::system_clock::now();
    writeAlbedoMap<Scalar>(mapPath,
                   Z,
                   optimisation.getLatentVariables(),
                   optimisation.getInverseMapping(),
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
