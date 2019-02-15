#include <chrono>
#include <cstdio>

#include "Parametrisation/types.h"
#include "Parametrisation/ParametrisationWithZ.h"
#include "Parametrisation/ParametrisationSmallStorage.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"
#include "Optimisation/OptiDataWriter.h"
#include "Optimisation/Albedo.h"


#define WRONG_USAGE 1


using namespace ChefDevr;
//using Scalar = boost::multiprecision::float128;
using Scalar = long double;

template <typename Scalar>
void writeBRDF(const std::string& path, const RowVector<Scalar>& brdf)
{
    std::ofstream file(path, ios::binary);
    
    if (!file.is_open()){
        std::cerr << "Could not create file \"" << path<< "\"" << std::endl;
    }
    
    int dims[3]{
        MERLReader::samplingResolution_thetaH,
        MERLReader::samplingResolution_thetaD,
        MERLReader::samplingResolution_phiD/2};
        
    file.write(reinterpret_cast<char*>(&dims[0]), sizeof(int));
    file.write(reinterpret_cast<char*>(&dims[1]), sizeof(int));
    file.write(reinterpret_cast<char*>(&dims[2]), sizeof(int));
    
    double conv;
    for (unsigned int i(0); i < brdf.cols(); ++i)
    {    
        conv = (double)brdf[i];
        file.write(reinterpret_cast<char*>(&conv), sizeof(double));
    }
}


static void show_usage(const char *name_program)
{
    std::cerr << "Usage: " << name_program << " <option> " << std::endl
              << "Options:\n"
              << "\t-h,--help\t\tShow the usage\n"
              << "\t--smallRam\t\tThe program will keep the Ram storage low but will take longer to execute" << std::endl;
}


int main(int numArguments, const char *argv[]) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double, std::milli> duration{};
    BRDFReader reader;
    BRDFReconstructor<Scalar> *reconstructor;
    OptimisationSolver<Scalar> *optimizer;
    const unsigned int dim = 2;
    const Scalar minStep = 0.0005;
    const char *brdfsDir = "../data";
    const std::string mapPath("../map.bmp"), optiDataPath("../paramtrzData");
    const unsigned int mapWidth(200), mapHeight(200), albedoSampling(16);
    const unsigned int reconstBRDFindex(0);
    const double latentSize(8.);
    RowVector<Scalar> meanBRDF;
    ChefDevr::Matrix<Scalar> Z;
    double r, g, b;
    long num_brdf;
    bool smallStorage = false;


    if (numArguments > 2) {
        std::cerr << "Too much arguments" << std::endl;
        show_usage(argv[0]);
        exit(WRONG_USAGE);
    } else if (numArguments == 2) {
        std::string argument = argv[1];

        if (argument == "-h" || argument == "--help") {
            show_usage(argv[0]);
            exit(WRONG_USAGE);
        } else if (argument == "--smallRam") {
            smallStorage = true;
        } else {
            std::cerr << argument << " is not a valid argument" << std::endl;
            show_usage(argv[0]);
            exit(WRONG_USAGE);
        }
    }


    if (smallStorage) {
        start = std::chrono::system_clock::now();
        auto ZZt = reader.createZZt_centered<Scalar>(brdfsDir, meanBRDF);
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Loading ZZt took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;

        num_brdf = ZZt.rows();

        optimizer = new OptimisationSolver<Scalar>(meanBRDF.cols(), minStep, ZZt, dim);
        start = std::chrono::system_clock::now();
        optimizer->optimizeMapping();
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Optimisation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;

        start = std::chrono::system_clock::now();
        reconstructor = new BRDFReconstructorSmallStorage<Scalar>(optimizer->getInverseMapping(),
                                                                  optimizer->getLatentVariables(), meanBRDF, dim, reader.getBRDFFilePaths());
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Reconstructor creation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;
    } else {
        start = std::chrono::system_clock::now();
        Z = reader.createZ<Scalar>(brdfsDir);
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Loading Z took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;

        centerMat(Z, meanBRDF);
        num_brdf = Z.rows();

        const auto ZZt = Z * Z.transpose();
        optimizer = new OptimisationSolver<Scalar>(meanBRDF.cols(), minStep, ZZt, dim);
        start = std::chrono::system_clock::now();
        optimizer->optimizeMapping();
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Optimisation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;

        start = std::chrono::system_clock::now();
        reconstructor = new BRDFReconstructorWithZ<Scalar>(Z, optimizer->getInverseMapping(),
                                                           optimizer->getLatentVariables(), meanBRDF, dim);
        end = std::chrono::system_clock::now();
        duration = end - start;
        std::cout << "Reconstructor creation took " << duration.count() * 0.001<< " seconds" << std::endl << std::endl;
    }

    writeParametrisationData<Scalar>(
        optiDataPath,
        reader.getBRDFFilenames(),
        optimizer->getLatentVariables(),
        optimizer->getInverseMapping(),
        dim);
    
    RowVector<Scalar> brdf_r(reconstructor->getBRDFCoeffNb());
    
    start = std::chrono::system_clock::now();
    reconstructor->reconstruct(brdf_r, optimizer->getLatentVariables().segment(reconstBRDFindex*dim,dim));
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout <<"Reconstruction of " << reader.getBRDFFilenames()[reconstBRDFindex] << " took " << duration.count()*0.001 << " seconds" << std::endl << std::endl;
    
    writeBRDF<Scalar>("../r_" + reader.getBRDFFilenames()[reconstBRDFindex], brdf_r);
    
    for (unsigned int i(0); i < std::min(static_cast<int>(num_brdf), 5); ++i)
    {
        std::cout << "Reconstruction error for " << reader.getBRDFFilenames()[i] <<  " : " << reconstructor->reconstructionError(i) << std::endl;
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
        mapHeight,
        latentSize,
        latentSize);
    
    end = std::chrono::system_clock::now();
    duration = end - start;
    std::cout << "Map computing took " << duration.count()*0.001 << " seconds" << std::endl;


    delete reconstructor;
    delete optimizer;
    exit(EXIT_SUCCESS);
}
