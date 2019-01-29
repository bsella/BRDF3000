#include <stdlib.h>
// #include <experimental/filesystem>

#include "Parametrisation/types.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"

using namespace ChefDevr;
using Scalar = double;


int main(int numArguments, const char *argv[]) {
	// Récupérer les arguments
	// Afficher le message d'aide si faut

    
    // test filesystem
    //auto path = std::experimental::filesystem::current_path();
    
    /*
	BRDFReader<Scalar> reader;
	const unsigned int dim = 2;
	const Scalar minStep = 0.5;
	const char *fileDirectory = nullptr;

	auto Z = reader.create_brdfSet(fileDirectory);
	OptimisationSolver<Scalar> optimisation{minStep, Z, dim};
    */
    
	exit(0);
}
