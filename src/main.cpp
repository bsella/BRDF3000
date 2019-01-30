#include <stdlib.h>

#include "Parametrisation/types.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"

using namespace ChefDevr;
using Scalar = double;


int main(int numArguments, const char *argv[]) {
	// Récupérer les arguments
	// Afficher le message d'aide si faut


	BRDFReader reader;
	//const unsigned int dim = 2;
	//const Scalar minStep = 0.5;
	const char *fileDirectory = "/home/slowlys/Documents/Cours/Master 2/ChefDoeuvre/BRDF/Code/BRDF3000/data/";

	auto Z = reader.createZ<Scalar>(fileDirectory);

	//OptimisationSolver<Scalar> optimisation{minStep, Z, dim};

    
	exit(0);
}
