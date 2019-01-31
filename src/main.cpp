#include <chrono>

#include "Parametrisation/types.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"

using namespace ChefDevr;
using Scalar = double;

void test(const Vector<Scalar>& v)
{
    std::cout << v << std::endl;
}

int main(int numArguments, const char *argv[]) {
	// Récupérer les arguments
	// Afficher le message d'aide si faut
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

	BRDFReader reader;
	const unsigned int dim = 2;
	const Scalar minStep = 0.5;
	const char *fileDirectory = "../data/";

	auto Z = reader.createZ<Scalar>(fileDirectory);
    test(Z.col(1).segment(1,5));
    
	//OptimisationSolver<Scalar> optimisation{minStep, Z, dim};
    
    start = std::chrono::system_clock::now();
    //optimisation.optimizeMapping();
    end = std::chrono::system_clock::now();
    
    int elapsed_ms(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());
    
    std::cout << "Optimisation took " << elapsed_ms << " ms" << std::endl;
	exit(0);
}
