#include <chrono>

#include "Parametrisation/types.h"
#include "BRDFReader/BRDFReader.h"
#include "Optimisation/OptimisationSolver.h"
#include "Parametrisation/Parametrisation.h"
using namespace ChefDevr;
using Scalar = double;

void test(
    Scalar* a,
    const Vector<Scalar>& X,
    const Vector<Scalar>& coordRef,
    const unsigned int dim,
    const unsigned int nb_data)
{
    std::cout << "a" << a << std::endl;
    std::cout << "X" << X << std::endl;
    std::cout << "coref" << coordRef << std::endl;
    std::cout << "dim" << dim << std::endl;
    std::cout << "njbdtata" << nb_data << std::endl;
}

void test2(const Vector<Scalar>& a)
{
    std::cout << "const vector access ok" << std::endl;
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
    Vector<Scalar> b(Z.col(1));
    test(Z.col(0).data(),b,b.segment(4,6),1,1);
	//OptimisationSolver<Scalar> optimisation{minStep, Z, dim};
    
    start = std::chrono::system_clock::now();
    //optimisation.optimizeMapping();
    end = std::chrono::system_clock::now();
    
    int elapsed_ms(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());
    
    std::cout << "Optimisation took " << elapsed_ms << " ms" << std::endl;
	exit(0);
}
