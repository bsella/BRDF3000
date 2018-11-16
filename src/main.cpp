#include <iostream>
#include <Eigen/Dense>

#include <Eigen/Core>
#include <GL/gl.h>
#include <GL/glu.h>
#include "BRDFReader/brdfreader.h"
#include <vector>
#include <string>
#include <iostream>

#include <stdio.h>

using namespace Eigen;

int main(int, char *[]) {
	std::vector<std::string> files={
		"../data/silver-metallic-paint2.binary",
		"../data/cherry-235.binary",
		"../data/chrome.binary",
		"../data/natural-209.binary",
		"../data/pink-fabric.binary"
	};

	const uint BRDFsize = 1458000*3; //16*16*64*64*3;

	MatrixXd merleData(BRDFsize,files.size());

	BRDFReader reader;
	for(uint i=0; i<files.size();i++){
		double *brdf = (double*)malloc(BRDFsize*sizeof(double));
		reader.read_brdf(files[i].c_str(), brdf);
		Map<VectorXd> vec(brdf,BRDFsize);
		merleData.col(i)=vec- Map<VectorXd>::Ones(BRDFsize)*vec.mean();
		free(brdf);
	}
//	std::cout << merleData<< std::endl;
	std::cout << merleData.transpose()*merleData << std::endl;
	return 0;
}


