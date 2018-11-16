#include <iostream>
#include <Eigen/Dense>
#include <GL/gl.h>
#include <GL/glu.h>
#include "BRDFReader/brdfreader.h"

using Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    /*MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;*/
    BRDFReader reader;
    double* brdf;
    reader.read_brdf("../datas/silver-metallic-paint2.binary", brdf);


    return 0;
}


