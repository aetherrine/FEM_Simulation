#include "library/eigen/Eigen/Eigen"
#include "particle.h"
using namespace Eigen;

class Triangle{

public:
    std::string idx[3]; // vertices

    Triangle(std::string v0, std::string v1, std::string v2){
        idx[0] = v0;
        idx[1] = v1;
        idx[2] = v2;
    }

};
