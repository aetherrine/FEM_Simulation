#include "library/eigen/Eigen/Eigen"
#include "particle.h"
using namespace Eigen;

class Tetrahedral{

public:
    Particle v[4]; // vertices

    Tetrahedral(Particle v0, Particle v1, Particle v2, Particle v3){
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
    }

};