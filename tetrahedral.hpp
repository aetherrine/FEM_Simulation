#include "library/eigen/Eigen/Eigen"
using namespace Eigen;

class Tetrahedral{

public:
    Vector3f v[4]; // vertices

    Tetrahedral(Vector3f v0, Vector3f v1, Vector3f v2, Vector3f v3){
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
    }

    Eigen::Vector3f x() const { return v[0]; }
    Eigen::Vector3f y() const { return v[1]; }
    Eigen::Vector3f z() const { return v[2]; }

};