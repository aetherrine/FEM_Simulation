#include "library/eigen/Eigen/Eigen"
using namespace Eigen;

class Particle{

public:
    Vector3f position;
    Vector3f last_position;
    Vector3f force;
    Vector3f velocity;
    Vector3f acceleration;

    int index;

    Particle();
    Particle(Vector3f pos, int idx){
        position = last_position = pos;
        index = idx;
        force = velocity = acceleration = Vector3f(0,0,0);
    }
};