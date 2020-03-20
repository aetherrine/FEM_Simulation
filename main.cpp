#define TINYOBJLOADER_IMPLEMENTATION
#include "library/tiny_obj_loader.h"
#include "tetrahedral.h"
#include <iostream>
#include <vector>
#include <fstream>

void outputOBJ(std::vector<Particle*> particle_list, std::string original_file, int idx){
    std::ifstream in(original_file);
    std::ofstream out("output/cube_"+std::to_string(idx)+".obj");
    int i = 0;
    if (in.is_open()){
        std::string line;
        while (getline(in, line)){
            if (line.substr(0,2) != "v "){
                out << line << "\n";
            }
            else{
                out << "v";
                out << " " << particle_list[i]->x();
                out << " " << particle_list[i]->y();
                out << " " << particle_list[i]->z() << "\n";
                i++;
            }
        }
        in.close();
        out.close();
    }
}

void precomputation(std::vector<Tetrahedral*> meshes, std::vector<Matrix3f>& B, std::vector<float>& W){
    for (int i=0; i<meshes.size(); i++){
        Matrix3f D_m;
        D_m << meshes[i]->v[0]->x()-meshes[i]->v[3]->x(), meshes[i]->v[1]->x()-meshes[i]->v[3]->x(), meshes[i]->v[2]->x()-meshes[i]->v[3]->x(),
               meshes[i]->v[0]->y()-meshes[i]->v[3]->y(), meshes[i]->v[1]->y()-meshes[i]->v[3]->y(), meshes[i]->v[2]->y()-meshes[i]->v[3]->y(),
               meshes[i]->v[0]->z()-meshes[i]->v[3]->z(), meshes[i]->v[1]->z()-meshes[i]->v[3]->z(), meshes[i]->v[2]->z()-meshes[i]->v[3]->z();

        B.push_back(1.0/6.0 * D_m.inverse());
        W.push_back(1.0/6.0 * D_m.determinant());
    }
}

Matrix3f VK_material(Matrix3f deform_grad){
    Matrix3f I = Matrix3f::Identity(3,3);
    Matrix3f energy = 0.5 * (deform_grad.transpose()*deform_grad - I);
    Matrix3f P = deform_grad * (2.0*0.1785*energy + 0.7141*energy.trace()*I);
    return P;
}

Matrix3f VK_material_differential(Matrix3f deform_grad, Matrix3f delta_deform_grad){
    Matrix3f I = Matrix3f::Identity(3,3);

    Matrix3f energy = 0.5 * (deform_grad.transpose()*deform_grad - I);
    Matrix3f delta_energy = 0.5 * (delta_deform_grad.transpose()*deform_grad + deform_grad.transpose()*delta_deform_grad);
    Matrix3f delta_p = delta_deform_grad*(2.0*0.1785*energy+0.7141*energy.trace()*I) + deform_grad*(2.0*0.1785*delta_energy+0.7141*delta_energy.trace()*I);
    return delta_p;
}

Matrix3f PK_stress_tensor_corotated(Matrix3f deform_grad){
    // JacobiSVD<Matrix3f> svd(deform_grad);

    Matrix3f I = Matrix3f::Identity(3,3);
    Affine3f t;
    t = deform_grad;
    Matrix3f R = t.rotation(); // rotation matrix in polar decomposition
    Matrix3f P = 2.0*0.1785*(deform_grad-R) + 0.7141*(R.transpose()*deform_grad-I).trace()*R;
    return P;
}

Matrix3f Neohookean(Matrix3f deform_grad){
    Matrix3f P = 0.1785*(deform_grad-0.1785*deform_grad.transpose().inverse()) + 0.7141*log(deform_grad.determinant())*deform_grad.transpose().inverse();
    return P;
}

void ComputeElasticForces(std::vector<Tetrahedral*> new_meshes, std::vector<Matrix3f>& B, std::vector<float>& W){
    for (int i=0; i<new_meshes.size(); i++){
        Matrix3f D_s;
        D_s << new_meshes[i]->v[0]->x()-new_meshes[i]->v[3]->x(), new_meshes[i]->v[1]->x()-new_meshes[i]->v[3]->x(), new_meshes[i]->v[2]->x()-new_meshes[i]->v[3]->x(),
               new_meshes[i]->v[0]->y()-new_meshes[i]->v[3]->y(), new_meshes[i]->v[1]->y()-new_meshes[i]->v[3]->y(), new_meshes[i]->v[2]->y()-new_meshes[i]->v[3]->y(),
               new_meshes[i]->v[0]->z()-new_meshes[i]->v[3]->z(), new_meshes[i]->v[1]->z()-new_meshes[i]->v[3]->z(), new_meshes[i]->v[2]->z()-new_meshes[i]->v[3]->z();

        Matrix3f deform_grad = D_s * B[i];

        Matrix3f P = VK_material(deform_grad);
        // Matrix3f P = PK_stress_tensor_corotated(deform_grad);
        // Matrix3f P = Neohookean(deform_grad); 
        Matrix3f H = -W[i] * P * B[i].transpose();

        new_meshes[i]->v[0]->force += H.col(0);
        new_meshes[i]->v[1]->force += H.col(1);
        new_meshes[i]->v[2]->force += H.col(2);
        new_meshes[i]->v[3]->force += (-H.col(0) - H.col(1) - H.col(2));
    }
}

void ComputeForceDifferentials(std::vector<Tetrahedral*> new_meshes, std::vector<Matrix3f>& B, std::vector<float>& W){
    for (int i=0; i<new_meshes.size(); i++){
        Matrix3f D_s;
        D_s << new_meshes[i]->v[0]->x()-new_meshes[i]->v[3]->x(), new_meshes[i]->v[1]->x()-new_meshes[i]->v[3]->x(), new_meshes[i]->v[2]->x()-new_meshes[i]->v[3]->x(),
               new_meshes[i]->v[0]->y()-new_meshes[i]->v[3]->y(), new_meshes[i]->v[1]->y()-new_meshes[i]->v[3]->y(), new_meshes[i]->v[2]->y()-new_meshes[i]->v[3]->y(),
               new_meshes[i]->v[0]->z()-new_meshes[i]->v[3]->z(), new_meshes[i]->v[1]->z()-new_meshes[i]->v[3]->z(), new_meshes[i]->v[2]->z()-new_meshes[i]->v[3]->z();
        
        Matrix3f D_s_delta;
        D_s_delta << new_meshes[i]->v[0]->displacement()[0]-new_meshes[i]->v[3]->displacement()[0],
                     new_meshes[i]->v[1]->displacement()[0]-new_meshes[i]->v[3]->displacement()[0],
                     new_meshes[i]->v[2]->displacement()[0]-new_meshes[i]->v[3]->displacement()[0],
                     new_meshes[i]->v[0]->displacement()[1]-new_meshes[i]->v[3]->displacement()[1],
                     new_meshes[i]->v[1]->displacement()[1]-new_meshes[i]->v[3]->displacement()[1],
                     new_meshes[i]->v[2]->displacement()[1]-new_meshes[i]->v[3]->displacement()[1],
                     new_meshes[i]->v[0]->displacement()[2]-new_meshes[i]->v[3]->displacement()[2],
                     new_meshes[i]->v[1]->displacement()[2]-new_meshes[i]->v[3]->displacement()[2],
                     new_meshes[i]->v[2]->displacement()[2]-new_meshes[i]->v[3]->displacement()[2];

        Matrix3f deform_grad = D_s * B[i];
        Matrix3f delta_deform_grad = D_s_delta * B[i];
        Matrix3f delta_p = VK_material_differential(deform_grad, delta_deform_grad);
        Matrix3f delta_h = -W[i] * delta_p * B[i].transpose();
        new_meshes[i]->v[0]->force += delta_h.col(0);
        new_meshes[i]->v[1]->force += delta_h.col(1);
        new_meshes[i]->v[2]->force += delta_h.col(2);
        new_meshes[i]->v[3]->force += (-delta_h.col(0) - delta_h.col(1) - delta_h.col(2));
    }
}

void resetForce(std::vector<Particle*> particles){
    for (auto p:particles)
        p->force = Vector3f(0,0,0);
}

void exertForce(std::vector<Particle*> particles){
    Vector3f gravity(0, -0.98, 0);
    for (auto p:particles){
        p->force += gravity * p->mass;
    }
}

void collision(std::vector<Particle*> particles, float dt){
    for (auto p:particles){
        if (p->position[1] <= 0){
            p->velocity = -p->velocity;
        }
    }
}

void forwardEuler(std::vector<Particle*> particles, float dt){
    collision(particles, dt);
    for (auto p:particles){
        p->velocity += p->force / p->mass * dt;
        p->position += p->velocity * dt;
    }
}

void backwardEuler(std::vector<Particle*> particles, float dt){
    collision(particles, dt);
}

bool compareVertex(Particle* p1, Particle* p2) {
    if (p1->x() != p2->x())
        return p1->x() < p2->x();
    if (p1->y() != p2->y())
        return p1->y() < p2->y();
    return p1->z() < p2->z();
}

bool compareIdx(Particle* p1, Particle* p2) {
    return p1->index < p2->index;
}

int main(){
    std::string OBJ_PATH = "../models/new444.obj";
    std::vector<Tetrahedral*> tetrahedral_list;
    std::vector<Particle*> particle_list;

    std::ifstream file(OBJ_PATH);
    if (file.is_open()) {
        std::string line;
        int i = 0;
        while (getline(file, line)) {
            std::vector<std::string> result;
            std::istringstream iss(line);
            for(std::string s; iss >> s; )
                result.push_back(s);

            if (result.size() == 0)
                continue;
            if (result[0] == "vt")
                break;
            if (result[0] != "v")
                continue;

            Vector3f vertex(std::stof(result[1]), std::stof(result[2])+5.0, std::stof(result[3]));
            Particle* p = new Particle(vertex, i++);
            particle_list.push_back(p);
        }
        file.close();
    }

    // insert all 5*8 volume tetrahedral meshes.
    sort(particle_list.begin(), particle_list.end(), compareVertex);
    int offset_x = 0;
    int offset_y = 0;
    int offset_z = 0;
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
            for (int k=0; k<4; k++){
                Tetrahedral* tetMesh1 = new Tetrahedral(particle_list[0+offset_x+offset_y+offset_z],particle_list[5+offset_x+offset_y+offset_z],particle_list[6+offset_x+offset_y+offset_z],particle_list[30+offset_x+offset_y+offset_z]);
                Tetrahedral* tetMesh2 = new Tetrahedral(particle_list[0+offset_x+offset_y+offset_z],particle_list[25+offset_x+offset_y+offset_z],particle_list[26+offset_x+offset_y+offset_z],particle_list[30+offset_x+offset_y+offset_z]);
                Tetrahedral* tetMesh3 = new Tetrahedral(particle_list[6+offset_x+offset_y+offset_z],particle_list[31+offset_x+offset_y+offset_z],particle_list[26+offset_x+offset_y+offset_z],particle_list[30+offset_x+offset_y+offset_z]);
                Tetrahedral* tetMesh4 = new Tetrahedral(particle_list[0+offset_x+offset_y+offset_z],particle_list[1+offset_x+offset_y+offset_z],particle_list[26+offset_x+offset_y+offset_z],particle_list[6+offset_x+offset_y+offset_z]);
                Tetrahedral* tetMesh5 = new Tetrahedral(particle_list[0+offset_x+offset_y+offset_z],particle_list[26+offset_x+offset_y+offset_z],particle_list[6+offset_x+offset_y+offset_z],particle_list[30+offset_x+offset_y+offset_z]);

                tetrahedral_list.push_back(tetMesh1);
                tetrahedral_list.push_back(tetMesh2);
                tetrahedral_list.push_back(tetMesh3);
                tetrahedral_list.push_back(tetMesh4);
                tetrahedral_list.push_back(tetMesh5);
                offset_x += 25;
            }
            offset_x = 0;
            offset_z += 5;
        }
        offset_x = offset_z = 0;
        offset_y += 1;
    }
    sort(particle_list.begin(), particle_list.end(), compareIdx);

    // Euler integration
    // deformed shape(Ds) = deformation gardient(F) * reference shape(Dm)
    std::vector<float> undeformed_vol;
    std::vector<Matrix3f> B_m;
    precomputation(tetrahedral_list, B_m, undeformed_vol);

    // for (int t=18; t<24; t++){
    //     particle_list[t]->force[1] -= 1;
    // }
    // particle_list[24]->force[1] -= 1;
    // particle_list[25]->force[1] -= 1;
    // particle_list[26]->force[1] -= 1;
    // for (auto vertex:particle_list)
    //     vertex->position += Vector3f(0,-1,0);

    float delta_t = 0.01;
    for (int i=0; i<100; i++){
        for (int j=0; j<5; j++){
            resetForce(particle_list);
            exertForce(particle_list);
            ComputeElasticForces(tetrahedral_list, B_m, undeformed_vol);
            forwardEuler(particle_list, delta_t);
            // backwardEuler(particle_list, delta_t);
        }
        outputOBJ(particle_list, OBJ_PATH, i);
    }





//---------------------------- TESTING ----------------------------
    // Test mesh loading
    // for (int i=0; i<4; i++){
    //     std::cout<< tetrahedral_list[20].v[i].x() << ", "
    //              << tetrahedral_list[20].v[i].y() << ", "
    //              << tetrahedral_list[20].v[i].z() << ", " <<std::endl;
    // }

    // Test output for Houdini
    // sort to reproduce .OBJ files. 
    // sort(particle_list.begin(), particle_list.end(), compareIdx);
    // for (int i=0; i<20; i++){
    //     for (auto& item : particle_list){
    //         item.position[1] -= 0.1;
    //     }
    //     outputOBJ(particle_list, OBJ_PATH, i);
    // }

    // Test polar decomposition
    // Matrix2f tmp;
    // tmp << 1.300, -0.375,
    //        0.750, 0.650;
    // Affine2f t;
    // t = tmp;
    // Matrix2f rot = t.rotation();
    // std::cout<<tmp<<std::endl;
    // std::cout<<"R: "<<rot<<std::endl;
    // std::cout<<"S: "<<rot.inverse()*tmp<<std::endl;
    // return 0;
//---------------------------- TESTING ----------------------------



//---------------------------- TESTING FEM functions ----------------------------

    // middle up puuuuuuuuuuull
    // sort(particle_list.begin(), particle_list.end(), compareIdx);
    // particle_list[21].position[1] += 1;
    // outputOBJ(particle_list, OBJ_PATH, 0);

    // ComputeForceDifferentials(tetrahedral_list, B_m, undeformed_vol);
    // std::cout<<particle_list[21].force<< std::endl;    

    
        // for (auto& item : particle_list){
        //     std::cout<< item.force;    
        //     // item.position[1] -= 0.1;
        // }
        // outputOBJ(particle_list, OBJ_PATH, i);
//---------------------------- TESTING END ----------------------------

    // hardcoded volume tetrahedral mesh
    // Vector3f origin(0,0,0);

    // for (int i=0; i<2; i++){
    //     for (int j=0; j<2; j++){
    //         for (int k=0; k<2; k++){
    //             Tetrahedral* tetMesh1 = new Tetrahedral(Vector3f(origin[0], origin[1], origin[2]),
    //                                                     Vector3f(origin[0], origin[1]+1, origin[2]),
    //                                                     Vector3f(origin[0], origin[1]+1, origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1]+1, origin[2]));

    //             Tetrahedral* tetMesh2 = new Tetrahedral(Vector3f(origin[0], origin[1], origin[2]),
    //                                                     Vector3f(origin[0]+1, origin[1], origin[2]),
    //                                                     Vector3f(origin[0]+1, origin[1], origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1]+1, origin[2]));

    //             Tetrahedral* tetMesh3 = new Tetrahedral(Vector3f(origin[0]+1, origin[1], origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1]+1, origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1]+1, origin[2]),
    //                                                     Vector3f(origin[0], origin[1]+1, origin[2]+1));

    //             Tetrahedral* tetMesh4 = new Tetrahedral(Vector3f(origin[0], origin[1], origin[2]),
    //                                                     Vector3f(origin[0], origin[1], origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1], origin[2]+1),
    //                                                     Vector3f(origin[0], origin[1]+1, origin[2]+1));

    //             Tetrahedral* tetMesh5 = new Tetrahedral(Vector3f(origin[0], origin[1], origin[2]),
    //                                                     Vector3f(origin[0]+1, origin[1], origin[2]+1),
    //                                                     Vector3f(origin[0], origin[1]+1, origin[2]+1),
    //                                                     Vector3f(origin[0]+1, origin[1]+1, origin[2]));

    //             tetrahedral_list.push_back(tetMesh1);
    //             tetrahedral_list.push_back(tetMesh2);
    //             tetrahedral_list.push_back(tetMesh3);
    //             tetrahedral_list.push_back(tetMesh4);
    //             tetrahedral_list.push_back(tetMesh5);
    //             origin[0] += 1;
    //         }
    //         origin[0] = 0;
    //         origin[1] += 1;
    //     }
    //     origin[0] = origin[1] = 0;
    //     origin[2] += 1;
    // }



    /*tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;
    std::vector<Tetrahedral*> tetrahedral_list;

    bool load_ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, OBJ_PATH.c_str());
    if (!load_ret)
        throw std::runtime_error(warn + err);

    for (size_t s = 0; s < shapes.size(); s++) { // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            int fv = shapes[s].mesh.num_face_vertices[f];

            for (size_t v = 0; v < fv; v++) { // Loop over vertices in the face.
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
                tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
                tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
                tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
                tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
                tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
                tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
                tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
                // Optional: vertex colors
                // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
                // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
                // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
            }
            index_offset += fv;
        }
    }*/
}
