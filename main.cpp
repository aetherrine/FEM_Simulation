#define TINYOBJLOADER_IMPLEMENTATION
#include "library/tiny_obj_loader.h"
#include "tetrahedral.h"
#include <iostream>
#include <vector>
#include <fstream>

// void precomputation(std::vector<Tetrahedral*> meshes, std::vector<Matrix3f>& B, std::vector<float>& W){
//     for (int i=0; i<meshes.size(); i++){
//         Matrix3f D;
//         D << meshes[i][0].x()-meshes[i][3].x(), meshes[i][1].x()-meshes[i][3].x(), meshes[i][2].x()-meshes[i][3].x(),
//              meshes[i][0].y()-meshes[i][3].y(), meshes[i][1].y()-meshes[i][3].y(), meshes[i][2].y()-meshes[i][3].y(),
//              meshes[i][0].z()-meshes[i][3].z(), meshes[i][1].z()-meshes[i][3].z(), meshes[i][2].z()-meshes[i][3].z();
//         B[i] = D.inverse();
//         W[i] = 1/6 * D.determinant();
//     }
// }

// // TODO || !TODO
// bool loader(std::string Path){
//     if (Path.substr(Path.size() - 4, 4) != ".obj")
//         return false;

//     std::ifstream file(Path);

//     if (!file.is_open())
//         return false;

//     std::string curline;
//     while (std::getline(file, curline)){

//     }
// }

int main(){
    std::string OBJ_PATH = "../models/yet_another_cube.obj";
    std::vector<Tetrahedral> tetrahedral_list;
    std::vector<Particle> particle_list;

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

            Vector3f vertex(std::stof(result[1]), std::stof(result[2]), std::stof(result[3]));
            Particle p(vertex, i);
            particle_list.push_back(p);
        }
        file.close();
    }
    // for (auto i : particle_list)
    //     std::cout<<i.position[0]<<", "<<i.position[1]<<", "<<i.position[2]<<std::endl;

    tetrahedral_list.push_back(Tetrahedral(particle_list[0], particle_list[0], particle_list[0], particle_list[0]));

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