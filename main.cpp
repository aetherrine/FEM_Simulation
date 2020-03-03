#define TINYOBJLOADER_IMPLEMENTATION
#include "library/tiny_obj_loader.h"
#include <iostream>

int main(){
    std::string OBJ_PATH = "../models/cube/cube.obj";
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

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
                std::cout<<vx<<", "<<vy<<", "<<vz<<std::endl;
            }
            index_offset += fv;
        }
    }
}