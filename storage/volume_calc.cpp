#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "definitions.h"
#include "geomtools.h"

#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

//struct Vertex {
//    int id;
//    double x, y, z;
//};
//
//struct Vector {
//    double x, y, z;
//};
//
//struct Face {
//    int fid;
//    std::vector<int> tri_vertices;
//};
//
//Vector Vertex_to_Vector (Vertex& vertex) {
//    Vector result;
//    result.x = vertex.x; result.y = vertex.y; result.z = vertex.z;
//    return result;
//}
//
//Vector Vector_Diff(Vector& vector1, Vector& vector2) {
//    Vector result;
//    result.x = vector1.x - vector2.x;
//    result.y = vector1.y - vector2.y;
//    result.z = vector1.z - vector2.z;
//    return result;
//}
//
//Vector Cross(Vector& vector1, Vector& vector2) {
//    Vector result;
//    result.x = vector1.y * vector2.z - vector1.z * vector2.y;
//    result.y = vector1.z * vector2.x - vector1.x * vector2.z;
//    result.z = vector1.x * vector2.y - vector1.y * vector2.x;
//    return result;
//}
//
//double Dot(Vector& vector1, Vector& vector2) {
//    double result_x, result_y, result_z, result;
//    result_x = vector1.x * vector2.x;
//    result_y = vector1.y * vector2.y;
//    result_z = vector1.z * vector2.z;
//    result = result_x + result_y + result_z;
//
//    return result;
//}
//
//double Volume_Calc(Face& face, std::map<int, Vertex> vertices) {
//    Vector vector_a, vector_b, vector_c;
//    Vector vector_ad, vector_bd, vector_cd;
//
//    vector_a = Vertex_to_Vector(vertices[face.tri_vertices[0]]);
//    vector_b = Vertex_to_Vector(vertices[face.tri_vertices[1]]);
//    vector_c = Vertex_to_Vector(vertices[face.tri_vertices[2]]);
//
//    Vector vector_out;
//    vector_out.x = 1000; vector_out.y = 1000; vector_out.z = 1000;
//
//    vector_ad = Vector_Diff(vector_a, vector_out);
//    vector_bd = Vector_Diff(vector_b, vector_out);
//    vector_cd = Vector_Diff(vector_c, vector_out);
//
//    Vector cross_bd_cd;
//    cross_bd_cd = Cross(vector_bd, vector_cd);
//
//    double volume;
//    Vector dot_a_bc;
//    volume = (1.0/6.0) * Dot(vector_ad, cross_bd_cd);
//
//    return volume;
//}
//

//int main(int argc, const char * argv[]) {
//
//    std::map<int, Vertex> vertices;
//    std::map<int, Face> faces;
//
//    // Reading the obj file
//    const char* inobj = (argc > 1) ? argv[1] : "../data/torus_out.obj";
//    std::ifstream input_stream;
//    input_stream.open(inobj);
//
//    std::string line;
//
//    //creating a vertex outside
//    Vector vector_out;
//    vector_out.x = 1000; vector_out.y = 1000; vector_out.z = 1000;
//
//
//    // Reading lines from obj
//    int i = 1;
//    while (std::getline(input_stream, line)) {
//        std::stringstream ss(line);
//        std::string prefix;
//        ss >> prefix;
//
//        if (prefix == "v") {  // reading the lines only with v
//            float x, y, z;
//            ss >> x >> y >> z;
////            std::cout << "vertex i: " << i << " x: " << x << " y: " << y << " z: " << z << std::endl;
//
//            vertices[i].id = i;
//            vertices[i].x = x;
//            vertices[i].y = y;
//            vertices[i].z = z;
//            i++;
//        }
//        // break if a new object is identified
//        if (prefix == "o")
//            break;
//    }
//    // checking if while is continuing to read the remaining lines
//    int faceid = 0;
//    while(std::getline(input_stream, line)) {
//        std::stringstream ss(line);
//        std::string prefix;
//        ss >> prefix;
//
//        if (prefix == "f") {
//            int v0, v1, v2;
//            ss >> v0 >> v1 >> v2;
//            faces[faceid].fid = faceid;
//            faces[faceid].tri_vertices.emplace_back(v0); faces[faceid].tri_vertices.emplace_back(v1); faces[faceid].tri_vertices.emplace_back(v2);
//        }
//        faceid++;
//
//        // break if a new object is identified
//        if (prefix == "o")
//            break;
//    }
//
//    // checking if vertices are loded or not
//    for (auto& [key, vertex] : vertices) {
//        std::cout << "vertex id: " << key << " - x: " << vertex.x << " y: " << vertex.y << " z: " << vertex.z << std::endl;
//    }
//    std::cout << "=============================" << std::endl;
//    // checking if all the faces are loaded or not
//    // calculating volume made by every triangle with outside vertex
//    double volume = 0;
//    for (auto& [key, f] : faces) {
//        std::cout << "face id: " << f.fid << " - v0: " << f.tri_vertices[0] << ", v1: " << f.tri_vertices[1] << ", v2: " << f.tri_vertices[2] << std::endl;
//
//        double face_vol;
//        face_vol = Volume_Calc(f, vertices);
//        volume += face_vol;
//        std::cout << "face id: " << f.fid << " - volume: " << face_vol << std::endl;
//        std::cout << "=============================" << std::endl;
//    }
//
//    std::cout << "Volume of Cube: " << volume;
//
//
//
//
//    return 0;
//}