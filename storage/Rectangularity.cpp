//
// Created by Sharath Chandra on 17-Mar-23.
//

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>
//#include <fstream>
#include <iostream>
#include <algorithm>
#include <CGAL/Optimal_bounding_box/Oriented_bounding_box_traits_3.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_3                                             Point;
typedef CGAL::Surface_mesh<Point>                              Surface_mesh;

struct Vertex {
    int id;
    double x, y, z;
};

struct Vector {
    double x, y, z;
};

struct Face {
    int fid;
    std::vector<int> tri_vertices;
};

//int main(int argc, char** argv) {
//
//    // Reading the obj file
////    const char* inobj = (argc > 1) ? argv[1] : "../data/cube_out.obj";
////    std::ifstream input_stream;
////    input_stream.open(inobj);
//
//    // creating a surface mesh
////    Surface_mesh sm;
////    PMP::IO::read_polygon_mesh(inobj, sm);
////    std::cout << "surface mesh: " << sm << std::endl;
//
////    std::array<Point, 8> obb_points;
////    for (auto& pt: obb_points) {
////        pt = Point(0, 0, 0);
////    }
//
//////    std::vector<Point> obb_points;
////    CGAL::oriented_bounding_box(sm, obb_points, CGAL::parameters::use_convex_hull(true));
////
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
//    // code starts here
//
//    std::vector<int> checker;
//    std::vector<Point> objects_points_3d;
//    for (auto& [key, face]: faces) {
//
//        int tv0, tv1, tv2;  // triangle vertices ids
//        tv0 = face.tri_vertices[0]; tv1 = face.tri_vertices[1]; tv2= face.tri_vertices[2];
//
//        if (std::find(checker.begin(), checker.end(), tv0) == checker.end()) {
//            // vertex 0 doesn't exist in checker, so we can add it
//            checker.push_back(tv0);
//            objects_points_3d.emplace_back(Point(vertices[tv0].x, vertices[tv0].y, vertices[tv0].z));
//        }
//        if (std::find(checker.begin(), checker.end(), tv1) == checker.end()) {
//            // vertex 0 doesn't exist in checker, so we can add it
//            checker.push_back(tv1);
//            objects_points_3d.emplace_back(Point(vertices[tv1].x, vertices[tv1].y, vertices[tv1].z));
//        }
//        if (std::find(checker.begin(), checker.end(), tv2) == checker.end()) {
//            // vertex 0 doesn't exist in checker, so we can add it
//            checker.push_back(tv2);
//            objects_points_3d.emplace_back(Point(vertices[tv2].x, vertices[tv2].y, vertices[tv2].z));
//        }
//
//    }
//
////    CGAL::oriented_bounding_box<K> obb(objects_points_3d.begin(), objects_points_3d.end());
//    std::array<Point, 8> obb_points;
//    CGAL::oriented_bounding_box(objects_points_3d, obb_points, CGAL::parameters::use_convex_hull(true));
//
//
//    return 0;
//}