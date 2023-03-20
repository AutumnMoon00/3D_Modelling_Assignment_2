//
// Created by Sharath Chandra on 17-Mar-23.
//

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <CGAL/Random.h>

#include "definitions.h"
#include "geomtools.h"

#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Optimal_bounding_box/Oriented_bounding_box_traits_3.h>

#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

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
    std::vector<Point3> v3_coors;
};

Vector Vertex_to_Vector (Vertex& vertex) {
    Vector result;
    result.x = vertex.x; result.y = vertex.y; result.z = vertex.z;
    return result;
}

Vector Vector_Diff(Vector& vector1, Vector& vector2) {
    Vector result;
    result.x = vector1.x - vector2.x;
    result.y = vector1.y - vector2.y;
    result.z = vector1.z - vector2.z;
    return result;
}

Vector Cross(Vector& vector1, Vector& vector2) {
    Vector result;
    result.x = vector1.y * vector2.z - vector1.z * vector2.y;
    result.y = vector1.z * vector2.x - vector1.x * vector2.z;
    result.z = vector1.x * vector2.y - vector1.y * vector2.x;
    return result;
}

double Dot(Vector& vector1, Vector& vector2) {
    double result_x, result_y, result_z, result;
    result_x = vector1.x * vector2.x;
    result_y = vector1.y * vector2.y;
    result_z = vector1.z * vector2.z;
    result = result_x + result_y + result_z;

    return result;
}

double Triangle_Area(Point2& p0, Point2& p1, Point2& p2) {
    double result;
    result = 0.5 *
            ( (p0.x()*p1.y() - p1.x()*p0.y()) +
            (p1.x()*p2.y() - p2.x()*p1.y()) +
            (p2.x()*p0.y() - p0.x()*p2.y()) );
    return result;
}

double Area_Calc(Face& face, std::map<int, Vertex> vertices) {

//    std::vector<Point3> face_vertices;
//    for (const auto &vertex: face.tri_vertices) {
//        face_vertices.emplace_back(Point3(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z));
//    }
//    Plane best_fitting_plane;
//    CGAL::linear_least_squares_fitting_3(face_vertices.begin(), face_vertices.end(), best_fitting_plane,
//                                         CGAL::Dimension_tag<0>());
//    std::cout << "a: " << best_fitting_plane.a() << ", b: " << best_fitting_plane.b() << ", c: " << best_fitting_plane.c()
//              << ", d: " << best_fitting_plane.d() << std::endl;
//
//    std::vector<Point2> proj_points {};
//    for (const auto& pt_3: face_vertices)
//        proj_points.emplace_back(best_fitting_plane.to_2d(pt_3));
//    proj_points.emplace_back(proj_points[0]);
//
//    Point2 point_out (0.0, 0.0);
//    double area_out_pt, face_area = 0;
//    for (auto it = proj_points.begin(); it != proj_points.end() - 1; ++it) {
//        area_out_pt = Triangle_Area(point_out, *it, *(it+1));
//        face_area += area_out_pt;
//    }

    double s1, s2, s3;
    double f0_x = face.v3_coors[0].x(), f0_y = face.v3_coors[0].y(), f0_z = face.v3_coors[0].z();
    double f1_x = face.v3_coors[1].x(), f1_y = face.v3_coors[1].y(), f1_z = face.v3_coors[1].z();
    double f2_x = face.v3_coors[2].x(), f2_y = face.v3_coors[2].y(), f2_z = face.v3_coors[2].z();
    s1 = pow(pow((f0_x - f1_x), 2) + pow((f0_y - f1_y), 2) + pow((f0_z - f1_z), 2), 0.5);
    s2 = pow(pow((f1_x - f2_x), 2) + pow((f1_y - f2_y), 2) + pow((f1_z - f2_z), 2), 0.5);
    s3 = pow(pow((f2_x - f0_x), 2) + pow((f2_y - f0_y), 2) + pow((f2_z - f0_z), 2), 0.5);

    std::cout << "f0: " << f0_x << ", " << f0_y << ", " << f0_z << std::endl;
    std::cout << "f1: " << f1_x << ", " << f1_y << ", " << f1_z << std::endl;
    std::cout << "f2: " << f2_x << ", " << f2_y << ", " << f2_z << std::endl;

    std::cout << "s1: " << s1 << ", s2: " << s2 << ", s3: " << s3 << std::endl;

    double s = (s1 + s2 + s3) / 2.0;
    std::cout << "S: " << s << std::endl;
    double Area = pow((s * (s - s1) * (s - s2) * (s - s3)), 0.5);
    return Area;
}

double Volume_Calc(Face& face, std::map<int, Vertex> vertices) {
    Vector vector_a, vector_b, vector_c;
    Vector vector_ad, vector_bd, vector_cd;

    vector_a = Vertex_to_Vector(vertices[face.tri_vertices[0]]);
    vector_b = Vertex_to_Vector(vertices[face.tri_vertices[1]]);
    vector_c = Vertex_to_Vector(vertices[face.tri_vertices[2]]);

    Vector vector_out;
    vector_out.x = 1000; vector_out.y = 1000; vector_out.z = 1000;

    vector_ad = Vector_Diff(vector_a, vector_out);
    vector_bd = Vector_Diff(vector_b, vector_out);
    vector_cd = Vector_Diff(vector_c, vector_out);

    Vector cross_bd_cd;
    cross_bd_cd = Cross(vector_bd, vector_cd);

    double volume;
    Vector dot_a_bc;
    volume = (1.0/6.0) * Dot(vector_ad, cross_bd_cd);

    return volume;
}

std::pair<double, double> SurfArea_Volume_Object(std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
    //creating a vertex outside
    Vector vector_out;
    vector_out.x = 1000;
    vector_out.y = 1000;
    vector_out.z = 1000;

    double surface_area = 0, volume = 0;
    for (auto &[key, f]: faces) {
        std::cout << "face id: " << f.fid << " - v0: " << f.tri_vertices[0] << ", v1: " << f.tri_vertices[1] << ", v2: "
                  << f.tri_vertices[2] << std::endl;

        double face_area, face_vol;
        face_area = Area_Calc(f, vertices);
        face_vol = Volume_Calc(f, vertices);
        surface_area += abs(face_area);
        volume += face_vol;
        std::cout << "face id: " << f.fid << " - area: " << face_area << std::endl;
        std::cout << "face id: " << f.fid << " - volume: " << face_vol << std::endl;
        std::cout << "=============================" << std::endl;
    }

    return std::make_pair(surface_area, volume);
}

std::vector<double> Rectangularity(std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
    std::vector<int> checker;
    std::vector<Point3> object_points_3d;
    for (auto& [key, face] : faces) {
//        Point3 tv0, tv1, tv2;  // triangle vertices ids
        int tv0, tv1, tv2;  // triangle vertices ids
        tv0 = face.tri_vertices[0]; tv1 = face.tri_vertices[1]; tv2= face.tri_vertices[2];
//        tv0 = face.v3_coors[0]; tv0 = face.v3_coors[1]; tv0 = face.v3_coors[2];

//        if (std::find(object_points_3d.begin(), object_points_3d.end(), tv0) == object_points_3d.end()) {
//            object_points_3d.emplace_back(tv0);
//        }
//        if (std::find(object_points_3d.begin(), object_points_3d.end(), tv1) == object_points_3d.end()) {
//            object_points_3d.emplace_back(tv1);
//        }
//        if (std::find(object_points_3d.begin(), object_points_3d.end(), tv2) == object_points_3d.end()) {
//            object_points_3d.emplace_back(tv2);
//        }

        if (std::find(checker.begin(), checker.end(), tv0) == checker.end()) {
            // vertex 0 doesn't exist in checker, so we can add it
            checker.push_back(tv0);
            object_points_3d.emplace_back(face.v3_coors[0]);
        }
        if (std::find(checker.begin(), checker.end(), tv1) == checker.end()) {
            // vertex 0 doesn't exist in checker, so we can add it
            checker.push_back(tv1);
            object_points_3d.emplace_back(face.v3_coors[1]);
        }
        if (std::find(checker.begin(), checker.end(), tv2) == checker.end()) {
            // vertex 0 doesn't exist in checker, so we can add it
            checker.push_back(tv2);
            object_points_3d.emplace_back(face.v3_coors[2]);
        }
    }

    int i = 1;
    for (auto& obj_pt : object_points_3d) {
        std::cout << "obj pt: " << i << ", coor: " << obj_pt << std::endl;
        i++;
    }

    std::array<Point3, 8> obb_pts;
    CGAL::oriented_bounding_box(object_points_3d, obb_pts, CGAL::parameters::use_convex_hull(true));

    int j = 1;
    for (auto& obb_pt : obb_pts) {
        std::cout << "oobb pt " << j << ", coors: " << obb_pt << std::endl;
        j++;
    }

    double length, width, height;
    length = pow((pow(obb_pts[0].x() - obb_pts[1].x(), 2) +
                    pow(obb_pts[0].y() - obb_pts[1].y(), 2) +
                    pow(obb_pts[0].z() - obb_pts[1].z(), 2)), 0.5);
    width = pow((pow(obb_pts[1].x() - obb_pts[2].x(), 2) +
                    pow(obb_pts[1].y() - obb_pts[2].y(), 2) +
                    pow(obb_pts[1].z() - obb_pts[2].z(), 2)), 0.5);
    height = pow((pow(obb_pts[5].x() - obb_pts[0].x(), 2) +
                     pow(obb_pts[5].y() - obb_pts[0].y(), 2) +
                     pow(obb_pts[5].z() - obb_pts[0].z(), 2)), 0.5);
    double vol_oobb = length * width * height;

    std::cout << "length: " << length << std::endl;
    std::cout << "width: " << width << std::endl;
    std::cout << "height: " << height << std::endl;
    std::vector<double> lwh ;
    lwh.emplace_back(length); lwh.emplace_back(width); lwh.emplace_back(height);
    return lwh;

}


double hemisphericality(const double& volume, const double& surface_area) {
    double hemisphericality;
    hemisphericality = (3.0 * pow(2.0 * M_PI, 0.5) ) * volume / pow(surface_area, 1.5);
    return hemisphericality;
}


double roughness_index(std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
    std::vector<Point3> object_points;
    std::vector<int> v_checker;
    int count_rand_pts = 1;  // random points per face
    int num_randPts_per_face = 2;
    std::vector<Point3> sample_pts_on_surface;
    CGAL::Random rand;
    for (auto& [key, face] : faces ) {
        for (auto& vert_coor : face.v3_coors) {
            auto it = std::find(object_points.begin(), object_points.end(), vert_coor);
            if (it == object_points.end()) {
                // Element does not exist, add it to the vector
                object_points.emplace_back(vert_coor);
                sample_pts_on_surface.emplace_back(vert_coor);
            }
        }

        // best fitting plane and 3 projected vertices a, b, c on to best_FP
        Plane best_FP;
        CGAL::linear_least_squares_fitting_3(face.v3_coors.begin(), face.v3_coors.end(), best_FP,CGAL::Dimension_tag<0>());
        Point2 a = best_FP.to_2d(face.v3_coors[0]);
        Point2 b = best_FP.to_2d(face.v3_coors[1]);
        Point2 c = best_FP.to_2d(face.v3_coors[2]);
        std::cout << "=============" << std::endl;
        std::cout << "face - " << key << std::endl;
        std::cout << "object points size: " << object_points.size() << std::endl;
        std::cout << "point 2D a: " << a << std::endl;
        std::cout << "point 3D a: " << face.v3_coors[0]<< std::endl;
        std::cout << "point 2D b: " << b << std::endl;
        std::cout << "point 3D b: " << face.v3_coors[1]<< std::endl;
        std::cout << "point 2D c: " << c << std::endl;
        std::cout << "point 3D c: " << face.v3_coors[2]<< std::endl;

        for (int rand_pt=0; rand_pt<num_randPts_per_face; rand_pt++) {
            double r1 = std::sqrt(rand.get_double());
            double r2 = rand.get_double();
            Point2 p ((1 - r1) * a.x() + (r1 * (1 - r2)) * b.x() + (r1 * r2) * c.x(),
                    (1 - r1) * a.y() + (r1 * (1 - r2)) * b.y() + (r1 * r2) * c.y());
            Point3 p_in3d {best_FP.to_3d(p)};  // projecting the sample point inside the triangle back to 3D
            sample_pts_on_surface.emplace_back(p_in3d);
            std::cout << "rand point 2D: " << p << std::endl;
            std::cout << "rand point 3D: " << p_in3d << std::endl;
        }
        std::cout << "=============" << std::endl;
    }

    std::cout << "surface sample points size: " << sample_pts_on_surface.size() << std::endl;


    double cx {0}, cy {0}, cz {0};
    for (auto& pt : object_points) {
        cx += pt.x(); cy += pt.y(); cz += pt.z();
    }

    double num_vertices = object_points.size();
    Point3 centroid_E (cx / num_vertices, cy / num_vertices, cz / num_vertices);

    double samples_dist {0}, dist {0}, squared_dist {0};
    int i {1};
    for (auto& surf_pt : sample_pts_on_surface) {
        squared_dist =  pow(centroid_E.x()-surf_pt.x(), 2) +
                        pow(centroid_E.y()-surf_pt.y(), 2) +
                        pow(centroid_E.z()-surf_pt.z(), 2);
        dist = pow(squared_dist, 0.5);
        samples_dist += dist;
        std::cout << "sample dist " << i << ": " << dist << std::endl;
        i++;
    }

    double mu {samples_dist/sample_pts_on_surface.size()};
    double rid {0};  // initiating roughness index to 0

    std::pair<double, double> SurfArea_Vol_obj = SurfArea_Volume_Object(faces, vertices);
    double surf_area = SurfArea_Vol_obj.first;
    double volume = SurfArea_Vol_obj.second;

    rid = 48.735 * pow(mu, 3) / (volume + pow(surf_area, 1.5));

    return rid;
}

int main(int argc, const char * argv[]) {

    std::map<int, Vertex> vertices;
    std::map<int, Face> faces;

    // Reading the obj file
    const char *inobj = (argc > 1) ? argv[1] : "../data/torus_out.obj";
    std::ifstream input_stream;
    input_stream.open(inobj);

    std::string line;

    //creating a vertex outside
    Vector vector_out {};
    vector_out.x = 1000;
    vector_out.y = 1000;
    vector_out.z = 1000;


    // Reading lines from obj
    int i = 1;
    while (std::getline(input_stream, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {  // reading the lines only with v
            float x, y, z;
            ss >> x >> y >> z;
//            std::cout << "vertex i: " << i << " x: " << x << " y: " << y << " z: " << z << std::endl;

            vertices[i].id = i;
            vertices[i].x = x;
            vertices[i].y = y;
            vertices[i].z = z;
            i++;
        }
        // break if a new object is identified
        if (prefix == "o")
            break;
    }
    // checking if while is continuing to read the remaining lines
    int faceid = 0;
    while (std::getline(input_stream, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "f") {
            int v0, v1, v2;
            ss >> v0 >> v1 >> v2;
            faces[faceid].fid = faceid;
            faces[faceid].tri_vertices.emplace_back(v0);
            faces[faceid].tri_vertices.emplace_back(v1);
            faces[faceid].tri_vertices.emplace_back(v2);
            faces[faceid].v3_coors.emplace_back(Point3(vertices[v0].x, vertices[v0].y, vertices[v0].z));
            faces[faceid].v3_coors.emplace_back(Point3(vertices[v1].x, vertices[v1].y, vertices[v1].z));
            faces[faceid].v3_coors.emplace_back(Point3(vertices[v2].x, vertices[v2].y, vertices[v2].z));
        }
        faceid++;

        // break if a new object is identified
        if (prefix == "o")
            break;
    }

    // checking if vertices are loded or not
    for (auto &[key, vertex]: vertices) {
        std::cout << "vertex id: " << key << " - x: " << vertex.x << " y: " << vertex.y << " z: " << vertex.z
                  << std::endl;
    }
    std::cout << "=============================" << std::endl;
    // checking if all the faces are loaded or not
    // calculating volume made by every triangle with outside vertex
    double surface_area, volume;
    std::pair<double, double> SurfArea_Vol_obj = SurfArea_Volume_Object(faces, vertices);
    surface_area = SurfArea_Vol_obj.first;
    volume = SurfArea_Vol_obj.second;
    std::cout << "Area of object: " << surface_area << std::endl;
    std::cout << "Volume of object: " << volume << std::endl;
    std::cout << "Hemisphericality of object: " << hemisphericality(volume, surface_area) << std::endl;


    // ROUGHNESS INDEX
    double rough_index;
    rough_index = roughness_index(faces, vertices);
//    std::cout << "rough_index of object: x - " << rough_index.x() << ", y - " << rough_index.y() << ", z - " << rough_index.z() << std::endl;
    std::cout << "roughness index of object: " << rough_index << std::endl;

    // Rectangularity Calc
    std::cout << "======================" << std::endl;
//    std::array<Point3, 8> rectangularity;
    std::vector<double> rectangularity;
//    std::vec<double rectangularity;
    rectangularity = Rectangularity(faces, vertices);
    int j = 1;
    for (auto& obb_pt: rectangularity) {
        std::cout << "side: " << obb_pt << std::endl;
        j++;
    }
//    std::cout << "volume_oobb: " << rectangularity << std::endl;
    return 0;
}