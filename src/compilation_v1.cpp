//
// Created by Sharath Chandra on 20-Mar-23.
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


std::map<int, Vertex> get_coordinates(const json& j, bool translate) {
    std::map<int, Vertex> lspts;
    std::vector<std::vector<int>> lvertices = j["vertices"];
    int i = 1;
    if (translate) {
        for (auto& vi : lvertices) {
            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
            lspts[i].id = i; lspts[i].x = x; lspts[i].y = y; lspts[i].z = z;
            i++;
        }
    } else {
        //-- do not translate, useful to keep the values low for downstream processing of data
        for (auto& vi : lvertices) {
            double x = (vi[0] * j["transform"]["scale"][0].get<double>());
            double y = (vi[1] * j["transform"]["scale"][1].get<double>());
            double z = (vi[2] * j["transform"]["scale"][2].get<double>());
            lspts[i].id = i; lspts[i].x = x; lspts[i].y = y; lspts[i].z = z;
            i++;
        }
    }
    return lspts;
}


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


double Rectangularity(double& obj_volume, std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
    std::vector<int> checker;
    std::vector<Point3> object_points_3d;
    for (auto& [key, face] : faces) {
        int tv0, tv1, tv2;  // triangle vertices ids
        tv0 = face.tri_vertices[0]; tv1 = face.tri_vertices[1]; tv2= face.tri_vertices[2];
        if (std::find(checker.begin(), checker.end(), tv0) == checker.end()) {
            // vertex 0 doesn't exist in checker, so we can add it
            checker.push_back(tv0);
            object_points_3d.emplace_back(face.v3_coors[0]);
        }
        if (std::find(checker.begin(), checker.end(), tv1) == checker.end()) {
            // vertex 1 doesn't exist in checker, so we can add it
            checker.push_back(tv1);
            object_points_3d.emplace_back(face.v3_coors[1]);
        }
        if (std::find(checker.begin(), checker.end(), tv2) == checker.end()) {
            // vertex 2 doesn't exist in checker, so we can add it
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
    double oobb_vol = length * width * height;

    std::cout << "length: " << length << std::endl;
    std::cout << "width: " << width << std::endl;
    std::cout << "height: " << height << std::endl;

    double rectangularity = obj_volume / oobb_vol;
    return rectangularity;
}




int main(int argc, const char * argv[]) {

    std::map<int, Vertex> vertices;
    std::map<int, Face> faces;

    // Reading the json file to check the number of objects
    const char* injson = (argc > 1) ? argv[1] : "../data/2b.city.json";
    std::ifstream input(injson);
    json j;
    input >> j; //-- store the content of the file in a nlohmann::json object
    input.close();

    //outfile
    std::string outfile = "../data/2b_v1.obj";
    std::ofstream ofile(outfile);

    // storing the vertices
    vertices = get_coordinates(j, true);
    std::vector<Point3> lspts;
    std::cout << "number of vertices: " << vertices.size() << std::endl;
    for (auto& [vid, vert] : vertices) {
//        std::cout << "v: " << vid << ", x: " << vert.x << ", y: " << vert.y << ", z: " << vert.z << std::endl;
        ofile << std::setprecision(5) << std::fixed << "v " << vert.x << " " << vert.y << " " << vert.z << " " << std::endl;
        lspts.emplace_back(Point3(vert.x, vert.y, vert.z));
    }

    std::vector<std::string> objects_ids;  // objects names which are solids - those which contain boundaries and semantics
    int i {1};
    for (auto& co : j["CityObjects"].items()) {
        std::cout << "City Object - " << ": " << co.key() << std::endl;
        faces.clear();
        int f = 1;
        int simple_roofs = 0;
        for (auto& obj_info :  co.value().items()) {
            std::cout << "\tkey: " << obj_info.key() << std::endl;
            if (obj_info.key() == "geometry") {

                for (auto& g : obj_info.value()) {

//                    for ( auto& geom_dict : g.items()) {
//
//                        // looping through all the geometry
//                        std::cout << "\t\t" << geom_dict.key() << std::endl;
//                        if (geom_dict.key() == "type" || geom_dict.key() == "lod") {
//                            std::cout << "\t\t\t" << geom_dict.value() << std::endl;
//                        }
//                        if (geom_dict.key() == "semantics") {
//                            for (auto& seman_dict : geom_dict.value().items()) {
//                                std::cout << "\t\t\t" << seman_dict.key() << std::endl;
//                                if (seman_dict.key() == "surfaces") {
//                                    std::cout << "\t\t\t\t" << seman_dict.value() << std::endl;
//                                }
//                                if (seman_dict.key() == "values") {
//                                    std::cout << "\t\t\t\tnum faces: " << seman_dict.value()[0].size() << std::endl;
//                                    std::cout << "\t\t\t\t" << seman_dict.value() << std::endl;
//                                    // roof surface has value 1
//                                    int surf_counter = 0; int searchValue = 1;
//                                    std::vector<int> roof_Surfaces_idx;
//                                    for (auto& surfValue : seman_dict.value()[0]) {
//                                        if (surfValue == searchValue) {
//                                            roof_Surfaces_idx.emplace_back(surf_counter);
//                                        }
//                                        surf_counter++;
//                                    }
//                                    std::cout << "\t\t\t\t\tnumber of roofs: " << roof_Surfaces_idx.size() << std::endl;
//                                    std::cout << "\t\t\t\t\troof surface idxs: ";
//                                    for (auto& roof_num : roof_Surfaces_idx) {
//                                        std:: cout << roof_num << " ";
//                                    }
//                                    std::cout << std::endl;
//                                }
//                            }
//                        }
//
//                        if (geom_dict.key() == "boundaries") {
//                            std::cout << "\t\t\tnum faces: " << geom_dict.value()[0].size() << std::endl;
//                        }
//
//                    }

                    if ( (g["type"] == "Solid") && (g["lod"] == "2.2")) {
                        ofile << "o " << co.key() << std::endl;
                        std::cout <<"o " << co.key() << std::endl;
                        std::cout << "num of boundaries: " << g["boundaries"][0].size() << std::endl;
                        for (int i = 0; i < g["boundaries"].size(); i++) {
                            std::cout << "i : " << i << std::endl;
//                            std::cout << "g[\"boundaries\": ]" << g["boundaries"][i] << std::endl;
                            for (int j = 0; j < g["boundaries"][i].size(); j++) {
//                                std:: cout << "i: " << i << " j: " << j << std::endl;
                                std::vector<std::vector<int>> gb = g["boundaries"][i][j];
                                std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
                                for (auto& tr : trs) {
//                                    ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
//                                    std::cout << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
                                    faces[f].fid = f;
                                    faces[f].tri_vertices.emplace_back(tr[0] + 1);
                                    faces[f].tri_vertices.emplace_back(tr[1] + 1);
                                    faces[f].tri_vertices.emplace_back(tr[2] + 1);
                                    faces[f].v3_coors.emplace_back(Point3 (vertices[tr[0] + 1].x, vertices[tr[0] + 1].y, vertices[tr[0] + 1].z));
                                    faces[f].v3_coors.emplace_back(Point3 (vertices[tr[1] + 1].x, vertices[tr[1] + 1].y, vertices[tr[1] + 1].z));
                                    faces[f].v3_coors.emplace_back(Point3 (vertices[tr[2] + 1].x, vertices[tr[2] + 1].y, vertices[tr[2] + 1].z));
                                    f++;
                                }
                            }
                        }
                        std::cout << "\tsemantic surfaces: " << g["semantics"]["surfaces"] << std::endl;
                        std::cout << "\tsemantic values: " << g["semantics"]["values"] << std::endl;
                        int surf_counter = 0; int searchValue = 1;
                        std::vector<int> roof_Surfaces_idx;
                        for (auto& surf_type_value : g["semantics"]["values"][0]) {
//                            std::cout << surf_type_value << std::endl;
                            if (surf_type_value == searchValue) {
                                roof_Surfaces_idx.emplace_back(surf_counter);
                            }
                            surf_counter++;
                        }

                        std::cout << "\t\tnum of roof surfaces: " << roof_Surfaces_idx.size() << std::endl;
                        std::cout << "\t\troof surfaces indexes: ";
                        for (auto& roof_surf_check : roof_Surfaces_idx) {
                            std::cout << roof_surf_check << " ";
                        }
                        std::cout << std::endl;

                        /////////
                        for (auto& roof_surf_id : roof_Surfaces_idx) {
                            std::cout << "\t\t\troof id: " << roof_surf_id << " - " << g["boundaries"][0][roof_surf_id] << std::endl;
                            // make the best fit plane and making a polygon to make sure the orientation is correct
                            Polygon2 pgn;
                            std::vector<Point3> Outer_ring_pt3;
                            for (auto& oring_id : g["boundaries"][0][roof_surf_id][0]) {
                                Outer_ring_pt3.emplace_back(Point3(vertices[oring_id].x, vertices[oring_id].y, vertices[oring_id].z));
//                                std::cout << "\t\t\t\tpoint3: " << Point3(vertices[oring_id].x, vertices[oring_id].y, vertices[oring_id].z) << std::endl;
                                }
                            std::cout << "\t\t\t\tnum outer vertices: " << Outer_ring_pt3.size() << std::endl;
                            Plane best_fit_plane;
                            CGAL::linear_least_squares_fitting_3(Outer_ring_pt3.begin(), Outer_ring_pt3.end(),
                                                                 best_fit_plane, CGAL::Dimension_tag<0>());
                            std::cout << "\t\t\t\ta, b, c, d: " << best_fit_plane << std::endl;
                            K::Vector_3 normal = best_fit_plane.orthogonal_vector();
                            std::cout << "\t\t\t\tnormal: " << normal << std::endl;

                            if (normal.z() < 0) {
                                std::cout << "\t\t\t\tNORMAL IS INVERTED" << std::endl;
                                normal = normal * -1;
                                std::cout << "\t\t\t\tnew normal: " << normal << std::endl;
                            }

//                            Point3 pt_above_surf (Outer_ring_pt3[0].x(), Outer_ring_pt3[0].y(), Outer_ring_pt3[0].z()-22.0);
//
//                            double normal_product;
//                            bool normal_orientation;
//                            normal_product = pt_above_surf.x() * normal.x() +
//                                                 pt_above_surf.y() * normal.y() +
//                                                 pt_above_surf.z() * normal.z();
//                            if (normal_product > 0)
//                                normal_orientation = true;
//                            else
//                                normal_orientation = false;
//                            std::cout << "\t\t\t\tnormal orientation: " << normal_orientation << std::endl;
//                            // pushing points to the polugon to make sure of its orientation is correct
//                            for (auto& vert3d : Outer_ring_pt3) {
//                                pgn.push_back(best_fit_plane.to_2d(vert3d));
//                            }
//
//                            pgn.push_back(best_fit_plane.to_2d(Outer_ring_pt3[0]));
//                            std::cout << "\t\t\t\tpolygon size and area: " << pgn.size() << " " << pgn.area() << std::endl;
//
////                            bool is_simple = CGAL::is_simple_2(pgn.vertices_begin(), pgn.vertices_end(), K());
//                            bool is_simple = pgn.is_simple();
//
//                            std::cout << "\t\t\t\tis simple: " << is_simple << std::endl;
//                            if (is_simple == true) {
//                                std::cout << "\t\t\t\torientation: " << pgn.orientation() << std::endl;
//                                simple_roofs++;
//                            }

//                            CGAL::Orientation roof_orientation = pgn.orientation();
//                            if (roof_orientation == CGAL::CLOCKWISE) {
//                                std::cout << "\t\t\t\tThe polygon has clockwise orientation." << std::endl;
//                            } else if (roof_orientation == CGAL::COUNTERCLOCKWISE) {
//                                std::cout << "\t\t\t\tThe polygon has counterclockwise orientation." << std::endl;
//                            } else {
//                                std::cout << "\t\t\t\tThe polygon is degenerate." << std::endl;
//                            }


//                            if (pgn.is_counterclockwise_oriented() == true)
//                                std::cout << "\t\t\t\torientation CCW" << std::endl;
//                            else
//                                std::cout << "\t\t\t\torientation CW" << std::endl;

                        }


                    }

                }

            }
            // else not required
//            else {
//                std::cout << "\t\tvalue: " << obj_info.value() << std::endl;
//            }
        }
        for (auto& [fkey, face] : faces) {
            ofile << "f " << face.tri_vertices[0] << " " << face.tri_vertices[1] << " " << face.tri_vertices[2] << std::endl;
        }

        std::cout << "\tnum of simple roofs: " << simple_roofs << std::endl;

    }
    ofile.close();

    // Reading the obj file
//    const char* inobj = (argc > 1) ? argv[1] : "../data/2b.obj";
//    std::ifstream input_stream;
//    input_stream.open(inobj);
//
//    std::string line;
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
//            std::cout << "vertex i: " << i << " x: " << x << " y: " << y << " z: " << z << std::endl;
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




    return 0;
}