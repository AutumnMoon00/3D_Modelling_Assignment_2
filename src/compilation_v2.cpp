//
// Created by Sharath Chandra on 22-Mar-23.
//

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <CGAL/Random.h>
#include <typeinfo>

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


json surf_adder(json& semantic_surfaces) {

    std::vector<json> new_surfaces;

    json ne_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "NorthEast"}
    };
    new_surfaces.emplace_back(ne_surf);

    json en_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "EastNorth"}
    };
    new_surfaces.emplace_back(en_surf);

    json es_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "EastSouth"}
    };
    new_surfaces.emplace_back(es_surf);

    json se_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "SouthEast"}
    };
    new_surfaces.emplace_back(se_surf);

    json sw_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "SouthWest"}
    };
    new_surfaces.emplace_back(sw_surf);

    json ws_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "WestSouth"}
    };
    new_surfaces.emplace_back(ws_surf);

    json wn_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "WestNorth"}
    };
    new_surfaces.emplace_back(wn_surf);

    json nw_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "NorthWest"}
    };
    new_surfaces.emplace_back(nw_surf);

    json horz_surf = {
            {"type", "RoofSurface"},
            {"Orientation", "Horizontal"}
    };
    new_surfaces.emplace_back(horz_surf);

    for (auto& surf : new_surfaces) {
        semantic_surfaces.emplace_back(surf);
    }

    return semantic_surfaces;
}


std::pair<double, double> roof_elevation_azimuth(K::Vector_3& normal) {
    double factor = std::pow(10.0, 1);
    double x = normal.x(), y = normal.y(), z = normal.z();

    double elevation_deg = atan2(z, sqrt(x*x + y*y)) * 180 / M_PI;
    elevation_deg = std::round (elevation_deg * factor) / factor;  // rounding to 1 decimal

    double azimuth_deg;
    if (x > 0) {
        azimuth_deg = atan2(x, y) * 180 / M_PI ;
        azimuth_deg = std::round (factor * azimuth_deg) / factor;
    } else if (x < 0) {
        azimuth_deg = (atan2(x, y) + 2 * M_PI) * 180 / M_PI;
        azimuth_deg = std::round (factor * azimuth_deg) / factor;
    } else if (y > 0) {
        azimuth_deg = 0.0;
    } else if (y < 0) {
        azimuth_deg = 180.0;
    } else {
        // vector is at the origin
        azimuth_deg = 0.0;
    }

    return std::make_pair(elevation_deg, azimuth_deg);
}

void orientation_setter(json& g, std::map<int, Vertex>& vertices) {
    g["semantics"]["surfaces"] = surf_adder(g["semantics"]["surfaces"]);
    std::cout << "\tsemantics surfces: " << g["semantics"]["surfaces"] << std::endl;
    std::cout << "\tsemantics values: " << g["semantics"]["values"] << std::endl;

    int surf_counter = 0; int searchValue = 1;
    for (auto& surf_value : g["semantics"]["values"][0]) {
        if (surf_value == 1) {
            std::vector<Point3> Outer_ring_pt3;
            Plane best_fit_plane;
            std::cout << "\t\tsurf oring size: " << g["boundaries"][0][surf_counter][0].size();
            for (auto& oring_id : g["boundaries"][0][surf_counter][0]) {
                json result = oring_id.get<int>() + 1;
                Outer_ring_pt3.emplace_back(Point3 (vertices[result].x, vertices[result].y, vertices[result].z));
            }
            CGAL::linear_least_squares_fitting_3(Outer_ring_pt3.begin(), Outer_ring_pt3.end(),
                                                 best_fit_plane, CGAL::Dimension_tag<0>());
            std::cout << "\t\troof surf id: " << surf_counter << ": size: " << Outer_ring_pt3.size() << std::endl;
            std::cout << "\t\t\ta, b, c, d: " << best_fit_plane << std::endl;
            K::Vector_3 normal = best_fit_plane.orthogonal_vector();
            std::cout << "\t\t\tnormal: " << normal << std::endl;
            if (normal.z() < 0) {
                std::cout << "\t\t\tNORMAL IS INVERTED" << std::endl;
                normal = normal * -1;
                std::cout << "\t\t\tnew normal: " << normal << std::endl;
            }
            std::pair<double, double> elevation_azimuth = roof_elevation_azimuth(normal);
            double roof_elevation = elevation_azimuth.first;
            double roof_azimuth = elevation_azimuth.second;
            std::cout << "\t\t\troof_elevation: " << roof_elevation << std::endl;
            std::cout << "\t\t\troof azimuth: " << roof_azimuth << std::endl;

            if (roof_elevation == 0.0)
                g["semantics"]["values"][0][surf_counter] = 12;
            else if (0 <= roof_azimuth && roof_azimuth < 45)  // northeast
                g["semantics"]["values"][0][surf_counter] = 4;
            else if (45 <= roof_azimuth && roof_azimuth < 90)  // eastnorth
                g["semantics"]["values"][0][surf_counter] = 5;
            else if (90 <= roof_azimuth && roof_azimuth < 135)  // eastsouth
                g["semantics"]["values"][0][surf_counter] = 6;
            else if (135 <= roof_azimuth && roof_azimuth < 180)  // southeast
                g["semantics"]["values"][0][surf_counter] = 7;
            else if (180 <= roof_azimuth && roof_azimuth < 225)  // southwest
                g["semantics"]["values"][0][surf_counter] = 8;
            else if (225 <= roof_azimuth && roof_azimuth < 270)  // westsouth
                g["semantics"]["values"][0][surf_counter] = 9;
            else if (270 <= roof_azimuth && roof_azimuth < 315)  // westnorth
                g["semantics"]["values"][0][surf_counter] = 10;
            else if (315 <= roof_azimuth && roof_azimuth < 360)  // northwest
                g["semantics"]["values"][0][surf_counter] = 11;

            std::cout << "\t\t\torientation: " << g["semantics"]["values"][0][surf_counter] << std::endl;
            int orient = g["semantics"]["values"][0][surf_counter];
            std::cout << "\t\t\troof orientation: " << g["semantics"]["surfaces"][orient] << std::endl;
//            std::cout << "\t\t\troof orientation: " << g["semantics"]["surfaces"][g["semantics"]["values"][0][surf_counter]] << std::endl;

        }
        surf_counter++;
    }
    std::cout << "\tsemantics values: " << g["semantics"]["values"] << std::endl;
}

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

//    std::cout << "f0: " << f0_x << ", " << f0_y << ", " << f0_z << std::endl;
//    std::cout << "f1: " << f1_x << ", " << f1_y << ", " << f1_z << std::endl;
//    std::cout << "f2: " << f2_x << ", " << f2_y << ", " << f2_z << std::endl;

//    std::cout << "s1: " << s1 << ", s2: " << s2 << ", s3: " << s3 << std::endl;

    double s = (s1 + s2 + s3) / 2.0;
//    std::cout << "S: " << s << std::endl;
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
//        std::cout << "face id: " << f.fid << " - v0: " << f.tri_vertices[0] << ", v1: " << f.tri_vertices[1] << ", v2: "
//                  << f.tri_vertices[2] << std::endl;

        double face_area, face_vol;
        face_area = Area_Calc(f, vertices);
        face_vol = Volume_Calc(f, vertices);
        surface_area += abs(face_area);
        volume += face_vol;
//        std::cout << "face id: " << f.fid << " - area: " << face_area << std::endl;
//        std::cout << "face id: " << f.fid << " - volume: " << face_vol << std::endl;
//        std::cout << "=============================" << std::endl;
    }

    return std::make_pair(surface_area, abs(volume));
}

//std::vector<double>
double Rectangularity(double& volume_object, std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
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

    std::array<Point3, 8> obb_pts;
    CGAL::oriented_bounding_box(object_points_3d, obb_pts, CGAL::parameters::use_convex_hull(true));

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

    std::cout << "\tlength: " << length << std::endl;
    std::cout << "\twidth: " << width << std::endl;
    std::cout << "\theight: " << height << std::endl;

    double rectangularity = volume_object / vol_oobb;
    return rectangularity;
}


double Roughness_Index(double& volume_object, double& surface_area_object, std::map<int, Face>& faces, std::map<int, Vertex>& vertices) {
    std::vector<Point3> object_points;
    std::vector<int> v_checker;
    int count_rand_pts = 1;  // random points per face
    int num_randPts_per_face = 3;
    std::vector<Point3> sample_pts_on_surface;
    CGAL::Random rand;
    std::map<Point3, double> face_centroid_area;  // face centroid point3, and its area
    for (auto& [key, face] : faces ) {
        double cent_x = 0, cent_y = 0, cent_z = 0;
        for (auto& vert_coor : face.v3_coors) {
            auto it = std::find(object_points.begin(), object_points.end(), vert_coor);
            if (it == object_points.end()) {
                // Element does not exist, add it to the vector
                object_points.emplace_back(vert_coor);
                sample_pts_on_surface.emplace_back(vert_coor);
            }
            cent_x += vert_coor.x(); cent_y += vert_coor.y(); cent_z += vert_coor.z();
        }

        cent_x = cent_x/3.0; cent_y = cent_y/3.0; cent_z = cent_z/3.0;
        Point3 centroid_face (cent_x, cent_y, cent_z);
        double face_area = Area_Calc(face, vertices);
        face_centroid_area[centroid_face] = face_area;
//        std::cout << "\t\tface centroid: " << centroid_face << std::endl;
//        std::cout << "\t\tface area: " << face_centroid_area[centroid_face] << std::endl;


        // best fitting plane and 3 projected vertices a, b, c on to best_FP
        Plane best_FP;
        CGAL::linear_least_squares_fitting_3(face.v3_coors.begin(), face.v3_coors.end(), best_FP,
                                             CGAL::Dimension_tag<0>());
        Point2 a = best_FP.to_2d(face.v3_coors[0]);
        Point2 b = best_FP.to_2d(face.v3_coors[1]);
        Point2 c = best_FP.to_2d(face.v3_coors[2]);

        for (int rand_pt = 0; rand_pt < num_randPts_per_face; rand_pt++) {
            double r1 = std::sqrt(rand.get_double());
            double r2 = rand.get_double();
            Point2 p((1 - r1) * a.x() + (r1 * (1 - r2)) * b.x() + (r1 * r2) * c.x(),
                     (1 - r1) * a.y() + (r1 * (1 - r2)) * b.y() + (r1 * r2) * c.y());
            Point3 p_in3d{best_FP.to_3d(p)};  // projecting the sample point inside the triangle back to 3D
            sample_pts_on_surface.emplace_back(p_in3d);
        }
    }

    double object_centroid_x = 0, object_centroid_y = 0, object_centroid_z = 0;
    double total_area = 0;
    for (auto& [cent_pt3, tri_area] : face_centroid_area) {
        object_centroid_x += cent_pt3.x() * tri_area;
        object_centroid_y += cent_pt3.y() * tri_area;
        object_centroid_z += cent_pt3.z() * tri_area;
        total_area += tri_area;
    }

    Point3 Object_Centriod (object_centroid_x/total_area, object_centroid_y/total_area, object_centroid_z/total_area);

    double samples_dist {0}, dist {0}, squared_dist {0};
    int i {1};
    for (auto& surf_pt : sample_pts_on_surface) {
        squared_dist =  pow(Object_Centriod.x()-surf_pt.x(), 2) +
                        pow(Object_Centriod.y()-surf_pt.y(), 2) +
                        pow(Object_Centriod.z()-surf_pt.z(), 2);
        dist = pow(squared_dist, 0.5);
        samples_dist += dist;
//        std::cout << "sample dist " << i << ": " << dist << std::endl;
        i++;
    }

    double mu {samples_dist/sample_pts_on_surface.size()};
    double rid = 48.735 * pow(mu, 3) / (volume_object + pow(surface_area_object, 1.5));
    return rid;

}


double Hemisphericality(const double& volume_object, const double& surface_area_object) {
    double hemisphericality;
    hemisphericality = (3.0 * pow(2.0 * M_PI, 0.5) ) * volume_object / pow(surface_area_object, 1.5);
    return hemisphericality;
}


int main(int argc, const char * argv[]) {

    double factor = std::pow(10.0, 1);

    std::map<int, Vertex> vertices;
    std::map<int, Face> faces;

    // Reading the json file to check the number of objects
    const char* injson = (argc > 1) ? argv[1] : "../data/2b.city.json";
    std::ifstream input(injson);
    json j;
    input >> j; //-- store the content of the file in a nlohmann::json object
    input.close();

    //output object file
    std::string outfile = "../output/2b.obj";
    std::ofstream ofile(outfile);

    // storing the vertices
    vertices = get_coordinates(j, true);
    std::vector<Point3> lspts;
    std::cout << "number of vertices: " << vertices.size() << std::endl;
    for (auto& [vid, vert] : vertices) {
        std::cout << "v: " << vid << ", x: " << vert.x << ", y: " << vert.y << ", z: " << vert.z << std::endl;
        ofile << std::setprecision(5) << std::fixed << "v " << vert.x << " " << vert.y << " " << vert.z << " " << std::endl;
        lspts.emplace_back(Point3(vert.x, vert.y, vert.z));
    }


    // starting with objects and faces
    std::vector<std::string> objects_ids;  // objects names which are solids - those which contain boundaries and semantics
//    int i {1};
    for (auto& co : j["CityObjects"].items()) {
        std::cout << "\nCity Object - " << ": " << co.key() << std::endl;
        faces.clear();
        int f = 1;  // starting face ids from 1 for every new object
        for (auto &g: co.value()["geometry"]) {
            if ((g["type"] == "Solid") && (g["lod"] == "2.2")) {   //-- LoD2.2 only!!!!!
                ofile << "o " << co.key() << std::endl;
//                std::cout << "o " << co.key() << std::endl;
                std::cout << "\tnum of boundaries: " << g["boundaries"][0].size() << std::endl;
                for (int i = 0; i < g["boundaries"].size(); i++) {
                    for (int j = 0; j < g["boundaries"][i].size(); j++) {
                        std::vector<std::vector<int>> gb = g["boundaries"][i][j];
                        std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
                        for (auto &tr: trs) {
                            ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
//                            std::cout << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
                            faces[f].fid = f;
                            faces[f].tri_vertices.emplace_back(tr[0] + 1);
                            faces[f].tri_vertices.emplace_back(tr[1] + 1);
                            faces[f].tri_vertices.emplace_back(tr[2] + 1);
                            faces[f].v3_coors.emplace_back(
                                    Point3(vertices[tr[0] + 1].x, vertices[tr[0] + 1].y, vertices[tr[0] + 1].z));
                            faces[f].v3_coors.emplace_back(
                                    Point3(vertices[tr[1] + 1].x, vertices[tr[1] + 1].y, vertices[tr[1] + 1].z));
                            faces[f].v3_coors.emplace_back(
                                    Point3(vertices[tr[2] + 1].x, vertices[tr[2] + 1].y, vertices[tr[2] + 1].z));
                            f++;
                        }
                    }
                }

                orientation_setter(g, vertices);

                std::pair<double, double> surf_area_volume = SurfArea_Volume_Object(faces, vertices);
                double surface_area_object = surf_area_volume.first;
                double volume_object = surf_area_volume.second;
                double rectangularity_object = Rectangularity(volume_object, faces, vertices);
                double hemisphericality_object = Hemisphericality(volume_object, surface_area_object);
                double roughness_index_object = Roughness_Index(volume_object, surface_area_object, faces, vertices);

                std::cout << "\tvolume of object: " << volume_object << std::endl;
                std::cout << "\tarea of object: " << surface_area_object << std::endl;
                std::cout << "\trectangularity of object: " << rectangularity_object << std::endl;
                std::cout << "\themisphericality of object: " << hemisphericality_object << std::endl;
                std::cout << "\troughness index of object: " << roughness_index_object << std::endl;

                std::cout << "\ttype id of object name: " << typeid(co.key()).name() << std::endl;
                std::string parent_obj_name = co.value()["parents"][0];
                std::cout << "\tparent name: " << parent_obj_name << std::endl;

                j["CityObjects"][parent_obj_name]["attributes"]["Volume"] = volume_object;
                j["CityObjects"][parent_obj_name]["attributes"]["Rectangularity"] = rectangularity_object;
                j["CityObjects"][parent_obj_name]["attributes"]["Hemisphericality"] = hemisphericality_object;
                j["CityObjects"][parent_obj_name]["attributes"]["RoughnessIndex"] = roughness_index_object;

                std::cout << "parent volume: " << j["CityObjects"][co.key().substr(0, co.key().size() - 2)]["attributes"]["Volume"] << std::endl;
            }

        }
    }

    ofile.close();

    std::string out_json_name = "../output/2b_updated.city.json";
    std::ofstream o(out_json_name);
    o << j.dump(2) << std::endl;
    o.close();

    std::cout << "=========================================" << std::endl;
    for (auto& co : j["CityObjects"].items()) {
        std::cout << "\nCity Object - " << ": " << co.key() << std::endl;
//        std::cout << "\tsemantics: " << co.value() << std::endl;
        std::cout << "\tvolume: " << co.value()["attributes"]["Volume"] << std::endl;
        std::cout << "\tRectangularity: " << co.value()["attributes"]["Rectangularity"] << std::endl;
        std::cout << "\tHemisphericality: " << co.value()["attributes"]["Hemisphericality"] << std::endl;
        std::cout << "\tRoughnessIndex: " << co.value()["attributes"]["RoughnessIndex"] << std::endl;
        std::cout << "\tOrientations: " << co.value()["geometry"][0]["semantics"]["values"];
    }

    return 0;
}