#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "definitions.h"
#include "geomtools.h"

#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

std::vector<Point3> get_coordinates(const json& j, bool translate = true);
void                save2obj(std::string filename, const json& j);
void                enrich_and_save(std::string filename, json& j);


//int main(int argc, const char * argv[]) {
//    std::cout << "hello world from test! " << std::endl;
//    //-- will read the file passed as argument or cube.city.json if nothing is passed
//    const char* filename = (argc > 1) ? argv[1] : "../data/cube.city.json";
//    std::cout << "Processing: " << filename << std::endl;
//    std::ifstream input(filename);
//    json j;
//    input >> j; //-- store the content of the file in a nlohmann::json object
////    std::cout << "json: " << j << std::endl;
//    input.close();
//
//    std::vector<std::vector<int>>  lvertices = j["vertices"];
//    int i {0};
//    for (const auto ver: lvertices) {
//        std::cout << "\nV " << i << ": ";
//        for (const auto coor: ver) {
//            std::cout << coor << " ";
//        }
//        i++;
//    }
//
////    std::cout << "\nCityObjects: " << j["CityObjects"] << std::endl;
////    std::cout << "\nCityObjects items: " << j["CityObjects"].items() << std::endl;
//
////    std::cout << "CityObjects items: " << j["CityObjects"]["id-1"] << std::endl;
////    std::cout << "CityObjects id-1 items: " << j["CityObjects"]["id-1"].items() << std::endl;
////    std::cout << "CityObjects id-1 items begin: " << j["CityObjects"]["id-1"].items().begin() << std::endl;
//
////    std::cout << "CityObjects geometry: " << j["CityObjects"]["id-1"]["geometry"] << std::endl;
////    std::cout << "CityObjects geometry: " << j["CityObjects"]["id-1"]["geometry"]["boundaruies"] << std::endl;
//    std::vector<Point3> lspts = get_coordinates(j, true);
////
//    for (auto& co: j["CityObjects"].items()) {
////        std::cout << "\nIm inside cityobjects heheheheheheheheheh" << std::endl;
////        std::cout << co << std::endl;
////        std::cout << co.key() << std::endl;
////        std::cout << co.value() << std::endl;
////        std::cout << co.value()["geometry"] << std::endl;
////        std::cout << co.value()["geometry"][0] << std::endl;
////        std::cout << co.value()["geometry"][0]["lod"] << std::endl;
////        std::cout << "=======================" << std::endl;
//        std::cout << std::endl;
//        for (auto& co: j["CityObjects"].items()) {
//            for (auto& g: co.value()["geometry"]) {
//                std::cout << g << std::endl;
//                if (g["type"] == "Solid") {
//                    if (g["lod"] == "1") {
////                        std::cout << "it is a solid and has lod 1" << std::endl;
//                        std::cout << "\no " << co.key() << std::endl;
//                        for (int i = 0; i < g["boundaries"].size(); i++) {
//                            for (int j = 0; j < g["boundaries"][i].size(); j++) {
//                                std::vector<std::vector<int>> gb = g["boundaries"][i][j];
//                                std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
//                                for (auto& tr: trs) {
//                                    std::cout << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
//                                }
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//
//
////        for (auto& g: co.value().items()) {
////            std::cout << "megh" << std::endl;
////            std::cout << g << std::endl;
////        }
////        for (auto& g: co.value().items()) {
////            std::cout << "megh" << std::endl;
////            std::cout << g << std::endl;
////        }
//    }
//
//    //-- convert each City Object in the file to OBJ and save to a file
////    save2obj("../data/out.obj", j);
//
//
//    return 0;
//}

//
//void save2obj(std::string filename, const json& j) {
//    std::ofstream ofile(filename);
//    //-- fetch all the vertices in real-world coordinates (so "transform" is applied)
//    std::vector<Point3> lspts = get_coordinates(j, true);
//    for (auto& p : lspts) {
//        ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//    }
//    //-- iterate over each object in the file and output the CDT
//    for (auto& co : j["CityObjects"].items()) {
//        for (auto& g : co.value()["geometry"]) {
//            if ( (g["type"] == "Solid") && (g["lod"] == "2.2") ) {   //-- LoD2.2 only!!!!!
//                ofile << "o " << co.key() << std::endl;
//                for (int i = 0; i < g["boundaries"].size(); i++) {
//                    for (int j = 0; j < g["boundaries"][i].size(); j++) {
//                        std::vector<std::vector<int>> gb = g["boundaries"][i][j];
//                        std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
//                        for (auto& tr : trs) {
//                            ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    ofile.close();
//    std::cout << "OBJ file written to disk: " << filename << std::endl;
//}
//
//std::vector<Point3> get_coordinates(const json& j, bool translate) {
//    std::vector<Point3> lspts;
//    std::vector<std::vector<int>> lvertices = j["vertices"];
//    if (translate) {
//        for (auto& vi : lvertices) {
//            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
//            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
//            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
//            lspts.push_back(Point3(x, y, z));
//        }
//    } else {
//        //-- do not translate, useful to keep the values low for downstream processing of data
//        for (auto& vi : lvertices) {
//            double x = (vi[0] * j["transform"]["scale"][0].get<double>());
//            double y = (vi[1] * j["transform"]["scale"][1].get<double>());
//            double z = (vi[2] * j["transform"]["scale"][2].get<double>());
//            lspts.push_back(Point3(x, y, z));
//        }
//    }
//    return lspts;
//}
