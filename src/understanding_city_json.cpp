//
// Created by Sharath Chandra on 18-Mar-23.
//

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <CGAL/Random.h>

#include "definitions.h"
#include "geomtools.h"

#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

//int main() {
//
//    // Reading the obj file
//    const char* filename = "../data/myfile.city.json";
//    std::cout << "Processing: " << filename << std::endl;
//    std::ifstream input_stream;
//    std::ifstream input(filename);
//    json j;
//    input >> j; //-- store the content of the file in a nlohmann::json object
//    input.close();
//
//    int i {1};
//    for (auto& co : j["CityObjects"].items()) {
//        std::cout << "City Object - " << i << ": " << co.key() << std::endl;
//        for (auto& obj_info :  co.value().items()) {
//            std::cout << "\tkey: " << obj_info.key() << std::endl;
//            std::cout << "\tvalue: " << obj_info.value() << std::endl;
//            if (obj_info.key() == "geometry") {
//                for (auto& geom_dict : obj_info.value()[0].items()) {
//                    std::cout << "\t\tgeom key - " << geom_dict.key() << std::endl;
//                    if (geom_dict.key() == "boundaries") {
//                        std::cout << "\t\t\tnum faces: " << geom_dict.value()[0].size() << std::endl;
//                    }
//                    if ( geom_dict.key() == "semantics") {
//                        std::cout << "\t\t\t" << geom_dict.value() << std::endl;
//                        for (auto& semantics_dict: geom_dict.value().items()) {
//                            std::cout << "\t\t\t" << semantics_dict.key() << std::endl;
//                        }
//                    }
//                }
//            }
//        }
//
//        i++;
////        if (i == 100)
////            break;
//    }
//
//    return 0;
//}