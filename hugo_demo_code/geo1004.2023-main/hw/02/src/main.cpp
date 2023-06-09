/*
  geo1004.2023
  hw02 help code
  Hugo Ledoux <h.ledoux@tudelft.nl>
  2023-03-01
*/

#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "definitions.h"
#include "geomtools.h"

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder
using json = nlohmann::json;

std::vector<Point3> get_coordinates(const json& j, bool translate = true);
void                save2obj(std::string filename, const json& j);
void                enrich_and_save(std::string filename, json& j);



int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or 2b.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../data/2b.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();
  
  //-- convert each City Object in the file to OBJ and save to a file
  save2obj("out.obj", j);

  //-- enrich with some attributes and save to a new CityJSON 
  enrich_and_save("out.city.json", j);

  return 0;
}


//-- write the OBJ file
void save2obj(std::string filename, const json& j) {
  std::ofstream ofile(filename);
  //-- fetch all the vertices in real-world coordinates (so "transform" is applied)
  std::vector<Point3> lspts = get_coordinates(j, true);
  for (auto& p : lspts) {
    ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
  }
  //-- iterate over each object in the file and output the CDT
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if ( (g["type"] == "Solid") && (g["lod"] == "2.2") ) {   //-- LoD2.2 only!!!!!
        ofile << "o " << co.key() << std::endl;
        for (int i = 0; i < g["boundaries"].size(); i++) {
          for (int j = 0; j < g["boundaries"][i].size(); j++) {
            std::vector<std::vector<int>> gb = g["boundaries"][i][j];
            std::vector<std::vector<int>> trs = construct_ct_one_face(gb, lspts);
            for (auto& tr : trs) {
              ofile << "f " << (tr[0] + 1) << " " << (tr[1] + 1) << " " << (tr[2] + 1) << std::endl;
            }
          }
        }
      }
    }
  }
  ofile.close();
  std::cout << "OBJ file written to disk: " << filename << std::endl;
}

//-- add a new attribute "volume" to each City Object and assign a random value
void enrich_and_save(std::string filename, json& j) {
  //-- seed to generate a random number 
  //-- https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  std::random_device rd;  
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<> dis(100.0, 200.0);
  //-- add an attribute "volume"
  for (auto& co : j["CityObjects"]) {
    if (co["type"] == "Building") {
      co["attributes"]["volume"] = dis(gen);
    }
  }
  //-- write to disk the modified city model (myfile.city.json)
  std::ofstream o(filename);
  o << j.dump(2) << std::endl;
  o.close();
  std::cout << "Enriched CityJSON file written to disk: " << filename << std::endl;
}


//-- get real-world coordinates from the "vertices" property
//-- https://www.cityjson.org/specs/#transform-object
//-- param translate is to use the translation in the "transform",
//-- it can be put to false to make the coords smaller (and better for computations)
std::vector<Point3> get_coordinates(const json& j, bool translate) {
  std::vector<Point3> lspts;
  std::vector<std::vector<int>> lvertices = j["vertices"];
  if (translate) {
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
      double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
      double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
      lspts.push_back(Point3(x, y, z));
    } 
  } else {
    //-- do not translate, useful to keep the values low for downstream processing of data
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>());
      double y = (vi[1] * j["transform"]["scale"][1].get<double>());
      double z = (vi[2] * j["transform"]["scale"][2].get<double>());
      lspts.push_back(Point3(x, y, z));
    }
  }
  return lspts;
}
