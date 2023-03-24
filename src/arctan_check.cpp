////
//// Created by Sharath Chandra on 24-Mar-23.
////
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <random>
//#include <cmath>
//#include "definitions.h"
//#include "geomtools.h"
//
//double roof_elevation_azimuth(K::Vector_3& normal) {
//    double factor = std::pow(10.0, 1);
//    double x = normal.x(), y = normal.y(), z = normal.z();
//
//    double elevation_deg = atan2(z, sqrt(x*x + y*y)) * 180 / M_PI;
//    double elevation_round = std::round (elevation_deg * factor) / factor;
//
//    double azimuth;
//    if (x > 0) {
//        azimuth = atan2(x, y) * 180 / M_PI;
//    } else if (x < 0) {
//        azimuth = (atan2(x, y) + 2 * M_PI) * 180 / M_PI;
//    } else if (y > 0) {
//        azimuth = 0.0;
//    } else if (y < 0) {
//        azimuth = 180.0;
//    } else {
//        // vector is at the origin
//        azimuth = 0.0;
//    }
//
//    return azimuth;
//}
//
//int main() {
//    K::Vector_3 normal (0.418156, -0.128336, 0.899264);
//    std::cout << "azimuth: " << roof_elevation_azimuth(normal) << std::endl;
//    return 0;
//}