#ifndef __geomtools__
#define __geomtools__

#include "definitions.h"


Plane                 get_best_fitted_plane(const std::vector<Point3> &lspts);
void                  mark_domains(CT& ct);
void                  mark_domains(CT& ct, CT::Face_handle start, int index, std::list<CT::Edge>& border);

std::vector<std::vector<int>>
construct_ct_one_face(const std::vector<std::vector<int>>& lsRings, 
                      const std::vector<Point3>& lspts);

#endif 