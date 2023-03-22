#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "Fracture3D.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "GeometryUtilities.hpp"

namespace XFEM_3D
{
namespace Utilities
{

    double signedDistanceFunction(Eigen::Vector3d point, Fracture3D& fracture);

    int heaviside(double x);

    std::vector<Gedim::GeometryUtilities::Polyhedron> splitTetrahedronInSubTetrahedra(Gedim::GeometryUtilities::Polyhedron elementAsPolyhedron,
                                                                                       Fracture3D *fracture,
                                                                                       Gedim::GeometryUtilities *geometryUtilities,
                                                                                       Gedim::MeshUtilities *meshUtilities);


}
}






#endif // UTILITIES_HPP
