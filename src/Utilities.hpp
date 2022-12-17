#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "Fracture3D.hpp"

namespace XFEM_3D
{
namespace Utilities
{

    double signedDistanceFunction(Eigen::Vector3d point, Fracture3D& fracture);

    int heaviside(double x);

    /// This Function computes the coefficient of the restriction of the basis lagrange function relative to node
    /// with nodeCoordinates on the Polyhedron element. The coefficients are A,B,C,D from the following expression
    ///          phi(x,y,z) = Ax + By + Cz + D
    /// The coeff vector returned by this function is to be given as second input to evaluateLagrange, evaluateLagrangeGrad for evaluating the correct basis
    /// function when integrating.
    Eigen::Vector4d lagrangeBasisCoeff(Eigen::Vector3d nodeCoordinates, Gedim::GeometryUtilities::Polyhedron element);

    double evaluateLagrange(Eigen::Vector3d x, Eigen::Vector4d coeff);

    Eigen::Vector3d evaluateLagrangeGrad(Eigen::Vector3d x, Eigen::Vector4d coeff);


}
}






#endif // UTILITIES_HPP
