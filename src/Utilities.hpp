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
    Eigen::Vector4d lagrangeBasisCoeff3D(Eigen::Vector3d nodeCoordinates, Gedim::GeometryUtilities::Polyhedron element);

    Eigen::Vector3d lagrangeBasisCoeff2D(Eigen::Vector3d nodeCoordinates, Eigen::MatrixXd elementPoints);

    double evaluate3DLagrange(Eigen::Vector3d x, Eigen::Vector4d coeff);

    Eigen::Vector3d evaluate3DLagrangeGrad(Eigen::Vector3d x, Eigen::Vector4d coeff);

    double evaluate2DLagrange(Eigen::Vector3d x, Eigen::Vector3d coeff);

    Eigen::Vector2d evaluate2DLagrangeGrad(Eigen::Vector3d x, Eigen::Vector3d coeff);


}
}






#endif // UTILITIES_HPP
