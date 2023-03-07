#include "Utilities.hpp"

namespace XFEM_3D
{
namespace Utilities
{


    double signedDistanceFunction(Eigen::Vector3d point, Fracture3D& fracture)
    {
        Eigen::Vector3d normal = fracture.getNormal();
        Eigen::Vector3d x = fracture.getOrigin();
        double a = normal(0);
        double b = normal(1);
        double c = normal(2);
        double d = - a*x(0) - b*x(1) - c*x(2);
        double res;

        res = (a*point(0) + b*point(1) + c*point(2) + d) / (sqrt(a*a + b*b + c*c));

        return res;

    }

    int heaviside (double input)
    {
        if (input >= 0)
            return 1;
        else
            return 0;
    }


    Eigen::Vector4d lagrangeBasisCoeff3D(Eigen::Vector3d nodeCoordinates, Gedim::GeometryUtilities::Polyhedron element)
    {

        Eigen::Vector4d coeff;

        double x_i = nodeCoordinates(0);
        double y_i = nodeCoordinates(1);
        double z_i = nodeCoordinates(2);

        Eigen::MatrixXd elementPoints = element.Vertices;
        Eigen::MatrixXd elementPointsExceptNode(3,3);

        // Determine which point in elementPoints is the node (has to be one of them). You are goign to exclude it
        // from the conditions of the linear system.

        unsigned int chosen_k = 42;
        unsigned int position = 0;
        for (unsigned int k = 0; k < 4; k++)
        {
            Eigen::Vector3d diff = elementPoints.col(k) - nodeCoordinates;
            double diffNorm = diff.norm();

            if (diffNorm < 1e-16)
                chosen_k = k;
            else
            {
                elementPointsExceptNode.col(position) = elementPoints.col(k);
                position++;
            }
        }

        if (chosen_k == 42)
        {
            std::cerr << "Something went wrong trying to determine the coefficient of the 3D lagrange basis functions!" << std::endl;
        }

        //

        Eigen::Vector3d x_e, y_e, z_e;
        for (unsigned int k = 0; k < 3; k++)
        {
            x_e(k) = elementPointsExceptNode.col(k).x();
            y_e(k) = elementPointsExceptNode.col(k).y();
            z_e(k) = elementPointsExceptNode.col(k).z();

        }

        // Construction of the 4x4 linear system matrix and right hand side.

        Eigen::Matrix4d A;
        Eigen::Vector4d b;

        A << x_i,    y_i,    z_i,    1,
                x_e(0), y_e(0), z_e(0), 1,
                x_e(1), y_e(1), z_e(1), 1,
                x_e(2), y_e(2), z_e(2), 1;

        b << 1, 0, 0, 0;

        // Solution of the 4x4 system

        coeff = A.colPivHouseholderQr().solve(b);

        return coeff;
    }



    Eigen::Vector3d lagrangeBasisCoeff2D(Eigen::Vector3d nodeCoordinates, Eigen::MatrixXd elementPoints)
    {

        Eigen::Vector3d coeff;

        double x_i = nodeCoordinates(0);
        double y_i = nodeCoordinates(1);

        Eigen::MatrixXd elementPointsExceptNode;

        // Determine which point in elementPoints is the node (has to be one of them). You are goign to exclude it
        // from the conditions of the linear system.

        unsigned int chosen_k = 42;
        unsigned int position = 0;
        for (unsigned int k = 0; k < 3; k++)
        {
            Eigen::Vector3d diff = elementPoints.col(k) - nodeCoordinates;

            if (diff.norm() < 1e-16)
                chosen_k = k;
            else
            {
                elementPointsExceptNode.col(position) = elementPoints.col(k);
                position++;
            }
        }

        if (chosen_k == 42)
        {
            std::cerr << "Something went wrong trying to determine the coefficient of the 2D lagrange basis functions!" << std::endl;
        }

        Eigen::Vector2d x_e, y_e;
        for (unsigned int k = 0; k < 2; k++)
        {
            x_e(k) = elementPointsExceptNode.col(k).x();
            y_e(k) = elementPointsExceptNode.col(k).y();
        }

        // Construction of the 4x4 linear system matrix and right hand side.

        Eigen::Matrix3d A;
        Eigen::Vector3d b;

        A << x_i,    y_i,    1,
                x_e(0), y_e(0), 1,
                x_e(1), y_e(1), 1,
                x_e(2), y_e(2), 1;

        b << 1, 0, 0;

        // Solution of the 3x3 system

        coeff = A.colPivHouseholderQr().solve(b);

        return coeff;

    }


    double evaluate3DLagrange(Eigen::Vector3d x, Eigen::Vector4d coeff)
    {
        return coeff(0)*x(0) + coeff(1)*x(1) + coeff(2)*x(2) + coeff(3);
    }

    Eigen::Vector3d evaluate3DLagrangeGrad(Eigen::Vector3d x, Eigen::Vector4d coeff)
    {
        Eigen::Vector3d grad;
        grad << coeff(0), coeff(1), coeff(2);
        return grad;
    }

    double evaluate2DLagrange(Eigen::Vector3d x, Eigen::Vector3d coeff)
    {
        return coeff(0)*x(0) + coeff(1)*x(1) + coeff(2);
    }

    Eigen::Vector2d evaluate2DLagrangeGrad(Eigen::Vector3d x, Eigen::Vector3d coeff)
    {
        Eigen::Vector2d grad;
        grad << coeff(0), coeff(1);
        return grad;
    }





}

}

