#include "discontinousTestProblem_1.h"


namespace XFEM_3D
{
namespace DiscontinousTestProblem_1
{


Eigen::VectorXd exactSolution(Eigen::MatrixXd meshPointsCoordinates,
                                              unsigned int numDofs3D_std,
                                              Eigen::MatrixXi pivot,
                                              Fracture3D* fracture)
{
    Eigen::VectorXd exactSolution(numDofs3D_std);
    Eigen::Vector3d point;

    for (unsigned int glob_id_node = 0; glob_id_node < meshPointsCoordinates.cols(); glob_id_node++)
    {
        int ii_std = pivot(glob_id_node, 0);

        if (ii_std < 0)
            continue;

        // Solita convenzione: i DOF sono numerati a partire da 1, ma C++ numera nelle matrici a partire da 0.
        ii_std--;

        point = meshPointsCoordinates.col(ii_std);
        double x = point.x(),
               y = point.y(),
               z = point.z();

        if (Utilities::signedDistanceFunction(point, *fracture) >= 0)
            exactSolution(ii_std) = -16 * x * y * z * (1-x) * (1-y) * (1-z);
        else
            exactSolution(ii_std) =  16 * x * y * z * (1-x) * (1-y) * (1-z);
    }

    return exactSolution;

}

Eigen::Vector3d exactSolutionGradient(Eigen::Vector3d point, Fracture3D* fracture)
{
    double x = point.x(),
           y = point.y(),
           z = point.z();

    Eigen::Vector3d gradient;

    if (Utilities::signedDistanceFunction(point, *fracture) >= 0)
    {
        gradient << -16 * y * z * (1-y) * (1-z) * (1-2*x),
                    -16 * x * z * (1-x) * (1-z) * (1-2*y),
                    -16 * x * y * (1-x) * (1-y) * (1-2*z);
    }
    else
    {
        gradient <<  16 * y * z * (1-y) * (1-z) * (1-2*x),
                     16 * x * z * (1-x) * (1-z) * (1-2*y),
                     16 * x * y * (1-x) * (1-y) * (1-2*z);
    }

    return gradient;
}


double forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture)
{
    double x = pointCoords(0), y = pointCoords(1), z = pointCoords(2);

    if (Utilities::signedDistanceFunction(pointCoords, fracture) >= 0)
       return -32*( x*y*(1-x)*(1-y) + x*z*(1-x)*(1-z) + y*z*(1-y)*(1-z) );

    else
        return  32*( x*y*(1-x)*(1-y) + x*z*(1-x)*(1-z) + y*z*(1-y)*(1-z) );
}


double g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                   unsigned int neumannBorderMarker,
                                   Gedim::MeshMatricesDAO blockMesh,
                                   Fracture3D* fracture)
{

    Eigen::Vector3d point1 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(0)),
                    point2 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(1)),
                    point3 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(2));

    Eigen::Vector3d solutionGradient_atPoint1 = exactSolutionGradient(point1, fracture),
                    solutionGradient_atPoint2 = exactSolutionGradient(point2, fracture),
                    solutionGradient_atPoint3 = exactSolutionGradient(point3, fracture);

    double averaged_gN, gN1, gN2, gN3;
    Eigen::Vector3d externalCubeFaceNormal;

    switch (neumannBorderMarker) {
    case 2: // {x=0}
        externalCubeFaceNormal << -1, 0, 0;
        break;
    case 4: // {x=1}
        externalCubeFaceNormal << 1, 0, 0;
        break;
    case 6: // {y=0}
        externalCubeFaceNormal << 0, -1, 0;
        break;
    case 8: // {y=1}
        externalCubeFaceNormal << 0, 1, 0;
        break;
    default:
        break;
    }

    gN1 = solutionGradient_atPoint1.transpose() * externalCubeFaceNormal;
    gN2 = solutionGradient_atPoint2.transpose() * externalCubeFaceNormal;
    gN3 = solutionGradient_atPoint3.transpose() * externalCubeFaceNormal;

    averaged_gN = (1 / 3) * (gN1 + gN2 + gN3);

    return averaged_gN;

}

}
}
