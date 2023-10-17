#include "testOneZero.hpp"


namespace XFEM_3D
{
namespace testOneZero
{


Eigen::VectorXd exactSolution(Eigen::MatrixXd meshPointsCoordinates,
                                              unsigned int numDofs3D_std,
                                              Eigen::MatrixXi pivot,
                                              Fracture3D* fracture)
{
    Eigen::VectorXd exactSolution(numDofs3D_std);

    for (unsigned int glob_id_node = 0; glob_id_node < meshPointsCoordinates.cols(); glob_id_node++)
    {
        Eigen::Vector3d point = meshPointsCoordinates.col(glob_id_node);

        int ii_std = pivot(glob_id_node, 0);

        if (ii_std < 0)
            continue;

        // Solita convenzione: i DOF sono numerati a partire da 1, ma C++ numera nelle matrici a partire da 0.
        ii_std--;

        bool point_on_positive_side_of_fracture = Utilities::signedDistanceFunction(point, *fracture);
        if (point_on_positive_side_of_fracture)
            exactSolution(ii_std) = 1.0;
        else
            exactSolution(ii_std) = 0.0;
    }

    return exactSolution;

}

Eigen::Vector3d exactSolutionGradient(Eigen::Vector3d point, Fracture3D* fracture)
{

    Eigen::Vector3d gradient;

    gradient << 0, 0, 0;

    return gradient;
}


double forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture)
{

    bool point_on_positive_side_of_fracture = Utilities::signedDistanceFunction(pointCoords, fracture);
    if (point_on_positive_side_of_fracture)
        return 0.0;
    else
        return 0.0;
}


double g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                   unsigned int neumannBorderMarker,
                                   Gedim::MeshMatricesDAO blockMesh,
                                   Fracture3D* fracture)
{
    return 0;
}

}
}
