#include "testAllZero.hpp"


namespace XFEM_3D
{

namespace testAllZero
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

        exactSolution(ii_std) = 0;
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
    return 0;
}


double g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                   unsigned int neumannBorderMarker,
                                   Gedim::MeshMatricesDAO blockMesh,
                                   Fracture3D* fracture)
{
    return 0;
}

setBoundaryConditionsResultType setBoundaryConditions(Gedim::MeshMatricesDAO& mesh)
{

    Eigen::MatrixXi pivot(mesh.Cell0DTotalNumber(), 2);
    std::vector<unsigned int> Neumann_triangles;
    unsigned int num_Dirichlet_3D = 0;
    unsigned int numDOF_3D_std = 0;
    for (unsigned int glob_id_point = 0; glob_id_point < mesh.Cell0DTotalNumber(); glob_id_point++)
    {
        Eigen::Vector3d point = mesh.Cell0DCoordinates(glob_id_point);
        double x = point.x(),
               y = point.y(),
               z = point.z();

        bool point_is_dirichlet_marker_1 = (z == 0.0), point_is_dirichlet_marker_3 = (z == 1.0);
        bool point_is_dirichlet = point_is_dirichlet_marker_1 || point_is_dirichlet_marker_3;

        if (point_is_dirichlet)
        {
            mesh.Cell0DSetMarker(glob_id_point, point_is_dirichlet_marker_1 * 1 + point_is_dirichlet_marker_3 * 3);
            num_Dirichlet_3D++;
            pivot(glob_id_point, 0) = -num_Dirichlet_3D;
            pivot(glob_id_point, 1) = -1;
        }
        else
        {
            mesh.Cell0DSetMarker(glob_id_point, 0);
            numDOF_3D_std++;
            pivot(glob_id_point, 0) = numDOF_3D_std;
            pivot(glob_id_point, 1) = -1;
        }
    }


    // Costruzione della struttura dati per i triangoli di Neumann
    for (unsigned int glob_id_triangle; glob_id_triangle < mesh.Cell2DTotalNumber(); glob_id_triangle++)
    {
        mesh.Cell2DSetMarker(glob_id_triangle, 0);

        std::vector<unsigned int> vertices_global_IDs = mesh.Cell2DVertices(glob_id_triangle);

        Eigen::Vector3d point1 = mesh.Cell0DCoordinates(vertices_global_IDs.at(0)),
                        point2 = mesh.Cell0DCoordinates(vertices_global_IDs.at(1)),
                        point3 = mesh.Cell0DCoordinates(vertices_global_IDs.at(2));

        double x1 = point1.x(), x2 = point2.x(), x3 = point3.x(),
               y1 = point1.y(), y2 = point2.y(), y3 = point3.y();

        bool triangle_lies_on_cube_face_with_marker_2  = (x1 == 0.0 && x2 == 0.0 && x3 == 0.0);
        bool triangle_lies_on_cube_face_with_marker_4  = (x1 == 1.0 && x2 == 1.0 && x3 == 1.0);
        bool triangle_lies_on_cube_face_with_marker_6  = (y1 == 0.0 && y2 == 0.0 && y3 == 0.0);
        bool triangle_lies_on_cube_face_with_marker_8  = (y1 == 1.0 && y2 == 1.0 && y3 == 1.0);


        // Faccia di Neumann con marker 2: {x=0}
        if (triangle_lies_on_cube_face_with_marker_2)
        {
            mesh.Cell0DSetMarker(vertices_global_IDs.at(0), 2);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(1), 2);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(2), 2);

            mesh.Cell2DSetMarker(glob_id_triangle, 2);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 4: {x=1}
        if (triangle_lies_on_cube_face_with_marker_4)
        {
            mesh.Cell0DSetMarker(vertices_global_IDs.at(0), 4);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(1), 4);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(2), 4);

            mesh.Cell2DSetMarker(glob_id_triangle, 4);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 6: {y=0}
        if (triangle_lies_on_cube_face_with_marker_6)
        {
            mesh.Cell0DSetMarker(vertices_global_IDs.at(0), 6);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(1), 6);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(2), 6);

            mesh.Cell2DSetMarker(glob_id_triangle, 6);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 8: {y=1}
        if (triangle_lies_on_cube_face_with_marker_8)
        {
            mesh.Cell0DSetMarker(vertices_global_IDs.at(0), 8);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(1), 8);
            mesh.Cell0DSetMarker(vertices_global_IDs.at(2), 8);

            mesh.Cell2DSetMarker(glob_id_triangle, 8);
            Neumann_triangles.push_back(glob_id_triangle);
        }
    }

    setBoundaryConditionsResultType result;
    result.pivot = pivot;
    result.Neumann_triangles = Neumann_triangles;

    return result;
}


}
}
