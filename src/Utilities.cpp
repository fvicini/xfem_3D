#include "Utilities.hpp"
#include<tuple>

namespace XFEM_3D
{
namespace Utilities
{


    double signedDistanceFunction(Eigen::Vector3d point, Fracture3D& fracture)
    {
        Eigen::Vector3d normal = fracture.getNormal();
        Eigen::Vector3d x = fracture.getOrigin();
        double a = normal(0),
               b = normal(1),
               c = normal(2);
        double d = - a*x(0) - b*x(1) - c*x(2);
        double res;

        res = (a*point(0) + b*point(1) + c*point(2) + d) / (sqrt(a*a + b*b + c*c));

        return res;
    }



    int heaviside (double input)
    {
        if (input > 0)
            return 1;
        else
            return -1;
    }



    std::vector<Gedim::GeometryUtilities::Polyhedron> splitTetrahedronInSubTetrahedra(Gedim::GeometryUtilities::Polyhedron element,
                                                                                       Fracture3D *fracture,
                                                                                       Gedim::GeometryUtilities *geometryUtilities,
                                                                                       Gedim::MeshUtilities *meshUtilities)
    {
        // Variabili geometriche ausiliarie -------------------------------------------------------------------------------------------------
        const std::vector<std::vector<bool>> polyhedronFaceEdgeDirections = geometryUtilities->PolyhedronFaceEdgeDirections(element.Vertices,
                                                                                                                            element.Edges,
                                                                                                                            element.Faces);
        const std::vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtilities->PolyhedronFaceVertices(element.Vertices,
                                                                                                              element.Faces);
        const Eigen::MatrixXd polyhedronEdgeTangents = geometryUtilities->PolyhedronEdgeTangents(element.Vertices, element.Edges);
        const std::vector<Eigen::MatrixXd> polyhedronFaceTangents = geometryUtilities->PolyhedronFaceEdgeTangents(element.Vertices,
                                                                                                                  element.Edges,
                                                                                                                  element.Faces,
                                                                                                                  polyhedronFaceEdgeDirections,
                                                                                                                  polyhedronEdgeTangents);
        const std::vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities->PolyhedronFaceTranslations(polyhedronFaceVertices);
        const std::vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities->PolyhedronFaceNormals(polyhedronFaceVertices);
        const std::vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities->PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                              polyhedronFaceNormals,
                                                                                                                              polyhedronFaceTranslations);
        // ------------------------------------------------------------------------------------------------------------------------------------

        // Suddivido il tetraedro della mesh in due poliedri
        Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult result = geometryUtilities->SplitPolyhedronWithPlane(element.Vertices,
                                                                                                                      element.Edges,
                                                                                                                      element.Faces,
                                                                                                                      polyhedronFaceVertices,
                                                                                                                      polyhedronFaceTangents,
                                                                                                                      polyhedronFaceTranslations,
                                                                                                                      polyhedronFaceRotationMatrices,
                                                                                                                      fracture->getNormal(),
                                                                                                                      fracture->getOrigin(),
                                                                                                                      fracture->getRotation(),
                                                                                                                      fracture->getTranslation());

        std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra;
        std::vector<Gedim::GeometryUtilities::Polyhedron> subPolyhedra = geometryUtilities->SplitPolyhedronWithPlaneResultToPolyhedra(result);


        for (unsigned int p = 0; p < subPolyhedra.size(); p++)
        {
            bool currentPolyhedronIsAlreadyTetrahedron = subPolyhedra.at(p).Vertices.cols() == 4;

            if (currentPolyhedronIsAlreadyTetrahedron)
                subTetrahedra.push_back(subPolyhedra.at(p));

            else
            {
                Gedim::MeshMatrices aux_meshData;
                Gedim::MeshMatricesDAO aux_mesh(aux_meshData);

                meshUtilities->CreateTetrahedralMesh(subPolyhedra.at(p).Vertices,
                                                     subPolyhedra.at(p).Edges,
                                                     subPolyhedra.at(p).Faces,
                                                     100000,
                                                     aux_mesh,
                                                     "Qpfezna");

                for (unsigned short int t = 0; t < aux_mesh.Cell3DTotalNumber(); t++)
                    subTetrahedra.push_back(meshUtilities->MeshCell3DToPolyhedron(aux_mesh, t));
            }
        }

        return subTetrahedra;
    }



    Eigen::Matrix<double, 8, 8> compute_k_enr(Eigen::Vector3d x_st,
                                              Gedim::GeometryUtilities::Polyhedron element,
                                              Fracture3D* fracture)
    {
        Eigen::Matrix<double, 8, 8> k_enr;
        k_enr  << 0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0;

        for (unsigned int i = 0; i < 8; i++)
        {
            for (unsigned int j = 0; j < 8; j++)
            {
                if (i < 4 && j < 4)
                {
                    k_enr(i,j) = 1;
                }

                if (i < 4 && j > 4)
                {
                    int j_tilde = j % 4;
                    Eigen::Vector3d x_j_tilde = element.Vertices.col(j_tilde);

                    k_enr(i,j) = Utilities::heaviside(Utilities::signedDistanceFunction(x_st, *fracture))
                                    - Utilities::heaviside(Utilities::signedDistanceFunction(x_j_tilde, *fracture));
                }

                if (i > 4 && j < 4)
                {
                    int i_tilde = i % 4;
                    Eigen::Vector3d x_i_tilde = element.Vertices.col(i_tilde);

                    k_enr(i,j) = Utilities::heaviside(Utilities::signedDistanceFunction(x_st, *fracture))
                                    - Utilities::heaviside(Utilities::signedDistanceFunction(x_i_tilde, *fracture));
                }

                if (i > 4 && j > 4)
                {
                    int i_tilde = i % 4, j_tilde = j % 4;

                    Eigen::Vector3d x_i_tilde = element.Vertices.col(i_tilde),
                                    x_j_tilde = element.Vertices.col(j_tilde);

                    k_enr(i,j) = (Utilities::heaviside(Utilities::signedDistanceFunction(x_st, *fracture))
                                     - Utilities::heaviside(Utilities::signedDistanceFunction(x_i_tilde, *fracture)))
                               * (Utilities::heaviside(Utilities::signedDistanceFunction(x_st, *fracture))
                                  - Utilities::heaviside(Utilities::signedDistanceFunction(x_j_tilde, *fracture)));
                }
            }
        }

        return k_enr;
    }
}

}

