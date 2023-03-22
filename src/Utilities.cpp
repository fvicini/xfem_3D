#include "Utilities.hpp"

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
        if (input >= 0)
            return 1;
        else
            return 0;
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
}

}

