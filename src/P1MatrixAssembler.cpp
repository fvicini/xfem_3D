#include "P1MatrixAssembler.hpp"
#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "MapTetrahedron.hpp"
#include "Utilities.hpp"

namespace XFEM_3D
{


// SETTERS ********************************************************************************************************

void P1MatrixAssembler::setHD_Mesh(Gedim::MeshMatricesDAO *newHD_Mesh)
{
    hD_Mesh = newHD_Mesh;
}

void P1MatrixAssembler::setHD_Pivot(Eigen::VectorXi *newHD_Pivot)
{
    hD_Pivot = newHD_Pivot;
}

void P1MatrixAssembler::setToEnrich_nodes(Eigen::SparseMatrix<unsigned int> *newToEnrich_nodes)
{
    toEnrich_nodes = newToEnrich_nodes;
}

void P1MatrixAssembler::setHF_Mesh(Gedim::MeshMatricesDAO *newHF_Mesh)
{
    hF_Mesh = newHF_Mesh;
}

void P1MatrixAssembler::setHF_Pivot(Eigen::VectorXi *newHF_Pivot)
{
    hF_Pivot = newHF_Pivot;
}

void P1MatrixAssembler::setPsiP_Mesh(Gedim::MeshMatricesDAO *newPsiP_Mesh)
{
    psiP_Mesh = newPsiP_Mesh;
}

void P1MatrixAssembler::setPsiP_Pivot(Eigen::VectorXi *newPsiP_Pivot)
{
    psiP_Pivot = newPsiP_Pivot;
}

void P1MatrixAssembler::setPsiM_Mesh(Gedim::MeshMatricesDAO *newPsiM_Mesh)
{
    psiM_Mesh = newPsiM_Mesh;
}

void P1MatrixAssembler::setPsiM_Pivot(Eigen::VectorXi *newPsiM_Pivot)
{
    psiM_Pivot = newPsiM_Pivot;
}

void P1MatrixAssembler::setPsiF_Mesh(Gedim::MeshMatricesDAO *newPsiF_Mesh)
{
    psiF_Mesh = newPsiF_Mesh;
}

void P1MatrixAssembler::setPsiF_Pivot(Eigen::VectorXi *newPsiF_Pivot)
{
    psiF_Pivot = newPsiF_Pivot;
}

void P1MatrixAssembler::setLambdaD_Mesh(Gedim::MeshMatricesDAO *newLambdaD_Mesh)
{
    lambdaD_Mesh = newLambdaD_Mesh;
}

void P1MatrixAssembler::setLambdaD_Pivot(Eigen::VectorXi *newLambdaD_Pivot)
{
    lambdaD_Pivot = newLambdaD_Pivot;
}

void P1MatrixAssembler::setLambdaF_Mesh(Gedim::MeshMatricesDAO *newLambdaF_Mesh)
{
    lambdaF_Mesh = newLambdaF_Mesh;
}

void P1MatrixAssembler::setLambdaF_Pivot(Eigen::VectorXi *newLambdaF_Pivot)
{
    lambdaF_Pivot = newLambdaF_Pivot;
}


// *****************************************************************************************************

void P1MatrixAssembler::assemble_AhD(Gedim::Eigen_SparseArray<>& AhD,
                       Gedim::Eigen_SparseArray<>& AhD_dirich,
                       Gedim::Eigen_Array<>& rightHandSide)
{
    Gedim::MeshMatricesDAO blockMesh = *this->hD_Mesh;
    Eigen::VectorXi pivot = *this->hD_Pivot;
    Eigen::SparseMatrix<unsigned int> toEnrich_nodes = *this->toEnrich_nodes;

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {
        // Cycle on the test function index
        for (unsigned short int i = 0; i < 3; i++)
        {

            int ii_std = pivot[blockMesh.Cell3DVertex(e, i)];

            if (ii_std > 0)
            {

                if (toEnrich_nodes.coeff(ii_std, 0))
                {

                    // Computing the index of the enriched DOF corresponding to (e,i) node.
                    unsigned int numEnrNodesUpToNow = 0;
                    for (unsigned int temp_ind = 0; temp_ind < ii_std; temp_ind++)
                    {
                        numEnrNodesUpToNow += toEnrich_nodes.coeff(temp_ind, 0);
                    }
                    unsigned int ii_enr = pivot.size() + numEnrNodesUpToNow;


                    // Cycling on the trial function index
                    for (unsigned short int j = 0; j < 3; j++)
                    {

                        int jj_std = pivot[blockMesh.Cell3DVertex(e, j)];

                        if (jj_std > 0)
                        {

                            if (toEnrich_nodes.coeff(jj_std, 0))
                            {

                                // Computing the index of the enriched DOF corresponding to (e,j) node.
                                numEnrNodesUpToNow = 0;
                                for (unsigned int temp_ind = 0; temp_ind < ii_std; temp_ind++)
                                {
                                    numEnrNodesUpToNow += toEnrich_nodes.coeff(temp_ind, 0);
                                }
                                unsigned int jj_enr = pivot.size() + numEnrNodesUpToNow;


                                /* In questo caso, devi aggiungere elementi ad A
                                 *   (std test vs std trial), B (enr test vs std
                                 *   trial), A_tilde (std test vs enr trial),
                                 *   B_tilde (enr test vs enr trial).
                                 */

                                this->constructElement_AhD(AhD, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, e, integrationType::enr_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, e, integrationType::std_enr);

                                this->constructElement_AhD(AhD, ii_enr, jj_enr, e, integrationType::enr_enr);

                            }
                            else
                            {
                                /* In questo caso, devi aggiungere elementi ad A
                                   (std test vs std trial) e B (enr test vs std
                                   trial). */

                                this->constructElement_AhD(AhD, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, e, integrationType::enr_std);
                            }

                        }
                        else  // jj_std è nodo di Dirichlet
                        {
                            /*  In questo caso, hai un nodo di Dirichlet su
                                jj_std. Hai fatto le cose in modo tale che i
                                nodi di Dirichlet non portino arricchimenti.
                                Quindi, metti nel rilevamento delle
                                condizioni di Dirichlet solo per i DOF
                                jj_std, che hanno però effetto sulle righe
                                ii_enr. */

                            this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, e, integrationType::std_std);

                            this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, e, integrationType::enr_std);
                        }
                    }

                    this->constructEl_RHS_A(rightHandSide, ii_std, *this->geometryUtilities);
                    this->constructEl_RHS_A(rightHandSide, ii_enr, *this->geometryUtilities);


                } else
                {
                    for (unsigned short int j = 0; j < 3; j++)
                    {

                        int jj_std = pivot[blockMesh.Cell3DVertex(e, j)];

                        if (jj_std > 0)
                        {

                            if (toEnrich_nodes.coeff(jj_std, 0))
                            {

                                // Computing the index of the enriched DOF corresponding to (e,j) node.
                                unsigned int numEnrNodesUpToNow = 0;
                                for (unsigned int temp_ind = 0; temp_ind < ii_std; temp_ind++)
                                {
                                    numEnrNodesUpToNow += toEnrich_nodes.coeff(temp_ind, 0);
                                }
                                unsigned int jj_enr = pivot.size() + numEnrNodesUpToNow;

                                /* In questo caso, devi aggiungere elementi ad A
                                   (std test vs std trial), A_tilde (std test vs enr trial). */

                                this->constructElement_AhD(AhD, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, e, integrationType::std_enr);
                            }
                            else
                            {
                                /* In questo caso, devi aggiungere elementi ad A
                                   (std test vs std trial). */

                                this->constructElement_AhD(AhD, ii_std, jj_std, e, integrationType::std_std);
                            }

                        }
                        else  // jj_std è nodo di Dirichlet
                        {
                            /*  In questo caso, hai un nodo di Dirichlet su
                                jj_std. Hai fatto le cose in modo tale che i
                                nodi di Dirichlet non portino arricchimenti.
                                Quindi, metti nel rilevamento delle
                                condizioni di Dirichlet solo per i DOF
                                jj_std. */

                            this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, e, integrationType::std_std);

                        }

                    }

                    this->constructEl_RHS_A(rightHandSide, ii_std, *this->geometryUtilities);

                }

            }
        }
    }

}


void P1MatrixAssembler::constructElement_AhD(Gedim::Eigen_SparseArray<>& M,
                                               const unsigned int row,
                                               const unsigned int col,
                                               const unsigned int elementIndex,
                                               const integrationType type)
{
    // We are assembling the A^hD matrix, therefore we use the mesh for the hD variable.
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;

    // Identification of the current tetrahedron
    const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(mesh, elementIndex);

    // Identification of the current test (i-th) basis function
    Eigen::Vector4d testLagrangeCoeff = Utilities::lagrangeBasisCoeff(mesh.Cell0DCoordinates(row), element);

    // Identification of the current trial (j-th) basis function
    Eigen::Vector4d trialLagrangeCoeff = Utilities::lagrangeBasisCoeff(mesh.Cell0DCoordinates(col), element);

    // Identification of "test node"'s coordinates
    Eigen::Vector3d x_jCoord = mesh.Cell0DCoordinates(row);

    // Identification of "trial node"'s coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(col);

    // ******************************************************************************************************************************
    // Integration: different integration routines for the different combinations of DOF types that may arise.
    double integralValue = 0.0;
    switch (type)
    {

    case integrationType::std_std:
    {

        // Preparation of the quadrature formula

        const unsigned int quadratureOrder = 1;
        Eigen::MatrixXd quadraturePointsRef;
        Eigen::VectorXd quadratureWeightsRef;
        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                    quadraturePointsRef,
                                                                    quadratureWeightsRef);

        Gedim::MapTetrahedron mapping(*geometryUtilities);

        // Mapping the reference quadrature nodes and weights on the tetrahedron

        const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(element.Vertices);

        Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                                 quadraturePointsRef);

        Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                    quadraturePointsRef).array().abs();

        // Computing the integral

        double integralValue = 0.0;

        Eigen::Vector3d GradPhi_i = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
        for (unsigned int q = 0; q < mappedPoints.size(); q++)
        {
            integralValue += mappedWeights(q) * GradPhi_i.transpose() * GradPhi_j;
        }

        break;
    }

    case integrationType::enr_std:
    {

        // Preparation of the quadrature formula

        const unsigned int quadratureOrder = 1;
        Eigen::MatrixXd quadraturePointsRef;
        Eigen::VectorXd quadratureWeightsRef;
        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                    quadraturePointsRef,
                                                                    quadratureWeightsRef);

        Gedim::MapTetrahedron mapping(*geometryUtilities);

        // Mapping the reference quadrature nodes and weights on the tetrahedron

        const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(element.Vertices);

        Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                                 quadraturePointsRef);

        Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                    quadraturePointsRef).array().abs();

        // Computing the integral

        Eigen::Vector3d GradPhi_i = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
        double scalarProduct = GradPhi_j.transpose() * GradPhi_i;

        double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(x_iCoord, *fracture));

        for (unsigned int q = 0; q < mappedPoints.size(); q++)
        {
            double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
            integralValue += mappedWeights(q) * ((Psi_q - Psi_i) * scalarProduct);
        }

        break;
    }

    case integrationType::std_enr:
    {

        // Preparation of the quadrature formula

        const unsigned int quadratureOrder = 1;
        Eigen::MatrixXd quadraturePointsRef;
        Eigen::VectorXd quadratureWeightsRef;
        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                    quadraturePointsRef,
                                                                    quadratureWeightsRef);

        Gedim::MapTetrahedron mapping(*geometryUtilities);

        // Mapping the reference quadrature nodes and weights on the tetrahedron

        const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(element.Vertices);

        Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                                 quadraturePointsRef);

        Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                    quadraturePointsRef).array().abs();

        // Computing the integral

        Eigen::Vector3d GradPhi_i = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
        double scalarProduct = GradPhi_j.transpose() * GradPhi_i;

        double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(x_jCoord, *fracture));

        for (unsigned int q = 0; q < mappedPoints.size(); q++)
        {
            double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
            integralValue += mappedWeights(q) * ((Psi_q - Psi_j) * scalarProduct);
        }

        break;
    }

    case integrationType::enr_enr:
    {

        // Auxiliary geometric variables
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

        // Splitting the mesh tetrahedron into sub-polyhedra
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

        std::vector<Gedim::GeometryUtilities::Polyhedron> subPolyhedra = geometryUtilities->SplitPolyhedronWithPlaneResultToPolyhedra(result);


        // Meshing each sub-polyhedron with tetrahedral elements in order to use a quadrature formula based on tetrahedral ref. element

        std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra;

        for (unsigned int p = 0; p < subPolyhedra.size(); p++)
        {
            if (subPolyhedra.at(p).Vertices.size() == 4) // the polyhedron is already a tetrahedron
            {
                subTetrahedra.push_back(subPolyhedra.at(p));
            }
            else                                         // the polyhedron is not a tetrahedron
            {
                Gedim::MeshMatrices meshData;
                Gedim::MeshMatricesDAO mesh(meshData);

                meshUtilities->CreateTetrahedralMesh(subPolyhedra.at(p).Vertices,
                                                    subPolyhedra.at(p).Edges,
                                                    subPolyhedra.at(p).Faces,
                                                    1,
                                                    mesh,
                                                    "Qpfezna");

                for (unsigned short int t = 0; t < mesh.Cell3DTotalNumber(); t++)
                {
                    subTetrahedra.push_back(meshUtilities->MeshCell3DToPolyhedron(mesh, t));
                }

            }

        }


        // Integration of function in every subtetrahedron
        // 1) Creation of the reference element (useful?)
        const Gedim::GeometryUtilities::Polyhedron refTet = geometryUtilities->CreateTetrahedronWithVertices(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                            Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                            Eigen::Vector3d(0.0, 1.0, 0.0),
                                                                                                            Eigen::Vector3d(0.0, 0.0, 1.0));
        // 2) Determination of quadrature nodes and weigths.
        const unsigned int quadratureOrder = 1;
        Eigen::MatrixXd quadraturePointsRef;
        Eigen::VectorXd quadratureWeightsRef;
        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                    quadraturePointsRef,
                                                                    quadratureWeightsRef);

        Gedim::MapTetrahedron mapping(*geometryUtilities);

        // 3) Computing the integral on each sub-tetrahedron
        for (unsigned short int t = 0; t < subTetrahedra.size(); t++)
        {
            // Mapping the reference quadrature nodes and weights on the t-th sub-tetrahedra
            const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(element.Vertices);

            Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                                     quadraturePointsRef);

            Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                  quadraturePointsRef).array().abs();

            Eigen::Vector3d GradPhi_i = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
            Eigen::Vector3d GradPhi_j = Utilities::evaluateLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
            double scalarProduct = GradPhi_j.transpose() * GradPhi_i;

            double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(x_iCoord, *fracture));
            double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(x_jCoord, *fracture));

            for (unsigned int q = 0; q < mappedPoints.size(); q++)
            {
                double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
                integralValue += mappedWeights(q) * ((Psi_q - Psi_i) * (Psi_q - Psi_j) * scalarProduct);
            }
        }


        break;
    }

    }
    // ******************************************************************************************************************************


    // Placing the result in the correct entrance of the stiffness matrix
    // TODO: check whether .Triplet sums the contribution of integralValue or replaces it.

    M.Triplet(row, col, integralValue);

}



void P1MatrixAssembler::constructEl_RHS_A(Gedim::Eigen_Array<> &rhs,
                                             unsigned int i,
                                             const Gedim::GeometryUtilities &geometryUtilities)
{

}

void P1MatrixAssembler::constructEl_RHS_B(Gedim::Eigen_Array<> &rhs,
                                             unsigned int i,
                                             const Gedim::GeometryUtilities &geometryUtilities)
{

}




}
