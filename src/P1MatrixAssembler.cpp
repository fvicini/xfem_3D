#include "P1MatrixAssembler.hpp"
#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "MapTetrahedron.hpp"
#include "MapTriangle.hpp"
#include "Utilities.hpp"

namespace XFEM_3D
{

// CONSTRUCTOR ***************************************************************************************************

P1MatrixAssembler::P1MatrixAssembler(Fracture3D* fracture,
                             Gedim::GeometryUtilities* geometryUtilities,
                             Gedim::MeshUtilities* meshutilities) {

        this->fracture = fracture;
        this->geometryUtilities = geometryUtilities;
        this->meshUtilities = meshutilities;

    }

// Initialization of the assembler ********************************************************************************

void P1MatrixAssembler::initialize()
{
    // This method has to be called after setting all the meshes.
    bool incompleteSettings = (this->hD_Mesh == NULL) || (this->hF_Mesh == NULL);

    if (incompleteSettings)
    {
        cerr << "Error: you have to set some meshes before initializing their coupling!";
        return;
    }

    // Initializing the enrichment information
    this->initializeEnrichmentInformation();

    // Initializing the coupling for hD, PsiP variables
    //this->initialize_2D_3DCoupling(fractureBorder::positive);

    // Initializing the coupling for hD, PsiM variables
    //this->initialize_2D_3DCoupling(fractureBorder::negative);

    // Initializing the coupling for hD, PsiF variables
    //this->initialize_2D_3DCoupling(fractureBorder::fracture);

}

void P1MatrixAssembler::initialize_2D_3DCoupling(fractureBorder type)
{
//    Gedim::MeshMatricesDAO* mesh3D;
//    Gedim::MeshMatricesDAO* mesh2D;
//    Eigen::SparseMatrix<unsigned int>* intersectionMatrix;

//    switch (type) {
//    case fractureBorder::positive:
//    {
//        mesh3D = this->hD_Mesh;
//        mesh2D = this->psiP_Mesh;
//        intersectionMatrix = this->hD_PsiP_MeshIntersections;

//        break;
//    }
//    case fractureBorder::negative:
//    {
//        mesh3D = this->hD_Mesh;
//        mesh2D = this->psiM_Mesh;
//        intersectionMatrix = this->hD_PsiM_MeshIntersections;

//        break;
//    }

//    case fractureBorder::fracture:
//    {
//        mesh3D = this->hD_Mesh;
//        mesh2D = this->psiF_Mesh;
//        intersectionMatrix = this->hD_PsiM_MeshIntersections;

//        break;
//    }
//    }


//    for (unsigned int e = 0; e < mesh3D->Cell3DTotalNumber(); e++)
//    {
//        if (!this->toEnrich_elements->coeff(e, 0))
//            continue;

//        // First column is global index in the mesh3D of the current element.
//        intersectionMatrix->coeffRef(e, 0) = e;

//        // Auxiliary geometric variables
//        Gedim::GeometryUtilities::Polyhedron elementAsPolyhedron = meshUtilities->MeshCell3DToPolyhedron(*mesh3D, e);
//        std::vector<Eigen::MatrixXd>  polyhedronFaceVertices = geometryUtilities->PolyhedronFaceVertices(elementAsPolyhedron.Vertices,
//                                                                                                         elementAsPolyhedron.Faces);
//        std::vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities->PolyhedronFaceTranslations(polyhedronFaceVertices);
//        std::vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities->PolyhedronFaceNormals(polyhedronFaceVertices);
//        std::vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities->PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
//                                                                                                                        polyhedronFaceNormals,
//                                                                                                                        polyhedronFaceTranslations);
//        std::vector<Eigen::MatrixXd> polyhedronFaceRotatedVertices = geometryUtilities->PolyhedronFaceRotatedVertices(polyhedronFaceVertices,
//                                                                                                                      polyhedronFaceTranslations,
//                                                                                                                      polyhedronFaceRotationMatrices);
//        Eigen::Vector3d pointInsidePolyhedron = geometryUtilities->PolyhedronBarycenter(elementAsPolyhedron.Vertices);
//        std::vector<bool> polyhedronFaceNormalDirections = geometryUtilities->PolyhedronFaceNormalDirections(polyhedronFaceVertices,
//                                                                                                             pointInsidePolyhedron,
//                                                                                                             polyhedronFaceNormals);


//        // Cycle on 2D mesh and find which elements are intersected
//        // It is enough to find a node of the 2D element lying inside of the current 3D element
//        for (unsigned int t = 0; t < mesh2D->Cell2DTotalNumber(); t++)
//        {
//            // Identification of the current triangle as a vector of point global Ids.
//            vector<unsigned int> elementPointsIds = mesh2D->Cell2DVertices(t);

//            // Checking whether the triangle is intersected by the current 3D element.
//            for (unsigned int k = 0; k < elementPointsIds.size(); k++)
//            {

//                Gedim::GeometryUtilities::PointPolyhedronPositionResult result = geometryUtilities->PointPolyhedronPosition(mesh2D->Cell0DCoordinates(elementPointsIds.at(k)),
//                                                                                                                            elementAsPolyhedron.Faces,
//                                                                                                                            polyhedronFaceVertices,
//                                                                                                                            polyhedronFaceRotatedVertices,
//                                                                                                                            polyhedronFaceNormals,
//                                                                                                                            polyhedronFaceNormalDirections,
//                                                                                                                            polyhedronFaceTranslations,
//                                                                                                                            polyhedronFaceRotationMatrices);

//                if (result.Type == Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Inside)
//                {
//                    // The triangle is chosen
//                    intersectionMatrix->coeffRef(e, t+1) = 1;
//                }
//            }

//        }
//    }

}

void P1MatrixAssembler::initializeEnrichmentInformation()
{


    Gedim::MeshMatricesDAO* mesh3D = this->hD_Mesh;

    this->toEnrich_elements = new Eigen::SparseMatrix<unsigned int>(mesh3D->Cell3DTotalNumber(), 1);
    this->toEnrich_nodes    = new Eigen::SparseMatrix<unsigned int>(mesh3D->Cell0DTotalNumber(), 1);

    for (unsigned int e = 0; e < mesh3D->Cell3DTotalNumber(); e++)
    {
        const Gedim::GeometryUtilities::Polyhedron elementAsPolyhedron = this->meshUtilities->MeshCell3DToPolyhedron(*mesh3D, e);

        if (this->fracture->intersects(elementAsPolyhedron))
        {
            // Set (e,f) entry of matrix to true, meaning element e is to enrich wrt fracture f.
            // There is only one
            this->toEnrich_elements->insert(e, 0) = 1;

            // Retrieve the global Id of the verteces of element
            unsigned int idNodo1, idNodo2, idNodo3, idNodo4;

            idNodo1 = mesh3D->Cell3DVertex(e, 0);
            idNodo2 = mesh3D->Cell3DVertex(e, 1);
            idNodo3 = mesh3D->Cell3DVertex(e, 2);
            idNodo4 = mesh3D->Cell3DVertex(e, 3);

            this->toEnrich_nodes->coeffRef(idNodo1, 0) = 1;
            this->toEnrich_nodes->coeffRef(idNodo2, 0) = 1;
            this->toEnrich_nodes->coeffRef(idNodo3, 0) = 1;
            this->toEnrich_nodes->coeffRef(idNodo4, 0) = 1;

        }
    }

}



// SETTERS ********************************************************************************************************

void P1MatrixAssembler::setHD_Mesh(Gedim::MeshMatricesDAO *newHD_Mesh)
{
    hD_Mesh = newHD_Mesh;
}

void P1MatrixAssembler::setHD_Pivot(Eigen::VectorXi *newHD_Pivot)
{
    hD_Pivot = newHD_Pivot;
}

void P1MatrixAssembler::setHF_Mesh(Gedim::MeshMatricesDAO *newHF_Mesh)
{
    hF_Mesh = newHF_Mesh;
}

void P1MatrixAssembler::setHF_Pivot(Eigen::VectorXi *newHF_Pivot)
{
    hF_Pivot = newHF_Pivot;
}

void P1MatrixAssembler::setPhysicalParameters(PhysicalParameters *newPhysicalParameters)
{
    physicalParameters = newPhysicalParameters;
}


// PDE Discretization *****************************************************************


void P1MatrixAssembler::assemble_hD_hD(Eigen::SparseMatrix<double>& AhD,
                                       Eigen::SparseMatrix<double>& AhD_dirich,
                                       Eigen::SparseMatrix<double>& GhD,
                                       Eigen::SparseMatrix<double>& GhD_dirich,
                                       Eigen::VectorXd&             rightHandSide)
{
    Gedim::MeshMatricesDAO blockMesh = *this->hD_Mesh;
    Eigen::VectorXi pivot = *this->hD_Pivot;
    Eigen::SparseMatrix<unsigned int> toEnrich_nodes = *this->toEnrich_nodes;

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {
        // Cycle on the test function index
        for (unsigned short int i = 0; i < 3; i++)
        {

            //int ii_std = pivot[blockMesh.Cell3DVertex(e, i)];
            int ii_std = blockMesh.Cell3DVertex(e, i);

            if (ii_std+1 > 0)
            {

                if (toEnrich_nodes.coeff(ii_std, 0)) // ii node is to be enriched
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

                        //int jj_std = pivot[blockMesh.Cell3DVertex(e, j)];
                        int jj_std = blockMesh.Cell3DVertex(e, j);

                        if (jj_std+1 > 0)
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

                                this->constructElement_AhD(AhD, ii_std, jj_std, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, ii_std, jj_std, e, integrationType::enr_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, ii_std, jj_std, e, integrationType::std_enr);

                                this->constructElement_AhD(AhD, ii_enr, jj_enr, ii_std, jj_std, e, integrationType::enr_enr);

                            }
                            else
                            {
                                /* In questo caso, devi aggiungere elementi ad A
                                   (std test vs std trial) e B (enr test vs std
                                   trial). */

                                this->constructElement_AhD(AhD, ii_std, jj_std, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, ii_std, jj_std, e,  integrationType::enr_std);
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

                            this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);

                            this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, ii_std, -jj_std, e, integrationType::enr_std);
                        }
                    }

                    this->constructElement_Rhs(rightHandSide, ii_std, ii_std, e, integrationType::std);

                    this->constructElement_Rhs(rightHandSide, ii_enr, ii_std, e, integrationType::enr);


                } else // ii node is NOT to be enriched
                {
                    for (unsigned short int j = 0; j < 3; j++)
                    {

                        //int jj_std = pivot[blockMesh.Cell3DVertex(e, j)];
                        int jj_std = blockMesh.Cell3DVertex(e, j);

                        if (jj_std+1 > 0)
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

                                this->constructElement_AhD(AhD, ii_std, jj_std, ii_std, jj_std, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, ii_std, jj_std, e, integrationType::std_enr);
                            }
                            else
                            {
                                /* In questo caso, devi aggiungere elementi ad A
                                   (std test vs std trial). */

                                this->constructElement_AhD(AhD, ii_std, jj_std, ii_std, jj_std, e, integrationType::std_std);
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

                            this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);

                        }

                    }

                    this->constructElement_Rhs(rightHandSide, ii_std, ii_std, e, integrationType::std);

                }

            }
        }
    }

}


void P1MatrixAssembler::constructElement_AhD(Eigen::SparseMatrix<double>& M,
                                               const unsigned int row,
                                               const unsigned int col,
                                               const unsigned int ii_std,
                                               const unsigned int jj_std,
                                               const unsigned int elementIndex,
                                               const integrationType type)
{    

    std::cout << "Entered constructElement function." << std::endl;


    if (row == 17 && col == 9)
    {
        std::cout << "Caso problematico" << std::endl;
    }

    // We are assembling the A^hD matrix, therefore we use the mesh for the hD variable.
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;

    // Identification of the current tetrahedron
    const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(mesh, elementIndex);

    // Identification of the current test (i-th) basis function
    std::cout << "Starting determination of lagrange test basis coeffs..." << std::endl;
    Eigen::Vector4d testLagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(ii_std), element);

    // Identification of the current trial (j-th) basis function
    std::cout << "Starting determination of lagrange trial basis coeffs..." << std::endl;
    Eigen::Vector4d trialLagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(jj_std), element);

    // Identification of "test node"'s coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(ii_std);

    // Identification of "trial node"'s coordinates
    Eigen::Vector3d x_jCoord = mesh.Cell0DCoordinates(jj_std);

    // ******************************************************************************************************************************
    // Integration: different integration routines for the different combinations of DOF types that may arise.
    double integralValue = 0.0;
    switch (type)
    {

    case integrationType::std_std:
    {

        std::cout << "Starting a std_std integration... position in matrix: " << row << ", " << col << "\n" << std::endl;

        // Preparation of the quadrature formula
        const unsigned int quadratureOrder = 2;
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

        Eigen::Vector3d GradPhi_i = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);

        double scalarProduct = (this->physicalParameters->getPermeabilityTensorVolume() * GradPhi_i).transpose() * GradPhi_j;

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            integralValue += mappedWeights(q) * scalarProduct;
        }

        break;
    }

    case integrationType::enr_std:
    {

        if(jj_std != col)
        {
            std::cout << "Qualcosa non va: jj_std è diverso da col. Dovrebbe essere uguale" << std::endl;
        }

        std::cout << "Starting a enr_std integration... position in matrix: " << row << ", " << col << "\n" << std::endl;


        // Preparation of the quadrature formula

        const unsigned int quadratureOrder = 2;
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

        Eigen::Vector3d GradPhi_i = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
        double scalarProduct = (this->physicalParameters->getPermeabilityTensorVolume() * GradPhi_j).transpose() * GradPhi_i;

        double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(x_iCoord, *fracture));

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
            integralValue += mappedWeights(q) * ((Psi_q - Psi_i) * scalarProduct);
        }

        break;
    }

    case integrationType::std_enr:
    {
        std::cout << "Starting a std_enr integration... position in matrix: " << row << ", " << col << "\n" << std::endl;


        // Preparation of the quadrature formula

        const unsigned int quadratureOrder = 2;
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

        Eigen::Vector3d GradPhi_i = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
        Eigen::Vector3d GradPhi_j = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
        double scalarProduct = (this->physicalParameters->getPermeabilityTensorVolume() * GradPhi_j).transpose() * GradPhi_i;

        double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(x_jCoord, *fracture));

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
            integralValue += mappedWeights(q) * ((Psi_q - Psi_j) * scalarProduct);
        }

        break;
    }

    case integrationType::enr_enr:
    {
        std::cout << "Starting a enr_enr integration...position in matrix: " << row << ", " << col << "\n" << std::endl;

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
        // By now, we only know that the two nodes are enriched, but the element may be just blending and
        // not reproducing. If it is blending, in the discontinous XFEM, we just have to integrate on it
        // like nothing happened.
        std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra;
        switch (result.Type)
        {
        // The element is blending: touches a reproducing element
        case Gedim::GeometryUtilities::SplitPolyhedronWithPlaneResult::Types::None:
        {
            subTetrahedra.push_back(element);

            break;
        }

        default: // The element is reproducing
        {
            std::vector<Gedim::GeometryUtilities::Polyhedron> subPolyhedra = geometryUtilities->SplitPolyhedronWithPlaneResultToPolyhedra(result);


            // Meshing each sub-polyhedron with tetrahedral elements in order to use a quadrature formula based on tetrahedral ref. element
            for (unsigned int p = 0; p < subPolyhedra.size(); p++)
            {
                if (subPolyhedra.at(p).Vertices.cols() == 4) // the polyhedron is already a tetrahedron
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

            break;
        }
        }


        // Integration of function in every subtetrahedron

        // 1) Determination of quadrature nodes and weigths.
        const unsigned int quadratureOrder = 2;
        Eigen::MatrixXd quadraturePointsRef;
        Eigen::VectorXd quadratureWeightsRef;
        Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                    quadraturePointsRef,
                                                                    quadratureWeightsRef);

        Gedim::MapTetrahedron mapping(*geometryUtilities);

        // 2) Computing the integral on each sub-tetrahedron
        for (unsigned short int t = 0; t < subTetrahedra.size(); t++)
        {
            // Mapping the reference quadrature nodes and weights on the t-th sub-tetrahedra
            const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(element.Vertices);

            Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                                     quadraturePointsRef);

            Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                        quadraturePointsRef).array().abs();

            Eigen::Vector3d GradPhi_i = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
            Eigen::Vector3d GradPhi_j = Utilities::evaluate3DLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);

            double scalarProduct = (this->physicalParameters->getPermeabilityTensorVolume() * GradPhi_j).transpose() * GradPhi_i;

            double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(x_iCoord, *fracture));
            double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(x_jCoord, *fracture));

            for (unsigned int q = 0; q < mappedPoints.cols(); q++)
            {
                double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
                integralValue += mappedWeights(q) * ((Psi_q - Psi_i) * (Psi_q - Psi_j) * scalarProduct);
            }
        }


        break;
    }
    default:
        break;

    }

    // Placing the result in the correct entrance of the stiffness matrix
    // TODO: check whether .Triplet sums the contribution of integralValue or replaces it.

    M.coeffRef(row, col) += integralValue;

}

void P1MatrixAssembler::constructElement_Rhs(Eigen::VectorXd&      rhs,
                                             const unsigned int    i,
                                             const unsigned int    ii_std,
                                             const unsigned int    elementIndex,
                                             const integrationType type)
{
    // The RHS is non zero only for the last row of the block linear system. Therefore, we use the hD mesh.
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;

    // Identification of the current tetrahedron
    const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(mesh, elementIndex);

    // Identification of the current (i-th) basis function
    Eigen::Vector4d lagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(ii_std), element);

    // Identification of the current node coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(ii_std);

    switch (type) {
    case integrationType::std:
    {

        // Preparation of the quadrature formula
        const unsigned int quadratureOrder = 2;
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

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            double f_xq = this->physicalParameters->forcingTerm(mappedPoints.col(q), *this->fracture);
            double phi_xq = Utilities::evaluate3DLagrange(mappedPoints.col(q), lagrangeCoeff);
            integralValue += mappedWeights(q) * f_xq * phi_xq;
        }


        rhs(i) += integralValue;

        break;

    }
    case integrationType::enr:
    {
        // Preparation of the quadrature formula
        const unsigned int quadratureOrder = 2;
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
        double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(x_iCoord, *fracture));

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            double Psi_q = Utilities::heaviside(Utilities::signedDistanceFunction(mappedPoints.col(q), *fracture));
            double f_xq = this->physicalParameters->forcingTerm(mappedPoints.col(q), *fracture);
            double phi_xq = Utilities::evaluate3DLagrange(mappedPoints.col(q), lagrangeCoeff);
            integralValue += mappedWeights(q) * f_xq * phi_xq * (Psi_q - Psi_i);
        }


        rhs(i) += integralValue;

        break;
    }
    default:
        break;
    }










}


void P1MatrixAssembler::assemble_h_h(Eigen::SparseMatrix<double>& Ah,
                                     Eigen::SparseMatrix<double>& Ah_dirich,
                                     Eigen::SparseMatrix<double>& Gh,
                                     Eigen::SparseMatrix<double>& Gh_dirich,
                                     Eigen::SparseMatrix<double>& EF,
                                     Eigen::SparseMatrix<double>& EF_dirich,
                                     Eigen::SparseMatrix<double>& GPsiF,
                                     Eigen::SparseMatrix<double>& GPsiF_dirich)
{
    Gedim::MeshMatricesDAO mesh = *this->hF_Mesh;
    Eigen::VectorXi pivot = *this->hF_Pivot;

    for (unsigned int e = 0; e < mesh.Cell2DTotalNumber(); e++)
    {
        for (unsigned short int i = 0; i < 2; i++)
        {
            unsigned int ii = pivot[mesh.Cell2DVertex(e, i)];

            if (ii > 0)
            {
                for (unsigned short int j = 0; j < 2; j++)
                {
                    unsigned int jj = pivot[mesh.Cell2DVertex(e, j)];

                    if (jj > 0)

                        this->constructElement_AhF(Ah, ii, jj, e);
                    else

                        this->constructElement_AhF(Ah_dirich, ii, -jj, e);

                }
            }
        }
    }
}



void P1MatrixAssembler::constructElement_AhF(Eigen::SparseMatrix<double>& M,
                                             const unsigned int row,
                                             const unsigned int col,
                                             const unsigned int elementIndex)
{
    // We are assembling the A^hF matrix, therefore we use the mesh for the hF variable.
    Gedim::MeshMatricesDAO mesh = *this->hF_Mesh;

    // Identification of the current triangle as a vector of point global Ids. Retrieval of their coordinates.
    vector<unsigned int> elementPointsIds = mesh.Cell2DVertices(elementIndex);
    Eigen::MatrixXd elementPointsCoords;
    for (unsigned int k = 0; k < elementPointsIds.size(); k++)
    {
        elementPointsCoords.col(k) = mesh.Cell0DCoordinates(elementPointsIds.at(k));
    }

    // Identification of the current test (i-th) basis function
    Eigen::Vector3d testLagrangeCoeff = Utilities::lagrangeBasisCoeff2D(mesh.Cell0DCoordinates(row), elementPointsCoords);

    // Identification of the current trial (j-th) basis function
    Eigen::Vector3d trialLagrangeCoeff = Utilities::lagrangeBasisCoeff2D(mesh.Cell0DCoordinates(row), elementPointsCoords);

    // Identification of "test node"'s coordinates
    Eigen::Vector3d x_jCoord = mesh.Cell0DCoordinates(row);

    // Identification of "trial node"'s coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(col);

    // Integration
    // Preparation of the quadrature formula

    const unsigned int quadratureOrder = 1;
    Eigen::MatrixXd quadraturePointsRef;
    Eigen::VectorXd quadratureWeightsRef;
    Gedim::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(quadratureOrder,
                                                             quadraturePointsRef,
                                                             quadratureWeightsRef);


    Gedim::MapTriangle mapping;

    // Mapping the reference quadrature nodes and weights on the tetrahedron

    const Gedim::MapTriangle::MapTriangleData mapData = mapping.Compute(elementPointsCoords);

    Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                             quadraturePointsRef);

    Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                quadraturePointsRef).array().abs();

    // Computing the integral
    double integralValue = 0.0;

    Eigen::Vector2d GradPhi_i = Utilities::evaluate2DLagrangeGrad(mappedPoints.col(0), testLagrangeCoeff);
    Eigen::Vector2d GradPhi_j = Utilities::evaluate2DLagrangeGrad(mappedPoints.col(0), trialLagrangeCoeff);
    double scalarProduct = (this->physicalParameters->getPermeabilityTensorFracture() * GradPhi_j).transpose() * GradPhi_i;

    for (unsigned int q = 0; q < mappedPoints.size(); q++)
    {
        double Phi_i_q = Utilities::evaluate2DLagrange(mappedPoints.col(q), testLagrangeCoeff);
        double Phi_j_q = Utilities::evaluate2DLagrange(mappedPoints.col(q), trialLagrangeCoeff);

        integralValue += mappedWeights(q) * (scalarProduct)
                + mappedWeights(q) * 2 * this->physicalParameters->getNormalTransmissivityFracture() * Phi_i_q * Phi_j_q;

    }


    M.coeffRef(row, col) += integralValue;
}


// **************************************************************************************









}
