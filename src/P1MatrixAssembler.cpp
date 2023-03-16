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

void P1MatrixAssembler::initialize(unsigned int numDOF_3D_std, unsigned int numDirich_3D)
{
    // This method has to be called after setting all the meshes.
    bool incompleteSettings = (this->hD_Mesh == NULL) || (this->hF_Mesh == NULL);

    if (incompleteSettings)
    {
        cerr << "Error: you have to set some meshes before initializing their coupling!";
        return;
    }

    // Initializing the enrichment information
    this->initializeEnrichmentInformation(numDOF_3D_std, numDirich_3D);

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

void P1MatrixAssembler::initializeEnrichmentInformation(unsigned int numDOF_3D_std, unsigned int numDirich_3D)
{
    /*
     * Questa funzione termina l'inizializzazione del vettore pivot per la variabile hD, legando i gradi di libertà standard
     * al corrispondente grado di libertà arricchito. Il tutto viene fatto trovando dapprima le celle tagliate dal piano della
     * frattura, per poi assegnare una numerazione ai DOF arricchiti sui punti di tali celle.
     * */

    Gedim::MeshMatricesDAO* mesh3D = this->hD_Mesh;
    Eigen::MatrixXi*       pivot3D = this->hD_Pivot;


    // Vettore di 0,1 (false, true). Numero di righe = numero di celle. Ci dice una data cella è tagliata
    // dalla frattura (quindi da arricchire) oppure no.
    this->toEnrich_elements = new Eigen::SparseMatrix<unsigned int>(mesh3D->Cell3DTotalNumber(), 1);


    // Ogni volta che troviamo una cella da arricchire, si aggiungono sempre esattamente 4 gradi di libertà di arricchimento.
    // Per numerare i DOF arricchiti è dunque sufficiente 'appendere' in coda alla lista di DOF gli arricchimenti.
    unsigned int global_counter_for_enriched_nodes = 0;


    for (unsigned int e = 0; e < mesh3D->Cell3DTotalNumber(); e++)
    {
        const Gedim::GeometryUtilities::Polyhedron elementAsPolyhedron = this->meshUtilities->MeshCell3DToPolyhedron(*mesh3D, e);

        if (this->fracture->intersects(elementAsPolyhedron))
        {
            // Setto l'elemento 'e' come da arricchire.
            this->toEnrich_elements->insert(e, 0) = 1;

            // Recupero l'Id globale dei punti del tetraedro corrente.
            unsigned int globIdNode1 = mesh3D->Cell3DVertex(e, 0),
                         globIdNode2 = mesh3D->Cell3DVertex(e, 1),
                         globIdNode3 = mesh3D->Cell3DVertex(e, 2),
                         globIdNode4 = mesh3D->Cell3DVertex(e, 3);

            // Calcolo l'indice del loro DOF di arricchimento.
            unsigned int nn1_enr, nn2_enr, nn3_enr, nn4_enr;

            bool node_1_still_to_be_numbered = (*pivot3D)(globIdNode1, 1) < 0,
                 node_2_still_to_be_numbered = (*pivot3D)(globIdNode2, 1) < 0,
                 node_3_still_to_be_numbered = (*pivot3D)(globIdNode3, 1) < 0,
                 node_4_still_to_be_numbered = (*pivot3D)(globIdNode4, 1) < 0;

            if (node_1_still_to_be_numbered)
            {
                global_counter_for_enriched_nodes++;
                nn1_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode1, 1) = nn1_enr;
            }

            if (node_2_still_to_be_numbered)
            {
                global_counter_for_enriched_nodes++;
                nn2_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode2, 1) = nn2_enr;
            }

            if (node_3_still_to_be_numbered)
            {
                global_counter_for_enriched_nodes++;
                nn3_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode3, 1) = nn3_enr;
            }

            if (node_4_still_to_be_numbered)
            {
                global_counter_for_enriched_nodes++;
                nn4_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode4, 1) = nn4_enr;
            }

        }
    }

    this->num_enrichments = global_counter_for_enriched_nodes;
}



// SETTERS ********************************************************************************************************

void P1MatrixAssembler::setHD_Mesh(Gedim::MeshMatricesDAO *newHD_Mesh)
{
    hD_Mesh = newHD_Mesh;
}

void P1MatrixAssembler::setHD_Pivot(Eigen::MatrixXi *newHD_Pivot)
{
    hD_Pivot = newHD_Pivot;
}

void P1MatrixAssembler::setHD_NeumannInfo(std::vector<unsigned int> *neumannInfo)
{
    this->neumannInfo = neumannInfo;
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
    Eigen::MatrixXi pivot = *this->hD_Pivot;

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {
        // Cycle on the test function index
        for (unsigned short int i = 0; i < 4; i++)
        {

            int globIdNode_i = blockMesh.Cell3DVertex(e, i);
            int ii_std = pivot(globIdNode_i, 0);

            bool i_node_is_DOF = ii_std > 0;
            bool i_node_is_to_enrich =  pivot(globIdNode_i, 1) > 0;

            if (i_node_is_DOF)
            {

                if (i_node_is_to_enrich)
                {
                    unsigned int ii_enr = pivot(globIdNode_i, 1);

                    // Cycling on the trial function index
                    for (unsigned short int j = 0; j < 4; j++)
                    {
                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);

                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;

                        if (j_node_is_DOF)
                        {

                            if (j_node_is_to_enrich)
                            {
                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                this->constructElement_AhD(AhD, ii_std, jj_std, globIdNode_i, globIdNode_j, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, globIdNode_i, globIdNode_j, e, integrationType::enr_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::std_enr);

                                this->constructElement_AhD(AhD, ii_enr, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::enr_enr);

                            }
                            else
                            {
                                this->constructElement_AhD(AhD, ii_std, jj_std, globIdNode_i, globIdNode_j, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, globIdNode_i, globIdNode_j, e, integrationType::enr_std);
                            }

                        }
                        else  // Il nodo (e,j) è di Dirichlet. Può ancora essere da arricchire e portare contributo coi suoi DOF enr.
                        {
                            // this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);

                            // this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, ii_std, -jj_std, e, integrationType::enr_std);

                            if (j_node_is_to_enrich)
                            {
                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::std_enr);

                                this->constructElement_AhD(AhD, ii_enr, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::enr_enr);
                            }
                        }
                    }

                    // Il nodo i è da arricchire. Quindi, costruiamo il RHS sia per il DOF std che per il DOF enr.

                    this->constructElement_Rhs(rightHandSide, ii_std, globIdNode_i, e, integrationType::std);
                    this->constructElement_Rhs(rightHandSide, ii_enr, globIdNode_i, e, integrationType::enr);

                }

                else // Il nodo (e,i) NON è da arricchire.
                {
                    for (unsigned short int j = 0; j < 4; j++)
                    {
                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);

                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;

                        if (j_node_is_DOF)
                        {
                            if (j_node_is_to_enrich)
                            {
                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                this->constructElement_AhD(AhD, ii_std, jj_std, globIdNode_i, globIdNode_j, e, integrationType::std_std);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::std_enr);
                            }
                            else
                            {
                                this->constructElement_AhD(AhD, ii_std, jj_std, globIdNode_i, globIdNode_j, e, integrationType::std_std);
                            }
                        }

                        else // Il nodo (e,j) è di Dirichlet
                        {
                            // this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);

                            if (j_node_is_to_enrich)
                            {
                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                this->constructElement_AhD(AhD, ii_std, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::std_enr);
                            }

                        }

                    }

                    // Il nodo i non è da arricchire. Costruiamo il RHS solo sul DOF std.
                    this->constructElement_Rhs(rightHandSide, ii_std, globIdNode_i, e, integrationType::std);

                }

            }

            else // Il nodo (e,i) non è un DOF (è una condizione di Dirichlet). Può essere da arricchire o no.
            {
                if (i_node_is_to_enrich)
                {
                    unsigned int ii_enr = pivot(globIdNode_i, 1);

                    // Cycling on the trial function index
                    for (unsigned short int j = 0; j < 4; j++)
                    {
                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);

                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;

                        if (j_node_is_DOF)
                        {

                            if (j_node_is_to_enrich)
                            {

                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                this->constructElement_AhD(AhD, ii_enr, jj_std, globIdNode_i, globIdNode_j, e, integrationType::enr_std);

                                this->constructElement_AhD(AhD, ii_enr, jj_enr, globIdNode_i, globIdNode_j, e, integrationType::enr_enr);

                            }
                            else
                            {
                                this->constructElement_AhD(AhD, ii_enr, jj_std, globIdNode_i, globIdNode_j, e, integrationType::enr_std);
                            }

                        }
                        else  // jj_std è nodo di Dirichlet
                        {
                            //  this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, -ii_std, -jj_std, e, integrationType::std_std);

                            //  this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, -ii_std, -jj_std, e, integrationType::enr_std);

                        }
                    }

                    // Il nodo i ha una condizione di Dirichlet ma è da arricchire. Costruiamo il RHS solo sul DOF enr.
                    this->constructElement_Rhs(rightHandSide, ii_enr, globIdNode_i, e, integrationType::enr);

                }
            }
        }
    }

}


void P1MatrixAssembler::constructElement_AhD(Eigen::SparseMatrix<double>& M,
                                              unsigned int row,
                                              unsigned int col,
                                              unsigned int globIdNode_i,
                                              unsigned int globIdNode_j,
                                             const unsigned int elementIndex,
                                             const integrationType type)
{    

    /*
     * Gli indici dei DOF/nodi Dirich. sono numerati in pivot seguendo la seguente convenzione:
     *      1,2,3,4,5...      per i DOF
     *      -1,-2,-3,-4,-5... per le condizioni di Dirichlet
     * Ora, voglio mettere il contributo del primo DOF nella prima riga della matrice di rigidezza, ovvero la riga 0.
     * Di conseguenza, decremento gli indici per avere la loro posizione nella matrice.
     * */
    row--;
    col--;

    #ifdef __DEBUG__
    std::cout << "Entered constructElement function." << std::endl;
    #endif
    // We are assembling the A^hD matrix, therefore we use the mesh for the hD variable.
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;

    // Identification of the current tetrahedron
    const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(mesh, elementIndex);

    // Identification of "test node"'s coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(globIdNode_i);

    // Identification of "trial node"'s coordinates
    Eigen::Vector3d x_jCoord = mesh.Cell0DCoordinates(globIdNode_j);

    // Identification of the current test (i-th) basis function
    #ifdef __DEBUG__
    std::cout << "Starting determination of lagrange test basis coeffs..." << std::endl;
    #endif
    Eigen::Vector4d testLagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(globIdNode_i), element);

    // Identification of the current trial (j-th) basis function
    #ifdef __DEBUG__
    std::cout << "Starting determination of lagrange trial basis coeffs..." << std::endl;
    #endif
    Eigen::Vector4d trialLagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(globIdNode_j), element);


    // ******************************************************************************************************************************
    // Integration: different integration routines for the different combinations of DOF types that may arise.
    double integralValue = 0.0;
    switch (type)
    {

    case integrationType::std_std:
    {
        #ifdef __DEBUG__
        std::cout << "Starting a std_std integration... position in matrix: " << row << ", " << col << "\n" << std::endl;
        #endif

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

        for (unsigned int q = 0; q < mappedPoints.cols(); q++)
        {
            integralValue += mappedWeights(q) * scalarProduct;
        }

        break;
    }

    case integrationType::enr_std:
    {
        #ifdef __DEBUG__
        std::cout << "Starting a enr_std integration... position in matrix: " << row << ", " << col << "\n" << std::endl;
        #endif

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
        #ifdef __DEBUG__
        std::cout << "Starting a std_enr integration... position in matrix: " << row << ", " << col << "\n" << std::endl;
        #endif

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
        #ifdef __DEBUG__
        std::cout << "Starting a enr_enr integration...position in matrix: " << row << ", " << col << "\n" << std::endl;
        #endif

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
                    Gedim::MeshMatrices aux_meshData;
                    Gedim::MeshMatricesDAO aux_mesh(aux_meshData);

                    meshUtilities->CreateTetrahedralMesh(subPolyhedra.at(p).Vertices,
                                                         subPolyhedra.at(p).Edges,
                                                         subPolyhedra.at(p).Faces,
                                                         1,
                                                         aux_mesh,
                                                         "Qpfezna");



                    for (unsigned short int t = 0; t < aux_mesh.Cell3DTotalNumber(); t++)
                    {
                        subTetrahedra.push_back(meshUtilities->MeshCell3DToPolyhedron(aux_mesh, t));
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
            const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(subTetrahedra.at(t).Vertices);

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
    // ******************************************************************************************************************************

    // Placing the result in the correct entrance of the stiffness matrix
    M.coeffRef(row, col) += integralValue;

}

void P1MatrixAssembler::constructElement_Rhs(Eigen::VectorXd&      rhs,
                                             unsigned int          row,
                                             const unsigned int    globIdNode_i,
                                             const unsigned int    elementIndex,
                                             const integrationType type)
{
    /*
     * Gli indici dei DOF/nodi Dirich. sono numerati in pivot seguendo la seguente convenzione:
     *      1,2,3,4,5...      per i DOF
     *      -1,-2,-3,-4,-5... per le condizioni di Dirichlet
     * Ora, voglio mettere il contributo del primo DOF nella prima riga della matrice di rigidezza, ovvero la riga 0.
     * Di conseguenza, decremento gli indici per avere la loro posizione nella matrice.
     * */
    row--;

    // The RHS is non zero only for the last row of the block linear system. Therefore, we use the hD mesh.
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;

    // Identification of the current tetrahedron
    const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(mesh, elementIndex);

    // Identification of the current (i-th) basis function
    Eigen::Vector4d lagrangeCoeff = Utilities::lagrangeBasisCoeff3D(mesh.Cell0DCoordinates(globIdNode_i), element);

    // Identification of the current node coordinates
    Eigen::Vector3d x_iCoord = mesh.Cell0DCoordinates(globIdNode_i);

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


        rhs(row) += integralValue;

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


        rhs(row) += integralValue;

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

void P1MatrixAssembler::addNeumann(Eigen::VectorXd &rhs)
{
    Gedim::MeshMatricesDAO mesh = *this->hD_Mesh;
    Eigen::MatrixXi pivot = *this->hD_Pivot;

    for (unsigned int t = 0; t < this->neumannInfo->size(); t++)
    {
        unsigned int triangle_glob_id = this->neumannInfo->at(t);
        unsigned int marker = mesh.Cell2DMarker(triangle_glob_id);
        std::vector<unsigned int> trianglePoints_Ids = mesh.Cell2DVertices(triangle_glob_id);
        Eigen::MatrixXd trianglePoints_coordinates = mesh.Cell2DVerticesCoordinates(triangle_glob_id);

        Eigen::MatrixXd rotatedTo2D_TrianglePoints_coordinates = this->geometryUtilities->RotatePointsFrom3DTo2D(trianglePoints_coordinates,
                                                                                                                 this->geometryUtilities->PlaneRotationMatrix(this->geometryUtilities->PolygonNormal(trianglePoints_coordinates)));

        double triangle_area = this->geometryUtilities->PolygonArea(rotatedTo2D_TrianglePoints_coordinates);

        double averaged_gN = DiscontinousTestProblem_1::g_Neumann_triangle_averaged(trianglePoints_Ids,
                                                                                    marker,
                                                                                    mesh,
                                                                                    this->fracture);

        for (unsigned int i = 0; i < 3; i++)
        {
            unsigned int current_point_id = trianglePoints_Ids.at(i);
            int ii = pivot(current_point_id, 0);

            if (ii > 0)
                rhs(ii) += (1 / 3) * averaged_gN * triangle_area;
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
