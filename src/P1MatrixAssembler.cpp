#include "P1MatrixAssembler.hpp"
#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "MapTetrahedron.hpp"
#include "Utilities.hpp"

#define DEBUG_CONSTANT 12
#define GEOMETRIC_TOLERANCE 1.0e-6

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

        /*
         * fatti ritornare da intersects il prodotto di min_SD * max_SD -> se è circa zero
         * elemento troppo vicino per arricchirlo. Togli i controlli su nodo_i_too_close che non hanno senso
        */
        double product = this->fracture->intersects(elementAsPolyhedron);
        bool fracture_intersects_element_not_too_close = fabs(product) >= GEOMETRIC_TOLERANCE && product < 0;



        if (e == 3)
        {
            std::cout << "stop" << std::endl;
        }

        if (fracture_intersects_element_not_too_close)
        {

            std::cout << "Elemento arricchito: " << e << std::endl;

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

            bool node_1_too_close_to_interface = false;//fabs(Utilities::signedDistanceFunction(mesh3D->Cell0DCoordinates(globIdNode1), *fracture)) < GEOMETRIC_TOLERANCE,
            bool node_2_too_close_to_interface = false;//fabs(Utilities::signedDistanceFunction(mesh3D->Cell0DCoordinates(globIdNode2), *fracture)) < GEOMETRIC_TOLERANCE,
            bool node_3_too_close_to_interface = false;//fabs(Utilities::signedDistanceFunction(mesh3D->Cell0DCoordinates(globIdNode3), *fracture)) < GEOMETRIC_TOLERANCE,
            bool node_4_too_close_to_interface = false;//fabs(Utilities::signedDistanceFunction(mesh3D->Cell0DCoordinates(globIdNode4), *fracture)) < GEOMETRIC_TOLERANCE;

            if (node_1_still_to_be_numbered && !node_1_too_close_to_interface)
            {
                global_counter_for_enriched_nodes++;
                nn1_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode1, 1) = nn1_enr;
            }

            if (node_2_still_to_be_numbered && !node_2_too_close_to_interface)
            {
                global_counter_for_enriched_nodes++;
                nn2_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode2, 1) = nn2_enr;
            }

            if (node_3_still_to_be_numbered && !node_3_too_close_to_interface)
            {
                global_counter_for_enriched_nodes++;
                nn3_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                (*pivot3D)(globIdNode3, 1) = nn3_enr;
            }

            if (node_4_still_to_be_numbered && !node_4_too_close_to_interface)
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

    Eigen::Matrix<double, 3, 4> GradPhi;

    GradPhi << -1, 1, 0, 0,
               -1, 0, 1, 0,
               -1, 0, 0, 1;

    Eigen::Matrix3d nu = this->physicalParameters->getPermeabilityTensorVolume();

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {

        bool elementCutByFracture = this->toEnrich_elements->coeff(e, 0);

        const Gedim::GeometryUtilities::Polyhedron element = meshUtilities->MeshCell3DToPolyhedron(blockMesh, e);
        Eigen::Vector3d tetrahedronBarycenter = geometryUtilities->PolyhedronBarycenter(element.Vertices);

        MappingFromReferenceTetrahedronInfo mapping = this->MapFromReferenceTetrahedron(element);
        double abs_detJ    = mapping.abs_detJ;
        Eigen::Matrix3d JJ = mapping.JJ;


        // Cycle on the test function index
        for (unsigned short int i = 0; i < 4; i++)
        {
            Eigen::Vector3d gradPhi_i = GradPhi.col(i);

            int globIdNode_i = blockMesh.Cell3DVertex(e, i);
            int ii_std = pivot(globIdNode_i, 0);

            bool i_node_is_DOF = ii_std > 0;
            bool i_node_is_to_enrich =  pivot(globIdNode_i, 1) > 0;

            if (i_node_is_DOF) // Comprende: 1A, 1B, 3A, 3B, 2, 4, 5, 6
            {

                if (i_node_is_to_enrich)
                {
                    unsigned int ii_enr = pivot(globIdNode_i, 1);

                    // Cycling on the trial function index
                    for (unsigned short int j = 0; j < 4; j++)                        
                    {
                        Eigen::Vector3d gradPhi_j = GradPhi.col(j);

                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);


                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;


                        if (j_node_is_DOF)
                        {
                            if (j_node_is_to_enrich)
                            {
                                if (elementCutByFracture)
                                {
                                    // CASISTICA 1A DELLA DOCUMENTAZIONE ============================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                           integralValue_se_Ahd = 0.0,
                                           integralValue_es_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    {
                                        double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                              *fracture)),
                                               Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                               *fracture));

                                        std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                     fracture,
                                                                                                                                                     this->geometryUtilities,
                                                                                                                                                     this->meshUtilities);
                                        for(auto& subTetrahedron : subTetrahedra)
                                        {
                                            Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                            double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                            Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                            subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                            element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                            element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                            1,                           1,                           1,                           1;

                                            double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                            integralValue_es_Ahd += (Psi - Psi_i) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                            integralValue_se_Ahd += (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                            integralValue_ee_Ahd += (Psi - Psi_i) * (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                        }

                                        AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                        AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;
                                        AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;
                                        AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;
                                    }
                                    // (FINE CASISTICA 1A) ===================================================================================================
                                }

                                else
                                {
                                    // CASISTICA 1B DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                           integralValue_se_Ahd = 0.0,
                                           integralValue_es_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;


                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                           Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture)),
                                           Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));


                                    integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_es_Ahd = (Psi - Psi_i) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_se_Ahd = (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_ee_Ahd = (Psi - Psi_i) * (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 1B) ===================================================================================================
                                }
                            }
                            else // Il nodo j non è da arricchire
                            {

                                if (elementCutByFracture)
                                {
                                    // CASISTICA 2A DELLA DOCUMENTAZIONE ==========================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                            integralValue_es_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                            Psi   = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter,
                                                                                                           *fracture));
                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);
                                    integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_es_Ahd += (Psi - Psi_i) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    }

                                    AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;

                                    // (FINE CASISTICA 2A) ========================================================================================================
                                }

                                else
                                {
                                    // CASISTICA 2B DELLA DOCUMENTAZIONE ==========================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                            integralValue_es_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                            Psi   = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter,
                                                                                                           *fracture));

                                    integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_es_Ahd = (Psi - Psi_i) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;

                                    // (FINE CASISTICA 2B) ========================================================================================================

                                }

                            }

                        }
                        else  // Il nodo j è di Dirichlet. Può ancora essere da arricchire e portare contributo coi suoi DOF enr.
                        {
                            /*
                             *   this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);
                             *
                             *   this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, ii_std, -jj_std, e, integrationType::enr_std);
                            */

                            if (j_node_is_to_enrich)
                            {
                                if (elementCutByFracture)
                                {
                                    // CASISTICA 3A DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_se_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                           Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);
                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_se_Ahd += (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                        integralValue_ee_Ahd += (Psi - Psi_i) * (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    }
                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 3A) ===================================================================================================
                                }
                                else
                                {
                                    // CASISTICA 3B DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_se_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                           Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));

                                    double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    integralValue_se_Ahd = (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_ee_Ahd = (Psi - Psi_i) * (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 3B) ===================================================================================================

                                }
                            }
                        }
                    }

                    // COSTRUZIONE RHS =======================================================================================================================
                    // ii_std
                    double averaged_f = this->physicalParameters->forcingTermAveragedOnTetrahedron(element.Vertices, *fracture);
                    rightHandSide(ii_std - 1) += averaged_f * abs_detJ * (1.0 / 24.0);

                    // ii_enr
                    double integralValue_rhs_enr = 0.0;
                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                          *fracture));
                    if (elementCutByFracture)
                    {
                        std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                     fracture,
                                                                                                                                     this->geometryUtilities,
                                                                                                                                     this->meshUtilities);
                        for(auto& subTetrahedron : subTetrahedra)
                        {
                            Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                            double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                            // TODO

                            Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                            subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                            element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                            element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                            1,                           1,                           1,                           1;

                            double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                            integralValue_rhs_enr += (Psi - Psi_i) * averaged_f * abs_detJ * (1.0 / 24.0);

                        }

                        rightHandSide(ii_enr - 1) = integralValue_rhs_enr;

                    }
                    else
                    {
                        double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i), *fracture));
                        double Psi   = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                        rightHandSide(ii_enr - 1) = (Psi - Psi_i) * averaged_f * abs_detJ * (1.0 / 24.0);
                    }
                    // (FINE COSTRUZIONE RHS) =================================================================================================================

                }

                else // Il nodo (e,i) NON è da arricchire.
                {
                    for (unsigned short int j = 0; j < 4; j++)
                    {
                        Eigen::Vector3d gradPhi_j = GradPhi.col(j);

                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);

                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;

                        if (j_node_is_DOF)
                        {
                            if (j_node_is_to_enrich)
                            {
                                unsigned int jj_enr = pivot(globIdNode_j, 1);

                                if (elementCutByFracture)
                                {
                                    // CASISTICA 4A DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                            integralValue_se_Ahd = 0.0;

                                    double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                          *fracture));

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);

                                    integralValue_ss_Ahd  = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_se_Ahd += (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    }

                                    AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;

                                    // (FINE CASISTICA 4A) ===================================================================================================
                                }

                                else
                                {
                                    // CASISTICA 4B DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_ss_Ahd = 0.0,
                                            integralValue_se_Ahd = 0.0;

                                    double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                          *fracture)),
                                            Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_se_Ahd = (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;
                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;

                                    // (FINE CASISTICA 4B) ===================================================================================================
                                }
                            }


                            else // Il nodo j non è da arricchire
                            {
                                // CASISTICA 5 DELLA DOCUMENTAZIONE =====================================================================================

                                double integralValue_ss_Ahd = 0.0;

                                integralValue_ss_Ahd = (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                AhD.coeffRef(ii_std - 1, jj_std - 1) += integralValue_ss_Ahd;

                                // (FINE CASISTICA 5) ===================================================================================================

                            }

                        }

                        else // Il nodo j è di Dirichlet
                        {
                            // this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, ii_std, -jj_std, e, integrationType::std_std);

                            if (j_node_is_to_enrich)
                            {

                                if (elementCutByFracture)
                                {
                                    // CASISTICA 6A DELLA DOCUMENTAZIONE =======================================================================================

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double integralValue_se_Ahd = 0.0;

                                    double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                          *fracture));

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);

                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_se_Ahd += (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    }

                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;

                                    // (FINE CASISTICA 6BA =====================================================================================================

                                }

                                else
                                {
                                    // CASISTICA 6B DELLA DOCUMENTAZIONE =======================================================================================

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double integralValue_se_Ahd = 0.0;

                                    double Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                          *fracture)),
                                            Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    integralValue_se_Ahd = (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_std - 1, jj_enr - 1) += integralValue_se_Ahd;

                                    // (FINE CASISTICA 6B) =====================================================================================================
                                }
                            }
                        }
                    }

                    // Costruzione RHS per il DOF ii_std
                    double averaged_f = this->physicalParameters->forcingTermAveragedOnTetrahedron(element.Vertices, *fracture);

                    rightHandSide(ii_std - 1) += averaged_f * abs_detJ * (1.0 / 24.0);

                }

            }

            else // Il nodo i non è un DOF. Comprende: 7A, 7B, 8, 9
            {
                if (i_node_is_to_enrich)
                {
                    unsigned int ii_enr = pivot(globIdNode_i, 1);

                    // Cycling on the trial function index
                    for (unsigned short int j = 0; j < 4; j++)
                    {
                        Eigen::Vector3d gradPhi_j = GradPhi.col(j);

                        int globIdNode_j = blockMesh.Cell3DVertex(e, j);
                        int jj_std = pivot(globIdNode_j, 0);

                        bool j_node_is_DOF = jj_std > 0;
                        bool j_node_is_to_enrich =  pivot(globIdNode_j, 1) > 0;

                        if (j_node_is_DOF)
                        {
                            if (j_node_is_to_enrich)
                            {
                                if (elementCutByFracture)
                                {

                                    // CASISTICA 7A DELLA DOCUMENTAZIONE ============================================================================================

                                    double integralValue_es_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                           Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);
                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_es_Ahd += (Psi - Psi_i) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                        integralValue_ee_Ahd += (Psi - Psi_i) * (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    }
                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 7A) ===================================================================================================
                                }

                                else
                                {
                                    // CASISTICA 7B DELLA DOCUMENTAZIONE ============================================================================================

                                    double integralValue_es_Ahd = 0.0,
                                           integralValue_ee_Ahd = 0.0;

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                           Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));
                                    double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    integralValue_es_Ahd = (Psi - Psi_i) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    integralValue_ee_Ahd = (Psi - Psi_i) * (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;
                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 7B) ===================================================================================================
                                }

                            }

                            else // Il nodo j non è da arricchire
                            {
                                if (elementCutByFracture)
                                {
                                    // CASISTICA 8A DELLA DOCUMENTAZIONE ==========================================================================================

                                    double integralValue_es_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture));

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);

                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_es_Ahd += (Psi - Psi_i) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    }
                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;

                                    // (FINE CASISTICA 8A) ========================================================================================================
                                }

                                else
                                {
                                    // CASISTICA 8B DELLA DOCUMENTAZIONE ==========================================================================================

                                    double integralValue_es_Ahd = 0.0;

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                            Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    integralValue_es_Ahd = (Psi - Psi_i) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_enr - 1, jj_std - 1) += integralValue_es_Ahd;

                                    // (FINE CASISTICA 8B) ========================================================================================================
                                }
                            }
                        }

                        else  // Anche j è nodo di Dirichlet.
                        {
                            //  this->constructElement_AhD(AhD_dirich, ii_std, -jj_std, -ii_std, -jj_std, e, integrationType::std_std);

                            //  this->constructElement_AhD(AhD_dirich, ii_enr, -jj_std, -ii_std, -jj_std, e, integrationType::enr_std);

                            if (j_node_is_to_enrich)
                            {
                                if (elementCutByFracture)
                                {
                                    // CASISTICA 9A DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_ee_Ahd = 0.0;

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                            Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));

                                    std::vector<Gedim::GeometryUtilities::Polyhedron> subTetrahedra = Utilities::splitTetrahedronInSubTetrahedra(element,
                                                                                                                                                 fracture,
                                                                                                                                                 this->geometryUtilities,
                                                                                                                                                 this->meshUtilities);
                                    for(auto& subTetrahedron : subTetrahedra)
                                    {
                                        Eigen::Vector3d subTetrahedronBarycenter = this->geometryUtilities->PolyhedronBarycenter(subTetrahedron.Vertices);
                                        double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(subTetrahedronBarycenter, *fracture));

                                        Eigen::Matrix<double, 4, 4> subTetrahedronVerticesMatrix;

                                        subTetrahedronVerticesMatrix << element.Vertices.col(0).x(), element.Vertices.col(1).x(), element.Vertices.col(2).x(), element.Vertices.col(3).x(),
                                                                        element.Vertices.col(0).y(), element.Vertices.col(1).y(), element.Vertices.col(2).y(), element.Vertices.col(3).y(),
                                                                        element.Vertices.col(0).z(), element.Vertices.col(1).z(), element.Vertices.col(2).z(), element.Vertices.col(3).z(),
                                                                        1,                           1,                           1,                           1;

                                        double subTetrahedronVolume = (1.0 / 6.0) * fabs(subTetrahedronVerticesMatrix.determinant());

                                        integralValue_ee_Ahd += (Psi - Psi_i) * (Psi - Psi_j) * subTetrahedronVolume * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;
                                    }

                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 9A) ===================================================================================================
                                }
                                else
                                {
                                    // CASISTICA 9B DELLA DOCUMENTAZIONE =====================================================================================

                                    double integralValue_ee_Ahd = 0.0;

                                    unsigned int jj_enr = pivot(globIdNode_j, 1);

                                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i),
                                                                                                          *fracture)),
                                            Psi_j = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_j),
                                                                                                           *fracture));

                                    double Psi = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                                    integralValue_ee_Ahd = (Psi - Psi_i) * (Psi - Psi_j) * (1.0 / 6.0) * abs_detJ * (nu * gradPhi_j).transpose() * JJ * gradPhi_i;

                                    AhD.coeffRef(ii_enr - 1, jj_enr - 1) += integralValue_ee_Ahd;

                                    // (FINE CASISTICA 9B) ===================================================================================================
                                }
                            }
                        }
                    }

                    // COSTRUZIONE RHS =======================================================================================================================

                    double averaged_f = this->physicalParameters->forcingTermAveragedOnTetrahedron(element.Vertices, *fracture);
                    double Psi_i = Utilities::heaviside(Utilities::signedDistanceFunction(blockMesh.Cell0DCoordinates(globIdNode_i), *fracture));
                    double Psi   = Utilities::heaviside(Utilities::signedDistanceFunction(tetrahedronBarycenter, *fracture));

                    rightHandSide(ii_enr - 1) += (Psi - Psi_i) * averaged_f * abs_detJ * (1.0 / 24.0);

                    // (FINE COSTRUZIONE RHS) =================================================================================================================


                }
            }
        }
    }

}


MappingFromReferenceTetrahedronInfo P1MatrixAssembler::MapFromReferenceTetrahedron(const Gedim::GeometryUtilities::Polyhedron element)
{
    MappingFromReferenceTetrahedronInfo result;

    double X0 = element.Vertices(0,0), X1 = element.Vertices(0,1), X2 = element.Vertices(0,2), X3 = element.Vertices(0,3),
           Y0 = element.Vertices(1,0), Y1 = element.Vertices(1,1), Y2 = element.Vertices(1,2), Y3 = element.Vertices(1,3),
           Z0 = element.Vertices(2,0), Z1 = element.Vertices(2,1), Z2 = element.Vertices(2,2), Z3 = element.Vertices(2,3);

    double detJ = X0*Y1*Z2 - X0*Y2*Z1 - X1*Y0*Z2 + X1*Y2*Z0 + X2*Y0*Z1 - X2*Y1*Z0 - X0*Y1*Z3 +
                  X0*Y3*Z1 + X1*Y0*Z3 - X1*Y3*Z0 - X3*Y0*Z1 + X3*Y1*Z0 + X0*Y2*Z3 - X0*Y3*Z2 -
                  X2*Y0*Z3 + X2*Y3*Z0 + X3*Y0*Z2 - X3*Y2*Z0 - X1*Y2*Z3 + X1*Y3*Z2 + X2*Y1*Z3 -
                  X2*Y3*Z1 - X3*Y1*Z2 + X3*Y2*Z1;

    Eigen::Matrix3d Jinv;

    Jinv << -(Y0*Z2-Y2*Z0-Y0*Z3+Y3*Z0+Y2*Z3-Y3*Z2), (X0*Z2-X2*Z0-X0*Z3+X3*Z0+X2*Z3-X3*Z2), -(X0*Y2-X2*Y0-X0*Y3+X3*Y0+X2*Y3-X3*Y2),
             (Y0*Z1-Y1*Z0-Y0*Z3+Y3*Z0+Y1*Z3-Y3*Z1), -(X0*Z1-X1*Z0-X0*Z3+X3*Z0+X1*Z3-X3*Z1), (X0*Y1-X1*Y0-X0*Y3+X3*Y0+X1*Y3-X3*Y1),
            -(Y0*Z1-Y1*Z0-Y0*Z2+Y2*Z0+Y1*Z2-Y2*Z1), (X0*Z1-X1*Z0-X0*Z2+X2*Z0+X1*Z2-X2*Z1), -(X0*Y1-X1*Y0-X0*Y2+X2*Y0+X1*Y2-X2*Y1);

    Jinv = (1.0 / detJ) * Jinv;

    Eigen::Matrix3d JJ = Jinv * Jinv.transpose();


    double abs_detJ = fabs(detJ);

    result.JJ = JJ;
    result.abs_detJ = abs_detJ;

    return result;

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

//                    if (jj > 0)

//                        this->constructElement_AhF(Ah, ii, jj, e);
//                    else

//                        this->constructElement_AhF(Ah_dirich, ii, -jj, e);

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
                rhs(ii - 1) += (1 / 3) * averaged_gN * triangle_area;
        }
    }
}

// **************************************************************************************









}
