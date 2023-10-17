#include "P1MatrixAssembler.hpp"
#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "MapTetrahedron.hpp"
#include "Utilities.hpp"
#include <string>

#include "discontinousTestProblem_1.h"

#define DEBUG_CONSTANT 42
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


    // Vettore di bool. Numero di righe = numero di celle. Valore: cella da tagliare si/no.
    this->toEnrich_elements = new Eigen::SparseMatrix<unsigned int>(mesh3D->Cell3DTotalNumber(), 1);


    unsigned int global_counter_for_enriched_nodes = 0;


    for (unsigned int e = 0; e < mesh3D->Cell3DTotalNumber(); e++)
    {
        const Gedim::GeometryUtilities::Polyhedron elementAsPolyhedron = this->meshUtilities->MeshCell3DToPolyhedron(*mesh3D, e);

        double product = this->fracture->intersects(elementAsPolyhedron);
        bool fracture_intersects_element_not_too_close = fabs(product) >= GEOMETRIC_TOLERANCE && product < 0;


        if (fracture_intersects_element_not_too_close)
        {

            std::cout << "Elemento arricchito: " << e << std::endl;

            // Setto l'elemento 'e' come da arricchire.
            this->toEnrich_elements->insert(e, 0) = 1;

            /* OLD
            // Numerazione dei DOF arricchiti
            for (unsigned int n = 0; n <= 3; n++)
            {
                unsigned int glob_id_node = mesh3D->Cell3DVertex(e, n);

                int nn_enr;

                bool node_still_to_be_numbered = (*pivot3D)(glob_id_node, 1) < 0;

                if (node_still_to_be_numbered)
                {
                    global_counter_for_enriched_nodes++;
                    nn_enr = numDOF_3D_std + global_counter_for_enriched_nodes;
                    (*pivot3D)(glob_id_node, 1) = nn_enr;
                }
            }
            */

            // -----------------------------------------------------------------------
            // NEW:
            // Marko con '1' i i nodi da arricchire nella seconda colonna di pivot
            for (unsigned int n = 0; n <= 3; n++)
            {
                unsigned int glob_id_node = mesh3D->Cell3DVertex(e, n);

                bool node_still_to_be_numbered = (*pivot3D)(glob_id_node, 1) < 0;

                if (node_still_to_be_numbered)
                    (*pivot3D)(glob_id_node, 1) = 1;
            }
            // -----------------------------------------------------------------------
        }
    }
    // -----------------------------------------------------------------------
    // NEW:
    int counter = 0;
    for (unsigned int n = 0; n < mesh3D->Cell0DTotalNumber(); n++)
    {
        if((*pivot3D)(n, 1) == 1)
        {
            (*pivot3D)(n, 1) = numDOF_3D_std + counter + 1;
            counter++;
        }
    }
    // -----------------------------------------------------------------------
    this->num_enrichments = counter;
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
                                       Eigen::VectorXd& rightHandSide)
{
    Gedim::MeshMatricesDAO blockMesh = *this->hD_Mesh;
    Eigen::MatrixXi pivot = *this->hD_Pivot;

    Eigen::Matrix<double, 3, 8> GradPhi;
    GradPhi << -1, 1, 0, 0, -1, 1, 0, 0,
               -1, 0, 1, 0, -1, 0, 1, 0,
               -1, 0, 0, 1, -1, 0, 0, 1;


    Eigen::Matrix3d nu = this->physicalParameters->getPermeabilityTensorVolume();

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {
        Gedim::GeometryUtilities::Polyhedron element_as_polyhedron = meshUtilities->MeshCell3DToPolyhedron(blockMesh, e);

        bool element_cut_by_fracture = this->toEnrich_elements->coeff(e, 0);

        MappingFromReferenceTetrahedronInfo mapping = this->MapFromReferenceTetrahedron(element_as_polyhedron);
        Eigen::Matrix3d JJ = mapping.JJ;

        // Se il tetraedro è tagliato dalla frattura, dovrò integrare su una sua sotto-partizione. Per semplicità, costruisco una lista
        // sub_tetrahedra, nella quale, se l'elemento è da tagliare, ci saranno i sottotetraedro; altrimenti (se non è da tagliare), ci
        // sarà l'elemento stesso.
        std::list<Gedim::GeometryUtilities::Polyhedron> sub_tetrahedra;
        std::vector<Eigen::Matrix<double, 8, 8>> k_enr_list;
        if (element_cut_by_fracture)
        {
            std::vector<Gedim::GeometryUtilities::Polyhedron> tmp = Utilities::splitTetrahedronInSubTetrahedra(element_as_polyhedron,
                                                                                                               this->fracture,
                                                                                                               this->geometryUtilities,
                                                                                                               this->meshUtilities);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(sub_tetrahedra));
        }
        else
            sub_tetrahedra.push_back(element_as_polyhedron);


        // Determine the mask matrix k_enr (one for each sub-tetrahedron)
        for (auto& sub_t : sub_tetrahedra)
        {
            Eigen::Vector3d x_st = geometryUtilities->PolyhedronBarycenter(sub_t.Vertices);
            Eigen::Matrix<double, 8, 8> k_enr = Utilities::compute_k_enr(x_st, element_as_polyhedron,
                                                                         this->fracture);
            k_enr_list.push_back(k_enr);
        }


        // Local matrix A_e construction
        Eigen::Matrix<double, 8, 8> A_e;
        Eigen::Vector<double, 8> rhs_e;
        A_e << 0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0;

        rhs_e << 0,0,0,0,0,0,0,0;

        for (unsigned int i = 0; i <= 7; i++)
        {
            Eigen::Vector3d gradPhi_i = GradPhi.col(i);

            for (unsigned int j = 0; j <= 7; j++)
            {
                Eigen::Vector3d gradPhi_j = GradPhi.col(j);

                unsigned int st = 0;
                for (auto& sub_t : sub_tetrahedra)
                {
                    Eigen::Matrix<double, 4, 4> vol_matrix;

                    vol_matrix << sub_t.Vertices.col(0).x(), sub_t.Vertices.col(1).x(), sub_t.Vertices.col(2).x(), sub_t.Vertices.col(3).x(),
                                  sub_t.Vertices.col(0).y(), sub_t.Vertices.col(1).y(), sub_t.Vertices.col(2).y(), sub_t.Vertices.col(3).y(),
                                  sub_t.Vertices.col(0).z(), sub_t.Vertices.col(1).z(), sub_t.Vertices.col(2).z(), sub_t.Vertices.col(3).z(),
                                  1, 1, 1, 1;

                    double sub_t_volume = (1.0 / 6.0) * fabs(vol_matrix.determinant());

                    Eigen::Matrix<double, 8, 8> k_enr = k_enr_list.at(st);

                    // Matrix
                    A_e(i, j) += sub_t_volume * (nu * gradPhi_i).dot(JJ * gradPhi_j) * k_enr(i, j);

                    // Right hand side
                    double averaged_f = this->physicalParameters->forcingTermAveragedOnTetrahedron(sub_t.Vertices,
                                                                                                   *this->fracture);
                    rhs_e(i) += (1.0 / 4.0) * sub_t_volume * averaged_f * k_enr(i, 0);

                    st++;
                }
            }

        }

        // # Put A_e local entries into global stiffness matrix + construct right hand side.
        for (unsigned int i = 0; i <= 7; i++)
        {
            bool is_enr_i = i > 3;
            unsigned int id_i = blockMesh.Cell3DVertex(e, i % 4);
            int ii = pivot(id_i, is_enr_i);

            for (unsigned int j = 0; j <= 7; j++)
            {
                bool is_enr_j = j > 3;
                unsigned int id_j = blockMesh.Cell3DVertex(e, j % 4);
                int jj = pivot(id_j, is_enr_j);

                if (ii > 0)
                {
                    if (jj > 0)
                    {
                        AhD.coeffRef(ii-1, jj-1) += A_e(i, j);
                    }
                    else
                    {
                        //AhD_dirich.coeffRef(ii-1, -jj-1) = A_e(i, j);
                    }

                    rightHandSide(ii-1) = rhs_e(i);
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
