#ifndef XFEMPOISSONASSEMBLER_HPP
#define XFEMPOISSONASSEMBLER_HPP

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include "Fracture3D.hpp"
#include "PhysicalParameters.hpp"

#include <map>

namespace XFEM_3D
{


enum struct integrationType
{
    std_std,
    std_enr,
    enr_std,
    enr_enr,
    std,  // for right hand side
    enr   // for right hand side
};

enum struct fractureBorder
{
    positive,
    negative,
    fracture
};

class P1MatrixAssembler
{
private:

    // 3D Mesh and pivot vector
    Gedim::MeshMatricesDAO* hD_Mesh;
    Eigen::MatrixXi* hD_Pivot;

    // Addtional boolean matrices for the enrichments and 2D/3D coupling
    Eigen::SparseMatrix<unsigned int>* toEnrich_elements;
    int num_enrichments = -1;

    Eigen::SparseMatrix<unsigned int>* hD_PsiP_MeshIntersections;
    Eigen::SparseMatrix<unsigned int>* hD_PsiM_MeshIntersections;
    Eigen::SparseMatrix<unsigned int>* hD_PsiF_MeshIntersections;

    // 2D Mesh and Pivot vector
    Gedim::MeshMatricesDAO* hF_Mesh;
    Eigen::VectorXi* hF_Pivot;

    // Fracture
    Fracture3D* fracture;

    // Physical parameters (permeabilities and normal trasmissivity of fracture)
    PhysicalParameters* physicalParameters;

    // Other utilities objects
    Gedim::GeometryUtilities* geometryUtilities;
    Gedim::MeshUtilities* meshUtilities;




public:

    P1MatrixAssembler(Fracture3D* fracture,
                        Gedim::GeometryUtilities* geometryUtilities,
                        Gedim::MeshUtilities* meshutilities);

    // Assemblers

    // This function is a public wrapper function for all the "assemble" functions
    void assembleXFEM();



    // Setters ******************************************************************************

     void setHD_Mesh(Gedim::MeshMatricesDAO *newHD_Mesh);
     void setHD_Pivot(Eigen::MatrixXi *newHD_Pivot);
     void setHF_Mesh(Gedim::MeshMatricesDAO *newHF_Mesh);
     void setHF_Pivot(Eigen::VectorXi *newHF_Pivot);
     void setPhysicalParameters(PhysicalParameters *newPhysicalParameters);

     // **************************************************************************************

     // Construction of _2D_3D_intersections and enrichment information

     void initialize(unsigned int numDOF, unsigned int numDirich);

private:

     // Right hand side *******************************************************************

     void constructElement_Rhs(Eigen::VectorXd&      rhs,
                               const unsigned int    i,
                               const unsigned int    ii_std,
                               const unsigned int    elementIndex,
                               const integrationType type);

     // PDE Discretization *****************************************************************
     // hD - hD
     void constructElement_AhD(Eigen::SparseMatrix<double>& M,
                               const unsigned int row,
                               const unsigned int col,
                               const unsigned int ii_std,
                               const unsigned int jj_std,
                               const unsigned int elementIndex,
                               const integrationType type);

     // h - h
     void constructElement_AhF(Eigen::SparseMatrix<double>& M,
                               const unsigned int i,
                               const unsigned int j,
                               const unsigned int elementIndex);

     // TODO...


     // Functional Discretization ************************************************************



     // TODO...

     // **************************************************************************************
public:

     // Funzioni assemble: queste saranno chiamate dentro assembleXFEM.
     void assemble_hD_hD(Eigen::SparseMatrix<double>& AhD,
                         Eigen::SparseMatrix<double>& AhD_dirich,
                         Eigen::SparseMatrix<double>& GhD,
                         Eigen::SparseMatrix<double>& GhD_dirich,
                         Eigen::VectorXd&             rightHandSide);

     void assemble_h_h(Eigen::SparseMatrix<double>& Ah,
                       Eigen::SparseMatrix<double>& Ah_dirich,
                       Eigen::SparseMatrix<double>& Gh,
                       Eigen::SparseMatrix<double>& Gh_dirich,
                       Eigen::SparseMatrix<double>& EF,
                       Eigen::SparseMatrix<double>& EF_dirich,
                       Eigen::SparseMatrix<double>& GPsiF,
                       Eigen::SparseMatrix<double>& GPsiF_dirich);

     void assemble_hD_h(Eigen::SparseMatrix<double>& BPsiF,
                        Eigen::SparseMatrix<double>& BPsiF_dirich);

     void assemble_hD_Psi(Eigen::SparseMatrix<double>& BPsiPlus,
                          Eigen::SparseMatrix<double>& BPsiPlus_dirich,
                          Eigen::SparseMatrix<double>& BPsiMinus,
                          Eigen::SparseMatrix<double>& BPsiMinus_dirich,
                          Eigen::SparseMatrix<double>& EPlus,
                          Eigen::SparseMatrix<double>& EPlus_dirich,
                          Eigen::SparseMatrix<double>& EMinus,
                          Eigen::SparseMatrix<double>& EMinus_dirich);

     void assemble_h_Psi(Eigen::SparseMatrix<double>& DPsiPlus,
                         Eigen::SparseMatrix<double>& DPsiPlus_dirich,
                         Eigen::SparseMatrix<double>& DPsiMinus,
                         Eigen::SparseMatrix<double>& DPsiMinus_dirich);

     void assemble_Psi_Psi(Eigen::SparseMatrix<double>& GPsiPlus,
                           Eigen::SparseMatrix<double>& GPsiPlus_dirich,
                           Eigen::SparseMatrix<double>& GPsiMinus,
                           Eigen::SparseMatrix<double>& GPsiMinus_dirich);


    // **************************************************************************************


     void initialize_2D_3DCoupling(fractureBorder type);

     void initializeEnrichmentInformation(unsigned int numDOF_3D, unsigned int numDirich_3D);

     inline unsigned int getNumberEnrichments() { return this->num_enrichments; };



};
}
#endif // XFEMPOISSONASSEMBLER_HPP
