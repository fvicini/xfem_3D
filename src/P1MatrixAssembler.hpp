#ifndef XFEMPOISSONASSEMBLER_HPP
#define XFEMPOISSONASSEMBLER_HPP

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshMatricesDAO.hpp"
#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include "Fracture3D.hpp"
#include "PhysicalParameters.hpp"

namespace XFEM_3D
{


enum struct integrationType
{
    std_std,
    std_enr,
    enr_std,
    enr_enr
};


class P1MatrixAssembler
{
private:

    // Mesh and pivot vector for the h^D variable
    Gedim::MeshMatricesDAO* hD_Mesh;
    Eigen::VectorXi* hD_Pivot;

    // Addtional boolean matrices for the enrichments and 2D/3D coupling
    Eigen::SparseMatrix<unsigned int>* toEnrich_nodes;
    Eigen::SparseMatrix<unsigned int>* toEnrich_elements;
    Eigen::SparseMatrix<unsigned int>* hD_PsiP_MeshIntersections;
    Eigen::SparseMatrix<unsigned int>* hD_PsiM_MeshIntersections;
    Eigen::SparseMatrix<unsigned int>* hD_PsiF_MeshIntersections;

    // Mesh and pivot vector for the h^F variable
    Gedim::MeshMatricesDAO* hF_Mesh;
    Eigen::VectorXi* hF_Pivot;

    // Mesh and pivot vector for the psi^p variable
    Gedim::MeshMatricesDAO* psiP_Mesh;
    Eigen::VectorXi* psiP_Pivot;

    // Mesh and pivot vector for the psi^m variable
    Gedim::MeshMatricesDAO* psiM_Mesh;
    Eigen::VectorXi* psiM_Pivot;

    // Mesh and pivot vector for the psi^F variable
    Gedim::MeshMatricesDAO* psiF_Mesh;
    Eigen::VectorXi* psiF_Pivot;

    // Mesh and pivot vector for the lambda^D Lagrange multiplier
    Gedim::MeshMatricesDAO* lambdaD_Mesh;
    Eigen::VectorXi* lambdaD_Pivot;

    // Mesh and pivot vector for the lambda^F Lagrange multiplier
    Gedim::MeshMatricesDAO* lambdaF_Mesh;
    Eigen::VectorXi* lambdaF_Pivot;

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

     void assemble_AhD(Gedim::Eigen_SparseArray<>& AhD,
                            Gedim::Eigen_SparseArray<>& AhD_dirich,
                            Gedim::Eigen_Array<>& rightHandSide);

     void assemble_AhF(Gedim::Eigen_SparseArray<>& AhF,
                            Gedim::Eigen_SparseArray<>& AhF_dirich,
                            Gedim::Eigen_Array<>& rightHandSide);

     void setHD_Mesh(Gedim::MeshMatricesDAO *newHD_Mesh);
     void setHD_Pivot(Eigen::VectorXi *newHD_Pivot);
     void setToEnrich_nodes(Eigen::SparseMatrix<unsigned int> *newToEnrich_nodes);
     void setToEnrich_elements(Eigen::SparseMatrix<unsigned int> *newToEnrich_elements);
     void setHF_Mesh(Gedim::MeshMatricesDAO *newHF_Mesh);
     void setHF_Pivot(Eigen::VectorXi *newHF_Pivot);
     void setPsiP_Mesh(Gedim::MeshMatricesDAO *newPsiP_Mesh);
     void setPsiP_Pivot(Eigen::VectorXi *newPsiP_Pivot);
     void setPsiM_Mesh(Gedim::MeshMatricesDAO *newPsiM_Mesh);
     void setPsiM_Pivot(Eigen::VectorXi *newPsiM_Pivot);
     void setPsiF_Mesh(Gedim::MeshMatricesDAO *newPsiF_Mesh);
     void setPsiF_Pivot(Eigen::VectorXi *newPsiF_Pivot);
     void setLambdaD_Mesh(Gedim::MeshMatricesDAO *newLambdaD_Mesh);
     void setLambdaD_Pivot(Eigen::VectorXi *newLambdaD_Pivot);
     void setLambdaF_Mesh(Gedim::MeshMatricesDAO *newLambdaF_Mesh);
     void setLambdaF_Pivot(Eigen::VectorXi *newLambdaF_Pivot);
     void setPhysicalParameters(PhysicalParameters *newPhysicalParameters);

     // Construction of _2D_3D_intersections.
     void initialize_2D_3DCoupling();

private:

    void constructElement_AhD(Gedim::Eigen_SparseArray<>& M,
                     const unsigned int i,
                     const unsigned int j,
                     const unsigned int elementIndex,
                     const integrationType type);

    void constructElement_AhF(Gedim::Eigen_SparseArray<>& M,
                     const unsigned int i,
                     const unsigned int j,
                     const unsigned int elementIndex);





    void constructEl_RHS_A(Gedim::Eigen_Array<> &rhs,
                           unsigned int i,
                           const Gedim::GeometryUtilities &geometryUtilities);
    void constructEl_RHS_B(Gedim::Eigen_Array<> &rhs,
                           unsigned int i,
                           const Gedim::GeometryUtilities &geometryUtilities);



};
}
#endif // XFEMPOISSONASSEMBLER_HPP
