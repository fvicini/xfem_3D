#include "Fracture3D.hpp"
#include "MeshMatricesDAO.hpp"


namespace XFEM_3D
{


typedef struct {

    Eigen::MatrixXi pivot;
    std::vector<unsigned int> Neumann_triangles;

} setBoundaryConditionsResultType;


namespace testAllZero
{
    // TODO: implement
    void constructBlockDomainAndMesh();

    // TODO: implement
    setBoundaryConditionsResultType setBoundaryConditions(Gedim::MeshMatricesDAO& mesh);

    /// Returns 0 for all points
    Eigen::VectorXd exactSolution(Eigen::MatrixXd meshPointsCoordinates,
                                                   unsigned int numDofs3D_std,
                                                   Eigen::MatrixXi pivot,
                                                   Fracture3D* fracture);

    /// Returns null vector for all points
    Eigen::Vector3d exactSolutionGradient(Eigen::Vector3d point, Fracture3D* fracture);

    /// Returns 0 for all points
    double forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture);

    /// Returns 0 for all points
    double testAllZero_g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                       unsigned int neumannBorderMarker,
                                       Gedim::MeshMatricesDAO blockMesh,
                                       Fracture3D* fracture);


}
}

