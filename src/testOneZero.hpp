

#include "Fracture3D.hpp"
#include "MeshMatricesDAO.hpp"
#include "Utilities.hpp"


namespace XFEM_3D
{

namespace testOneZero
{
    // TODO: implement
    void constructBlockDomainAndMesh();

    // TODO: implement
    void setBoundaryConditions();

    /// Returns 0, 1 depending on which side of fracture
    Eigen::VectorXd exactSolution(Eigen::MatrixXd meshPointsCoordinates,
                                                    unsigned int numDofs3D_std,
                                                    Eigen::MatrixXi pivot,
                                                    Fracture3D* fracture);

    /// Returns 0, 1 depending on which side of fracture
    Eigen::Vector3d exactSolutionGradient(Eigen::Vector3d point, Fracture3D* fracture);

    /// Returns 0, 1 depending on which side of fracture
    double forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture);

    /// Returns 0, 1 depending on which side of fracture
    double testAllZero_g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                       unsigned int neumannBorderMarker,
                                       Gedim::MeshMatricesDAO blockMesh,
                                       Fracture3D* fracture);


}
}

