

#include "Fracture3D.hpp"
#include "MeshMatricesDAO.hpp"
#include "Utilities.hpp"


namespace XFEM_3D
{

namespace testZero
{



    void constructBlockDomainAndMesh();

    void setBoundaryConditions();

    Eigen::VectorXd exactSolution(Eigen::MatrixXd meshPointsCoordinates,
                                                    unsigned int numDofs3D_std,
                                                    Eigen::MatrixXi pivot,
                                                    Fracture3D* fracture);

    double forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture);

    double g_Neumann_triangle_averaged(std::vector<unsigned int> vertices_global_IDs,
                                       unsigned int neumannBorderMarker,
                                       Gedim::MeshMatricesDAO blockMesh,
                                       Fracture3D* fracture);


}
}

