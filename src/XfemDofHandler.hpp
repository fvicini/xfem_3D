#ifndef XFEMDOFHANDLER_HPP
#define XFEMDOFHANDLER_HPP

#include "MeshMatricesDAO.hpp"

namespace XFEM_3D
{

class XfemDofHandler
{
public:
    Eigen::VectorXd boundaryConditions;
    Eigen::VectorXd pivot;


public:
    XfemDofHandler();
    XfemDofHandler(Gedim::MeshMatricesDAO& mesh);
};

}

#endif // XFEMDOFHANDLER_HPP
