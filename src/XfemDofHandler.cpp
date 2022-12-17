#include "XfemDofHandler.hpp"

namespace XFEM_3D
{

XfemDofHandler::XfemDofHandler()
{
}

XfemDofHandler::XfemDofHandler(Gedim::MeshMatricesDAO& blockMesh)
{
    // Constructs boundary condition vector.
    // IDEA: vector with same size of point vector in blockMesh.
    // contains:   0 if internal DOF
    //             1 if boundary Dirichlet node
    //             2 if boundary DOF (Neumann node)


    for (unsigned int n = 0; n < blockMesh.Cell0DTotalNumber(); n++)
    {
        this->boundaryConditions[n] = 0;
        this->pivot[n] = blockMesh.Cell0DMarker(n);
    }

}

}
