#include "EllipticProblem.hpp"

#include "Configurations.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_CholeskySolver.hpp"
#include "math.h"

using namespace std;
using namespace Eigen;

namespace XFEM_3D
{
// ***************************************************************************
EllipticProblem_ProgramConfiguration::EllipticProblem_ProgramConfiguration()
{
    // Export parameters
    Gedim::Configurations::AddProperty("ExportFolder",
                                       "./Run",
                                       "Folder where to export data (Default: ./Export)");
    // Geometric parameters
    Gedim::Configurations::AddProperty("GeometricTolerance",
                                       numeric_limits<double>::epsilon(),
                                       "Geometric tolerance to perform 1D operations (Default: machine epsilon)");

    // Mesh parameters
    Gedim::Configurations::AddProperty("MeshMaximumTetrahedronVolume",
                                       0.1,
                                       "Mesh 3D maximum tetrahedron volume (Default: 0.1)");
}
// ***************************************************************************
EllipticProblem::EllipticProblem(const EllipticProblem_ProgramConfiguration& config) :
    config(config)
{
}
EllipticProblem::~EllipticProblem()
{
}
// ***************************************************************************
void EllipticProblem::Run()
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = config.GeometricTolerance();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    /// Create folders
    const string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const string exportCsvFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    const string logFolder = exportFolder + "/Log";

    /// Export Configuration of the following Run
    Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini",
                                       false);

    /// Set Log folder
    Gedim::Output::CreateFolder(logFolder);
    Gedim::LogFile::LogFolder = logFolder;

    /// Get problem data
    // Create domain
    const Gedim::GeometryUtilities::Polyhedron block = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    // Create fracture network
    const unsigned int num_fractures = 1;
    vector<Fracture3D*> fractureNetwork;

    const MatrixXd f3D = geometryUtilities.CreateParallelogram(Vector3d(-0.5, -0.5, 0.5),
                                                               Vector3d(2.0, 0.0, 0.0),
                                                               Vector3d(0.0, 2.0, 0.0));
    const Vector3d fNormal = geometryUtilities.PolygonNormal(f3D);
    const Vector3d fTransl = geometryUtilities.PolygonTranslation(f3D);
    const Matrix3d fRot = geometryUtilities.PolygonRotationMatrix(f3D,
                                                                  fNormal,
                                                                  fTransl);
    const MatrixXd fracture2D = geometryUtilities.RotatePointsFrom3DTo2D(f3D,
                                                                         fRot.transpose(),
                                                                         fTransl);

    fractureNetwork.push_back(new Fracture3D(f3D, fNormal, fTransl, fRot));


    // Export Fracture
    for (unsigned int f = 0; f < num_fractures; f++)
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(fractureNetwork.at(f)->getPolygon());
        vtkUtilities.Export(exportVtuFolder + "/Facture" + to_string(f) + ".vtu");
    }

    // Export block
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolyhedron(block.Vertices,
                                   block.Edges,
                                   block.Faces);
        vtkUtilities.Export(exportVtuFolder + "/Block.vtu");
    }


    // Create block mesh
    Gedim::Output::PrintGenericMessage("Create Block Mesh...", true);

    Gedim::MeshMatrices blockMeshData;
    Gedim::MeshMatricesDAO blockMesh(blockMeshData);

    meshUtilities.CreateTetrahedralMesh(block.Vertices,
                                        block.Edges,
                                        block.Faces,
                                        config.MeshMaximumTetrahedronVolume(),
                                        blockMesh);


    // Export the block mesh
    meshUtilities.ExportMeshToVTU(blockMesh,
                                  exportVtuFolder,
                                  "Block_Mesh");

    Gedim::Output::PrintStatusProgram("Create Block Mesh");

    Gedim::Output::PrintGenericMessage("Create Fracture Mesh...", true);

    Gedim::MeshMatrices fractureMeshData;
    Gedim::MeshMatricesDAO fractureMesh(fractureMeshData);

    meshUtilities.CreateTriangularMesh(fracture2D,
                                       0.1,
                                       fractureMesh);


    // Export the block mesh

    meshUtilities.ExportMeshToVTU(fractureMesh,
                                  exportVtuFolder,
                                  "Fracture_Mesh");

    Gedim::Output::PrintStatusProgram("Create Fracture Mesh");

    Gedim::Output::PrintGenericMessage("Compute block geometric properties...", true);

    Gedim::MeshUtilities::MeshGeometricData3D meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                        blockMesh);

    string exportTest = "./TEST";
    Gedim::Output::CreateFolder(exportTest);


    // Determine the effective number of DOFs (enrichment included).
    // Construction of boolean matrices
    double numDofs = blockMeshData.NumberCell0D; // Not so true: Dirichlet conditions?

    Gedim::Eigen_SparseArray<> toEnrich_elements;
    Gedim::Eigen_SparseArray<> toenrich_nodes;

    toEnrich_elements.Create();
    toenrich_nodes.Create();

    toEnrich_elements.SetSize(blockMesh.Cell3DTotalNumber(), num_fractures);
    toenrich_nodes.SetSize(blockMeshData.NumberCell0D, num_fractures);

    for (unsigned int e = 0; e < blockMesh.Cell3DTotalNumber(); e++)
    {
        const Gedim::GeometryUtilities::Polyhedron element = meshUtilities.MeshCell3DToPolyhedron(blockMesh, e);

        for (unsigned int f = 0; f < num_fractures; f++)
        {
            if (FractureItersectsElement(fractureNetwork.at(f), element))
            {
                // Set (e,f) entry of matrix to true, meaning element e is to enrich wrt fracture f.
                toEnrich_elements.Triplet(e, f, 1);

                // Retrieve the global Id of the verteces of element -> possible?




                // Export reproducing element to Paraview for visualization
                {
                    const vector<double> cellVolume(1, meshGeometricData.Cell3DsVolumes[e]);
                    Gedim::VTKUtilities vtkUtilities;
                    vtkUtilities.AddPolyhedron(element.Vertices,
                                               element.Edges,
                                               element.Faces
                                               /*{
                                                   {
                                                       "Volume",
                                                       Gedim::VTPProperty::Cells,
                                                       static_cast<unsigned int>(cellVolume.size()),
                                                       cellVolume.data()
                                                   }
                                               }*/);
                    vtkUtilities.Export(exportTest +
                                        "/ReproducingElement_" +
                                        to_string(e) +
                                        ".vtu");
                }
            }
            else // This is just to export the non-reproducing FEs to Paraview and compare them
            {
                {
                    const vector<double> cellVolume(1, meshGeometricData.Cell3DsVolumes[e]);
                    Gedim::VTKUtilities vtkUtilities;
                    vtkUtilities.AddPolyhedron(element.Vertices,
                                               element.Edges,
                                               element.Faces
                                               /*{
                                                   {
                                                       "Volume",
                                                       Gedim::VTPProperty::Cells,
                                                       static_cast<unsigned int>(cellVolume.size()),
                                                       cellVolume.data()
                                                   }
                                               }*/);
                    vtkUtilities.Export(exportTest +
                                        "/StandardElement_" +
                                        to_string(e) +
                                        ".vtu");
                }
            }

        }
    }
























    for (unsigned int c = 0; c < blockMesh.Cell3DTotalNumber(); c++)
    {
        const MatrixXd vertices = blockMesh.Cell3DVerticesCoordinates(c);

        const Gedim::GeometryUtilities::Polyhedron cellToPolyhedron = meshUtilities.MeshCell3DToPolyhedron(blockMesh,
                                                                                                           c);

//        {
//            const vector<double> cellVolume(1, meshGeometricData.Cell3DsVolumes[c]);
//            Gedim::VTKUtilities vtkUtilities;
//            vtkUtilities.AddPolyhedron(cellToPolyhedron.Vertices,
//                                       cellToPolyhedron.Edges,
//                                       cellToPolyhedron.Faces,
//                                       {
//                                           {
//                                               "Volume",
//                                               Gedim::VTPProperty::Cells,
//                                               static_cast<unsigned int>(cellVolume.size()),
//                                               cellVolume.data()
//                                           }
//                                       });
//            vtkUtilities.Export(exportTest +
//                                "/Cell_" +
//                                to_string(c) +
//                                ".vtu");
//        }

    }

    const double h = *max_element(std::begin(meshGeometricData.Cell3DsVolumes),
                                  std::end(meshGeometricData.Cell3DsVolumes));

    Gedim::Output::PrintStatusProgram("Compute block geometric properties");

    /// Assemble System
    Gedim::Output::PrintGenericMessage("Assemble System FEM...", true);


    Gedim::Eigen_SparseArray<> globalMatrixA;
    Gedim::Eigen_Array<> rightHandSide;
    Gedim::Eigen_Array<> solution;

    if (numDofs > 0)
    {
        MatrixXd cellMatrix = MatrixXd::Zero(3, 5);
        VectorXd cellVector = VectorXd::Zero(5);

        cellMatrix(0, 1) = 3.4;
        cellVector[1] = 3;
        cerr<< "cellMatrix\n"<< cellMatrix<< endl;
        cerr<< "cellVector\n"<< cellVector<< endl;
        cerr<< "result\n"<< cellMatrix * cellVector<< endl;

        globalMatrixA.SetSize(numDofs, numDofs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
        rightHandSide.SetSize(numDofs);
        solution.SetSize(numDofs);

        cerr<< rightHandSide<< endl;
        rightHandSide[0] = 10.4;
        cerr<< rightHandSide<< endl;
        rightHandSide.AddValue(0, 12.2);
        cerr<< rightHandSide<< endl;
        rightHandSide.Zeros();
        cerr<< rightHandSide<< endl;

        globalMatrixA.Triplet(1, 0, 12.4);
        globalMatrixA.Triplet(1, 0, 23.4);

        globalMatrixA.Create();
        cerr<< globalMatrixA<< endl;
        globalMatrixA.SetSize(numDofs, numDofs, Gedim::ISparseArray::SparseArrayTypes::Symmetric); // reset
        cerr<< globalMatrixA<< endl;
    }

    Gedim::Output::PrintStatusProgram("Assemble System");

    /// Solve
    Gedim::Output::PrintGenericMessage("Solve...", true);

    if (numDofs > 0)
    {
        Gedim::Eigen_CholeskySolver<> choleskySolver;
        choleskySolver.Initialize(globalMatrixA,
                                  rightHandSide,
                                  solution);
        choleskySolver.Solve();
        cerr<< solution<< endl;
    }

    Gedim::Output::PrintStatusProgram("Solve");
}


bool EllipticProblem::FractureItersectsElement(Fracture3D *fracture,
                                               Gedim::GeometryUtilities::Polyhedron element)
{

    Vector3d normal = fracture->getNormal();
    Vector3d x = fracture->getOrigin();
    double a = normal(0);
    double b = normal(1);
    double c = normal(2);
    double d = -a*x(0) - b*x(1) - c*x(2);
    vector<double> signed_distances;

    for (unsigned int i = 0; i <= 3; i++)
    {
        double sdf;
        Vector3d p = element.Vertices(seq(0,2), i);
        sdf = (a*p(0) + b*p(1) + c*p(2) + d) / (sqrt(a*a + b*b + c*c));
        signed_distances.push_back(sdf);
    }

    double max_sd = *max_element(signed_distances.begin(),
                                 signed_distances.end());
    double min_sd = *min_element(signed_distances.begin(),
                                 signed_distances.end());

    if (max_sd * min_sd < 0)
    {
        // element is cut by the fracture
        return true;
    }

    return false;

}

// *************************************************************************** Class Fracture3D

Fracture3D::Fracture3D(Eigen::MatrixXd pol,
                       Eigen::Vector3d norm,
                       Eigen::Vector3d transl,
                       Eigen::Matrix3d rot)
{
    this->polygon     = pol;
    this->normal      = norm;
    this->translation = transl;
    this->rotation    = rot;

}


}
