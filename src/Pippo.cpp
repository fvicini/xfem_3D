#include "EllipticProblem.hpp"

#include "Configurations.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "Eigen_SparseArray.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_CholeskySolver.hpp"

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
    const Gedim::GeometryUtilities::Polyhedron block = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                              1.0);
    const Eigen::MatrixXd fracture3D = geometryUtilities.CreateParallelogram(Eigen::Vector3d(-0.5, -0.5, 0.5),
                                                                             Eigen::Vector3d(2.0, 0.0, 0.0),
                                                                             Eigen::Vector3d(0.0, 2.0, 0.0));
    const Eigen::Vector3d fractureNormal = geometryUtilities.PolygonNormal(fracture3D);
    const Eigen::Vector3d fractureTranslation = geometryUtilities.PolygonTranslation(fracture3D);
    const Eigen::Matrix3d fractureRotationMatrix = geometryUtilities.PolygonRotationMatrix(fracture3D,
                                                                                           fractureNormal,
                                                                                           fractureTranslation);
    const Eigen::MatrixXd fracture2D = geometryUtilities.RotatePointsFrom3DTo2D(fracture3D,
                                                                                fractureRotationMatrix.transpose(),
                                                                                fractureTranslation);

    // Export block
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolyhedron(block.Vertices,
                                   block.Edges,
                                   block.Faces);
        vtkUtilities.Export(exportVtuFolder + "/Block.vtu");
    }

    // Export fracture3D
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(fracture3D);
        vtkUtilities.Export(exportVtuFolder + "/Facture.vtu");
    }

    /// Create block mesh
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

    for (unsigned int c = 0; c < blockMesh.Cell3DTotalNumber(); c++)
    {
        const Eigen::MatrixXd vertices = blockMesh.Cell3DVerticesCoordinates(c);

        const Gedim::GeometryUtilities::Polyhedron cellToPolyhedron = meshUtilities.MeshCell3DToPolyhedron(blockMesh,
                                                                                                           c);

        {
            const vector<double> cellVolume(1, meshGeometricData.Cell3DsVolumes[c]);
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolyhedron(cellToPolyhedron.Vertices,
                                       cellToPolyhedron.Edges,
                                       cellToPolyhedron.Faces,
                                       {
                                           {
                                               "Volume",
                                               Gedim::VTPProperty::Cells,
                                               static_cast<unsigned int>(cellVolume.size()),
                                               cellVolume.data()
                                           }
                                       });
            vtkUtilities.Export(exportTest +
                                "/Cell_" +
                                to_string(c) +
                                ".vtu");
        }

    }

    const double h = *max_element(std::begin(meshGeometricData.Cell3DsVolumes),
                                  std::end(meshGeometricData.Cell3DsVolumes));

    Gedim::Output::PrintStatusProgram("Compute block geometric properties");

    /// Assemble System
    Gedim::Output::PrintGenericMessage("Assemble System FEM...", true);

    double numDofs = 3;

    Gedim::Eigen_SparseArray<> globalMatrixA;
    Gedim::Eigen_Array<> rightHandSide;
    Gedim::Eigen_Array<> solution;

    if (numDofs > 0)
    {
        Eigen::MatrixXd cellMatrix = Eigen::MatrixXd::Zero(3, 5);
        Eigen::VectorXd cellVector = Eigen::VectorXd::Zero(5);

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
// ***************************************************************************
}
