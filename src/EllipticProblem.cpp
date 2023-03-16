#include "EllipticProblem.hpp"
#include "Fracture3D.hpp"
#include "P1MatrixAssembler.hpp"

#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCholesky>

#include "Configurations.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
//#include "Eigen_SparseArray.hpp"
//#include "Eigen_Array.hpp"
//#include "Eigen_LUSolver.hpp"
//#include "Eigen_CholeskySolver.hpp"
#include "math.h"
#include "discontinousTestProblem_1.h"

#include "Quadrature_Gauss3D_Tetrahedron.hpp"
#include "MapTetrahedron.hpp"

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

void EllipticProblem::TestQuadrature(const Gedim::GeometryUtilities &geometryUtilities) const
{
    string exportFolder = "./TEST_QUAD";
    Gedim::Output::CreateFolder(exportFolder);

    Eigen::MatrixXd vertices;
    vertices.setZero(3, 4);
    vertices.row(0) << +5.0, +0.0, -8.0, +0.0;
    vertices.row(1) << -5.0, +8.0, -4.0, +0.0;
    vertices.row(2) << +0.0, +0.0, +3.0, +12.0;

    {
        const Gedim::GeometryUtilities::Polyhedron tet = geometryUtilities.CreateTetrahedronWithVertices(vertices.col(0),
                                                                                                         vertices.col(1),
                                                                                                         vertices.col(2),
                                                                                                         vertices.col(3));
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolyhedron(tet.Vertices,
                                   tet.Edges,
                                   tet.Faces
                                   );
        vtkUtilities.Export(exportFolder +
                            "/TET_" +
                            ".vtu");
    }

    {
        const Gedim::GeometryUtilities::Polyhedron tet = geometryUtilities.CreateTetrahedronWithVertices(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                         Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                         Eigen::Vector3d(0.0, 1.0, 0.0),
                                                                                                         Eigen::Vector3d(0.0, 0.0, 1.0));
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolyhedron(tet.Vertices,
                                   tet.Edges,
                                   tet.Faces
                                   );
        vtkUtilities.Export(exportFolder +
                            "/TET_REF_" +
                            ".vtu");
    }


    const unsigned int quadratureOrder = 2;
    Eigen::MatrixXd quadraturePointsRef;
    Eigen::VectorXd quadratureWeightsRef;
    Gedim::Quadrature_Gauss3D_Tetrahedron::FillPointsAndWeights(quadratureOrder,
                                                                quadraturePointsRef,
                                                                quadratureWeightsRef);

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPoints(quadraturePointsRef);
        vtkUtilities.Export(exportFolder +
                            "/POINT_REF_" +
                            ".vtu");
    }

    Gedim::MapTetrahedron mapping(geometryUtilities);
    const Gedim::MapTetrahedron::MapTetrahedronData mapData = mapping.Compute(vertices);

    Eigen::MatrixXd mappedPoints = mapping.F(mapData,
                                             quadraturePointsRef);


    Eigen::VectorXd mappedWeights = quadratureWeightsRef.array() * mapping.DetJ(mapData,
                                                                                quadraturePointsRef).array().abs();

    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPoints(mappedPoints);
        vtkUtilities.Export(exportFolder +
                            "/POINT_" +
                            ".vtu");
    }

    const double volume = mappedWeights.sum();
    double int_x_sqrd = 0.0;
    for (unsigned int q = 0; q < quadraturePointsRef.cols(); q++)
      int_x_sqrd += mappedPoints.col(q).x() * mappedPoints.col(q).y() * mappedWeights[q];

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

result_for_error_estimate EllipticProblem::Run(double max_volume_tetrahedra)
{

    // CONFIGURAZIONE INIZIALE *****************************************************************************

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = config.GeometricTolerance();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    //TestQuadrature(geometryUtilities);

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

    // Create fracture network (here only one is appended to the vector)
    const unsigned int num_fractures = 1;
    vector<Fracture3D*> fractureNetwork;
    vector<MatrixXd>    fractureNetwork2D;

    fractureNetwork.push_back(new Fracture3D(Vector3d(-0.5, -0.5, 0.5),
                                             Vector3d(2.0, 0.0, 0.0),
                                             Vector3d(0.0, 2.0, 0.0),
                                             geometryUtilities));


    // For every fracture, create the corresponding 2D domain.
    for (unsigned int f = 0; f < num_fractures; f++)
    {
        const MatrixXd fracture2D = geometryUtilities.RotatePointsFrom3DTo2D(fractureNetwork.at(f)->getPolygon(),
                                                                             fractureNetwork.at(f)->getRotation().transpose(),
                                                                             fractureNetwork.at(f)->getTranslation());
        fractureNetwork2D.push_back(fracture2D);
    }


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
                                        max_volume_tetrahedra,
                                        blockMesh,
                                        "Qpqfenza");

// al posto di max_volume_tetrahedra: config.MeshMaximumTetrahedronVolume().

    // Export the block mesh
    meshUtilities.ExportMeshToVTU(blockMesh,
                                  exportVtuFolder,
                                  "Block_Mesh");


    Gedim::Output::PrintStatusProgram("Create Block Mesh");


    // Create mesh for every 2D fracture in fractureNetwork2D
    Gedim::Output::PrintGenericMessage("Create Fracture Mesh...", true);

    Gedim::MeshMatrices fractureMeshData;
    Gedim::MeshMatricesDAO fractureMesh(fractureMeshData);

    meshUtilities.CreateTriangularMesh(fractureNetwork2D.at(0),
                                        0.1,
                                        fractureMesh);


    // Export the block mesh

    //    meshUtilities.ExportMeshToVTU(fractureMesh,
    //                                  exportVtuFolder,
    //                                  "Fracture_Mesh");

    //    Gedim::Output::PrintStatusProgram("Create Fracture Mesh");

    Gedim::Output::PrintGenericMessage("Compute block geometric properties...", true);

    Gedim::MeshUtilities::MeshGeometricData3D meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                        blockMesh);

    string exportTest = "./TEST";
    Gedim::Output::CreateFolder(exportTest);

    // ***************************************************************************************************



    // COSTRUZIONE DELLA MATRICE PIVOT PER LA MESH PER LA VARIABILE hD ***********************************
    /*
     * La prima colonna di pivot ('0') contiene la numerazione 1,2,3,4,5... (se DOF) e -1,-2,-3... (se Dirichlet).
     * La seconda colonna ('1') per il momento è tutta a -1. In fase di inizializzazione, l'assembler provvederà a riempirla
     * con l'indice del DOF di arricchimento corrispondente al dato DOF standard contenuto alla stessa riga, prima colonna.
     * */

    Eigen::MatrixXi pivot(blockMesh.Cell0DTotalNumber(), 2);
    std::vector<unsigned int> Neumann_triangles;
    unsigned int num_Dirichlet_3D = 0;
    unsigned int numDOF_3D_std = 0;
    for (unsigned int glob_id_point = 0; glob_id_point < blockMesh.Cell0DTotalNumber(); glob_id_point++)
    {
        Eigen::Vector3d point = blockMesh.Cell0DCoordinates(glob_id_point);
        double x = point.x(),
               y = point.y(),
               z = point.z();

        bool point_is_dirichlet = (z == 0.0) || (z == 1.0);

        if (point_is_dirichlet)
        {
            blockMesh.Cell0DSetMarker(glob_id_point, 2 * num_Dirichlet_3D + 1);
            num_Dirichlet_3D++;
            pivot(glob_id_point, 0) = -num_Dirichlet_3D;
            pivot(glob_id_point, 1) = -1;
        }
        else
        {
            blockMesh.Cell0DSetMarker(glob_id_point, 0);
            numDOF_3D_std++;
            pivot(glob_id_point, 0) = numDOF_3D_std;
            pivot(glob_id_point, 1) = -1;
        }
    }

    // Costruzione della struttura dati per i triangoli di Neumann
    for (unsigned int glob_id_triangle; glob_id_triangle < blockMesh.Cell2DTotalNumber(); glob_id_triangle++)
    {
        blockMesh.Cell2DSetMarker(glob_id_triangle, -1);

        std::vector<unsigned int> vertices_global_IDs = blockMesh.Cell2DVertices(glob_id_triangle);

        Eigen::Vector3d point1 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(0)),
                        point2 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(1)),
                        point3 = blockMesh.Cell0DCoordinates(vertices_global_IDs.at(2));

        double x1 = point1.x(), x2 = point2.x(), x3 = point3.x(),
               y1 = point1.y(), y2 = point2.y(), y3 = point3.y();

        bool triangle_lies_on_cube_face_with_marker_2  = (x1 == 0.0 && x2 == 0.0 && x3 == 0.0);
        bool triangle_lies_on_cube_face_with_marker_4  = (x1 == 1.0 && x2 == 1.0 && x3 == 1.0);
        bool triangle_lies_on_cube_face_with_marker_6  = (y1 == 0.0 && y2 == 0.0 && y3 == 0.0);
        bool triangle_lies_on_cube_face_with_marker_8  = (y1 == 1.0 && y2 == 1.0 && y3 == 1.0);


        // Faccia di Neumann con marker 2: {x=0}
        if (triangle_lies_on_cube_face_with_marker_2)
        {
            blockMesh.Cell2DSetMarker(glob_id_triangle, 2);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 4: {x=1}
        if (triangle_lies_on_cube_face_with_marker_4)
        {
            blockMesh.Cell2DSetMarker(glob_id_triangle, 4);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 6: {y=0}
        if (triangle_lies_on_cube_face_with_marker_6)
        {
            blockMesh.Cell2DSetMarker(glob_id_triangle, 6);
            Neumann_triangles.push_back(glob_id_triangle);
        }

        // Faccia di Neumann con marker 8: {y=1}
        if (triangle_lies_on_cube_face_with_marker_8)
        {
            blockMesh.Cell2DSetMarker(glob_id_triangle, 8);
            Neumann_triangles.push_back(glob_id_triangle);
        }
    }
    // ***************************************************************************************************

    // COSTRUZIONE DEL VETTORE PIVOT PER LA MESH PER LA VARIABILE h **************************************
    Eigen::VectorXi fracturePivot(fractureMesh.Cell0DTotalNumber());
    for (unsigned int n = 0; n < fractureMesh.Cell0DTotalNumber(); n++)
    {
        fracturePivot[n] = fractureMesh.Cell0DMarker(n);
    }
    // ***************************************************************************************************


    const double h = *max_element(std::begin(meshGeometricData.Cell3DsVolumes),
                                  std::end(meshGeometricData.Cell3DsVolumes));


    Gedim::Output::PrintStatusProgram("Compute block geometric properties");



    /// Assemble System
    Gedim::Output::PrintGenericMessage("Assemble System FEM...", true);


    // Physical parameters of the problem
    PhysicalParameters* params = new PhysicalParameters();
    params->setIdrostaticPermeabilityTensorOnVolume(1);
    params->setIdrostaticPermeabilityTensorOnFracture(1);
    params->setNormalTransmissivityFracture(1);

    // Creation and setting of the assembler object
    P1MatrixAssembler *assembler = new P1MatrixAssembler(fractureNetwork.at(0),
                                                         &geometryUtilities,
                                                         &meshUtilities);
    assembler->setHD_Mesh(&blockMesh);
    assembler->setHF_Mesh(&fractureMesh);
    assembler->setHD_Pivot(&pivot);
    assembler->setHD_NeumannInfo(&Neumann_triangles);
    assembler->setHF_Pivot(&fracturePivot);
    assembler->setPhysicalParameters(params);
    assembler->initialize(numDOF_3D_std, num_Dirichlet_3D);

    // Determine the effective number of DOFs (enrichment included).
    unsigned int numEnrichements = assembler->getNumberEnrichments();
    assert(numEnrichements != -1);

    unsigned int numDofs3D = blockMesh.Cell0DTotalNumber();
    numDofs3D -= num_Dirichlet_3D;
    numDofs3D += numEnrichements;

    // System matrices definition
    Eigen::SparseMatrix<double> AhD(numDofs3D, numDofs3D),
                                AhD_dirich(numDofs3D, num_Dirichlet_3D),
                                GhD,
                                GhD_dirich;
    Eigen::VectorXd rightHandSide(numDofs3D);
    Eigen::VectorXd solution(numDofs3D), exactSolution;

    // Per evitare che ci siano valori inizializzati a 'nan'.
    rightHandSide.setZero();

    Gedim::Output::PrintStatusProgram("Assemble System");

    if (numDofs3D > 0)
    {
        assembler->assemble_hD_hD(AhD, AhD_dirich, GhD, GhD_dirich, rightHandSide);
        assembler->addNeumann(rightHandSide);
    }


    // EXPORT AHD MATRIX TO .txt FILE and SAVE IN TRIPLETS FORM *********************************
    const string exportMatrFolder = exportFolder + "/Matrices";
    Gedim::Output::CreateFolder(exportMatrFolder);

    const string matrixFile = exportMatrFolder + "/AhD.txt";
    ofstream fw(matrixFile, std::ofstream::out);

    if (fw.is_open())
    {
        fw << "\n\n";

        Eigen::MatrixXd AhD_dense = AhD.toDense();
        for(unsigned int i = 0; i < AhD_dense.rows(); i++)
        {
            for(unsigned int j = 0; j < AhD_dense.cols(); j++)
            {
                fw << AhD_dense(i,j) << " ";
            }

            fw << "\n";
        }

        fw.close();
    }
    else std::cout << "Problem with opening file";

    // ***************************************************************************************



    // SOLUZIONE DEL SISTEMA LINEARE *********************************************************

    if (numDofs3D > 0)
    {
        SimplicialLLT<Eigen::SparseMatrix<double>> choleskySolver;

        choleskySolver.compute(AhD);

        if(choleskySolver.info() != Success)
          throw std::runtime_error("The Cholesky decomposition has failed.");

        solution = choleskySolver.solve(rightHandSide);

        exactSolution = DiscontinousTestProblem_1::exactSolution(blockMesh.Cell0DsCoordinates(),
                                                                    numDOF_3D_std,
                                                                    pivot,
                                                                    fractureNetwork.at(0));

        if(choleskySolver.info() != Success)
          throw std::runtime_error("The solution of the linear system has failed.");

    }
    // ***************************************************************************************

    Gedim::Output::PrintStatusProgram("Solve");

    result_for_error_estimate result;
    result.solution = solution(seq(0, numDOF_3D_std - 1));
    result.exactSolution = exactSolution;
    result.h = h;

    return result;
}

}
