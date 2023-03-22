#include "Configurations.hpp"
#include "IOUtilities.hpp"

#include "EllipticProblem.hpp"

#include <list>

using namespace std;

typedef struct
{
    double h;
    double error;
} error_result;

// ***************************************************************************
int main(int argc, char** argv)
{
  // Register program configuration
  const XFEM_3D::EllipticProblem_ProgramConfiguration discontinousTestProblem_configuration;

  // Import Parameters
  if (!Gedim::Output::FileExists("./Parameters.ini"))
    Gedim::Configurations::ExportToIni("./Parameters.ini",
                                       false);
  else
    Gedim::Configurations::InitializeFromIni("./Parameters.ini");
  Gedim::Configurations::Initialize(argc, argv);


  // ************************************************************************************************
  // Problema test discontinuo con il quale testiamo l'andamento dell'errore previsto dalla teoria.

  double max_volumes[] = {0.05, 0.007, 0.001};
  list<error_result> results;

  for (double v : max_volumes)
  {
      XFEM_3D::EllipticProblem discontinousTestProblem(discontinousTestProblem_configuration);
      XFEM_3D::result_for_error_estimate XFEM_result =  discontinousTestProblem.Run(v);

      error_result res;

      string exportResultsFolder = "/home/Scrivania/tesi_matlab_visual/errore_VS_parametroMesh";
      string resultFile = exportResultsFolder + "/data.txt";
      ofstream f(resultFile);

      Eigen::VectorXd error = XFEM_result.exactSolution - XFEM_result.solution;
      res.error = error.norm();
      res.h = XFEM_result.h;

      if (f.is_open())
      {
          f << XFEM_result.h << " " << res.error << "n";

          f.close();

      } else cout << "Problem with opening file";

      results.push_back(res);

  }

  // Close Program
  Gedim::Configurations::Reset();
  // ************************************************************************************************


  return 0;
}
// ***************************************************************************
