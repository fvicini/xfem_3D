#include "Configurations.hpp"
#include "IOUtilities.hpp"

#include "EllipticProblem.hpp"

using namespace std;

// ***************************************************************************
int main(int argc, char** argv)
{
  // Register program configuration
  const XFEM_3D::EllipticProblem_ProgramConfiguration programConfiguration;

  // Import Parameters
  if (!Gedim::Output::FileExists("./Parameters.ini"))
    Gedim::Configurations::ExportToIni("./Parameters.ini",
                                       false);
  else
    Gedim::Configurations::InitializeFromIni("./Parameters.ini");
  Gedim::Configurations::Initialize(argc, argv);

  XFEM_3D::EllipticProblem program(programConfiguration);
  program.Run();

  // Close Program
  Gedim::Configurations::Reset();

  return 0;
}
// ***************************************************************************
