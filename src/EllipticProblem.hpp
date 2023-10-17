#ifndef __EllipticProblem_H
#define __EllipticProblem_H

#include "Configurations.hpp"

#include "GeometryUtilities.hpp"

#include <string>

namespace XFEM_3D
{

  typedef struct {

    Eigen::VectorXd solution;
    Eigen::VectorXd exactSolution;
    double h;

  } result_for_error_estimate;


  class EllipticProblem_ProgramConfiguration final
  {
    public:
      EllipticProblem_ProgramConfiguration();

      inline std::string ExportFolder() const
      { return Gedim::Configurations::GetPropertyValue<string>("ExportFolder"); }
      inline double GeometricTolerance() const
      { return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance"); }
      inline double MeshMaximumTetrahedronVolume() const
      { return Gedim::Configurations::GetPropertyValue<double>("MeshMaximumTetrahedronVolume"); }
  };


  class EllipticProblem final
  {
    private:
      const EllipticProblem_ProgramConfiguration& config;

      void TestQuadrature(const Gedim::GeometryUtilities& geometryUtilities) const;

    public:
      EllipticProblem(const EllipticProblem_ProgramConfiguration& config);
      ~EllipticProblem();


      result_for_error_estimate Run(double max_volume_tetrahedra);

      void RunTestXFEM();

};


}

#endif
