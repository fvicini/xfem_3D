#ifndef __EllipticProblem_H
#define __EllipticProblem_H

#include "Configurations.hpp"

#include "GeometryUtilities.hpp"

#include <string>

namespace XFEM_3D
{

  class Fracture3D
  {
    // Attributes
    private:
      Eigen::MatrixXd polygon;
      Eigen::Vector3d normal;
      Eigen::Vector3d translation;
      Eigen::Matrix3d rotation;

    // Methods
    public:
      Fracture3D();
      ~Fracture3D();
      Fracture3D(Eigen::MatrixXd, Eigen::Vector3d, Eigen::Vector3d, Eigen::Matrix3d);

      inline Eigen::MatrixXd getPolygon()     { return this->polygon; }
      inline Eigen::Vector3d getNormal()      { return this->normal; }
      inline Eigen::Vector3d getTranslation() { return this->translation; }
      inline Eigen::Matrix3d getRotation()    { return this->rotation; }

      inline Eigen::Vector3d getOrigin() { return this->polygon(Eigen::seq(0,2), 0); }


  };


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
      inline double Pippo() const
      { return Gedim::Configurations::GetPropertyValue<int>("Pippo"); }
  };


  class EllipticProblem final
  {
    private:
      const EllipticProblem_ProgramConfiguration& config;

    public:
      EllipticProblem(const EllipticProblem_ProgramConfiguration& config);
      ~EllipticProblem();

      void Run();

    private:
      bool FractureItersectsElement(Fracture3D*, Gedim::GeometryUtilities::Polyhedron);
};


}

#endif
