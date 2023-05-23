#ifndef FRACTURE3D_HPP
#define FRACTURE3D_HPP

#endif // FRACTURE3D_HPP
#include "GeometryUtilities.hpp"
#pragma once

namespace XFEM_3D {

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
    Fracture3D() {};
    ~Fracture3D() {};
    Fracture3D(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Gedim::GeometryUtilities);

    inline Eigen::MatrixXd getPolygon()     { return this->polygon; }
    inline Eigen::Vector3d getNormal()      { return this->normal; }
    inline Eigen::Vector3d getTranslation() { return this->translation; }
    inline Eigen::Matrix3d getRotation()    { return this->rotation; }

    inline Eigen::Vector3d getOrigin() { return this->polygon(Eigen::seq(0,2), 0); }
    double intersects(Gedim::GeometryUtilities::Polyhedron element);


};

}
