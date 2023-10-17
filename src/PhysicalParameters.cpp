#include "PhysicalParameters.hpp"
#include "discontinousTestProblem_1.h"

namespace XFEM_3D {

PhysicalParameters::PhysicalParameters()
{

}

Eigen::Matrix3d PhysicalParameters::getPermeabilityTensorVolume() const
    {
        return permeabilityTensorVolume;
    }

void PhysicalParameters::setPermeabilityTensorVolume(const Eigen::Matrix3d &newPermeabilityTensorVolume)
    {
        permeabilityTensorVolume = newPermeabilityTensorVolume;
    }

Eigen::Matrix2d PhysicalParameters::getPermeabilityTensorFracture() const
    {
        return permeabilityTensorFracture;
    }

void PhysicalParameters::setPermeabilityTensorFracture(const Eigen::Matrix2d &newPermeabilityTensorFracture)
    {
        permeabilityTensorFracture = newPermeabilityTensorFracture;
    }

void PhysicalParameters::setIdrostaticPermeabilityTensorOnVolume(const double p)
    {
        this->permeabilityTensorVolume << p, 0, 0,
                                          0, p, 0,
                                          0, 0, p;
    }

void PhysicalParameters::setIdrostaticPermeabilityTensorOnFracture(const double p)
    {
        this->permeabilityTensorFracture << p, 0,
                                            0, p;
    }

void PhysicalParameters::setNormalTransmissivityFracture(double newNormalTransmissivityFracture)
    {
        normalTransmissivityFracture = newNormalTransmissivityFracture;
    }

double PhysicalParameters::getNormalTransmissivityFracture() const
    {
        return normalTransmissivityFracture;
    }

double PhysicalParameters::forcingTerm(Eigen::Vector3d pointCoords, Fracture3D fracture)
    {
        return DiscontinousTestProblem_1::forcingTerm(pointCoords, fracture);
    }

double PhysicalParameters::forcingTermAveragedOnTetrahedron(Eigen::MatrixXd tetrahedronVertices, Fracture3D fracture)
{
    double res = 0.0;
    for(unsigned int i = 0; i < 3; i++)
    {
        res += this->forcingTerm(tetrahedronVertices.col(i), fracture);
    }
    return res / 3;
}

}
