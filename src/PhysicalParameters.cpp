#include "PhysicalParameters.hpp"



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
