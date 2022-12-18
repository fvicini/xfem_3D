#ifndef PHYSICALPARAMETERS_H
#define PHYSICALPARAMETERS_H

#include "GeometryUtilities.hpp"

class PhysicalParameters
{
private:
    Eigen::Matrix3d permeabilityTensorVolume;
    Eigen::Matrix2d permeabilityTensorFracture;
    double normalTransmissivityFracture;
public:
    PhysicalParameters();

    Eigen::Matrix3d getPermeabilityTensorVolume() const;
    void setPermeabilityTensorVolume(const Eigen::Matrix3d &newPermeabilityTensorVolume);
    Eigen::Matrix2d getPermeabilityTensorFracture() const;
    void setPermeabilityTensorFracture(const Eigen::Matrix2d &newPermeabilityTensorFracture);

    void setIdrostaticPermeabilityTensorOnVolume(const double p);
    void setIdrostaticPermeabilityTensorOnFracture(const double p);
    void setNormalTransmissivityFracture(double newNormalTransmissivityFracture);
    double getNormalTransmissivityFracture() const;
};

#endif // PHYSICALPARAMETERS_H
