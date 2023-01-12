#ifndef PHYSICALPARAMETERS_H
#define PHYSICALPARAMETERS_H

#include "GeometryUtilities.hpp"
#include "Utilities.hpp"
#include "Fracture3D.hpp"

namespace XFEM_3D {



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

    inline double forcingTerm(const Eigen::Vector3d pointCoords, Fracture3D fracture)
    {
        double x = pointCoords[0], y = pointCoords[1], z = pointCoords[2];

        if (Utilities::signedDistanceFunction(pointCoords, fracture) >= 0)
            return -32*( x*y*(1-x)*(1-y) + x*z*(1-x)*(1-z) + y*z*(1-y)*(1-z) );

        else
            return  32*( x*y*(1-x)*(1-y) + x*z*(1-x)*(1-z) + y*z*(1-y)*(1-z) );

    };
};

}
#endif // PHYSICALPARAMETERS_H
