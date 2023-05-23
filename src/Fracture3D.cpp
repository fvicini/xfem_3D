
#include "Fracture3D.hpp"
#include "GeometryUtilities.hpp"

using namespace Eigen;

namespace XFEM_3D {


Fracture3D::Fracture3D(Vector3d p1, Vector3d p2, Vector3d p3, Gedim::GeometryUtilities geometryUtilities)
{
    this->polygon     = geometryUtilities.CreateParallelogram(p1, p2, p3);
    this->normal      = geometryUtilities.PolygonNormal(this->polygon);
    this->translation = geometryUtilities.PolygonTranslation(this->polygon);
    this->rotation    = geometryUtilities.PolygonRotationMatrix(this->polygon,
                                                                this->normal,
                                                                this->translation);
};

double Fracture3D::intersects(Gedim::GeometryUtilities::Polyhedron element)
{
    Vector3d normal = this->normal;
    Vector3d x = this->getOrigin();
    double a = normal(0);
    double b = normal(1);
    double c = normal(2);
    double d = -a*x(0) - b*x(1) - c*x(2);
    std::vector<double> signed_distances;

    for (unsigned int i = 0; i <= 3; i++)
    {
        double sdf;
        Vector3d p = element.Vertices(seq(0,2), i);
        sdf = (a*p(0) + b*p(1) + c*p(2) + d) / (sqrt(a*a + b*b + c*c));
        signed_distances.push_back(sdf);
    }

    double max_sd = *max_element(signed_distances.begin(),
                                 signed_distances.end());
    double min_sd = *min_element(signed_distances.begin(),
                                 signed_distances.end());

//    if (max_sd * min_sd < 0)
//    {
//        // element is cut by the fracture
//        return true;
//    }

    return max_sd * min_sd;
}



}
