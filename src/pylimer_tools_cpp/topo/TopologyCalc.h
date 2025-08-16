#ifndef TOPOLOGY_CALC_H
#define TOPOLOGY_CALC_H

#include <Eigen/Dense>

namespace pylimer_tools {
namespace topo {

  static inline Eigen::Vector3d sampleIntersectionPoint(
    const Eigen::Vector3d origin1,
    const double radius1,
    const Eigen::Vector3d origin2,
    const double radius2,
    double theta = 0.0)
  {
    // after https://gamedev.stackexchange.com/a/75775
    double d = (origin2 - origin1).norm();
    if (d < 1e-10) {
      return origin1;
    }
    if (d > radius1 + radius2) {
      throw std::runtime_error(
        "No intersection findable between the two spheres.");
    }
    double h = 0.5 + (radius1 * radius1 - radius2 * radius2) / (2 * d * d);
    Eigen::Vector3d intersectionCenter = origin1 + h * (origin2 - origin1);
    double intersectionRadius = sqrt(radius1 * radius1 - h * h * d * d);
    Eigen::Vector3d intersectionNormal = (origin2 - origin1) / d;

    // find non-parallel axis to the intersection normal
    Eigen::Index minRow;
    intersectionNormal.cwiseAbs().minCoeff(&minRow);

    Eigen::Vector3d nonParallelAxis = Eigen::Vector3d::Zero();
    nonParallelAxis[minRow] = 1.0;

    Eigen::Vector3d intersectionTangent =
      nonParallelAxis.cross(intersectionNormal).normalized();
    Eigen::Vector3d intersectionBinormal =
      intersectionTangent.cross(intersectionNormal);

    return intersectionCenter +
           intersectionRadius * (intersectionTangent * std::cos(theta) +
                                 intersectionBinormal * std::sin(theta));
  }

  static inline bool segmentIntersectsTriangle(
    const Eigen::Vector3d rayOrigin,
    const Eigen::Vector3d rayTarget,
    const Eigen::Vector3d vertex0,
    const Eigen::Vector3d vertex1,
    const Eigen::Vector3d vertex2,
    Eigen::Vector3d& outIntersectionPoint,
    const std::function<Eigen::Vector3d(Eigen::Vector3d)>& pbc,
    const double EPSILON = 1e-6)
  {
    // TODO: check & fix PBC usage
    // Möller-Trumbore intersection algorithm, see
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    // Eigen::Vector3d edge1 = pbc(vertex1 - vertex0);
    // Eigen::Vector3d edge2 = pbc(vertex2 - vertex0);
    Eigen::Vector3d edge1 = (vertex1 - vertex0);
    Eigen::Vector3d edge2 = (vertex2 - vertex0);

    // assume target and origin are in the same periodic image
    Eigen::Vector3d direction = rayTarget - rayOrigin;
    Eigen::Vector3d dir_norm = direction.normalized();

    Eigen::Vector3d h = dir_norm.cross(edge2);
    double a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON) {
      return false; // This ray is parallel to this triangle.
    }

    Eigen::Vector3d s = (rayOrigin - vertex0);
    // Eigen::Vector3d s = pbc(rayOrigin - vertex0);
    double f = 1.0 / a;
    double u = f * s.dot(h);

    if (u < 0.0 || u > 1.0) {
      return false;
    }

    Eigen::Vector3d q = s.cross(edge1);
    double v = f * dir_norm.dot(q);

    if (v < 0.0 || u + v > 1.0) {
      return false;
    }

    // At this stage we can compute t to find out where the intersection point
    // is on the line.
    double t = f * edge2.dot(q);

    if (t > EPSILON && t < sqrt(direction.dot(
                             direction))) { // here, adaption for segment, as of
      // https://stackoverflow.com/a/59475111/3909202
      outIntersectionPoint = rayOrigin + dir_norm * t;
      return true;
    }
    // This means that there is a line intersection but not a ray
    // intersection.
    return false;
  };

  static inline bool segmentIntersectsTriangle(
    const Eigen::Vector3d rayOrigin,
    const Eigen::Vector3d rayTarget,
    const Eigen::Vector3d vertex0,
    const Eigen::Vector3d vertex1,
    const Eigen::Vector3d vertex2,
    Eigen::Vector3d& outIntersectionPoint,
    const double EPSILON = 1e-6)
  {
    return segmentIntersectsTriangle(
      rayOrigin,
      rayTarget,
      vertex0,
      vertex1,
      vertex2,
      outIntersectionPoint,
      [](Eigen::Vector3d vec) { return vec; },
      EPSILON);
  }

}
}

#endif
