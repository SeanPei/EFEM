#ifndef EFEM_MATH_MATH_H
#define EFEM_MATH_MATH_H

#include "Eigen/Dense"
#include "math/MathType.h"
#include <iomanip>
#include <memory>
#include <mutex>
#include <sstream>
#include <vector>

namespace EFEM {
#define S_PI 3.14159265358979323846

class Math {
public:
    void singularValueDecomp(Eigen::Matrix3d& U, Eigen::Vector3d& SIGMA, Eigen::Matrix3d& V, const Eigen::Matrix3d& A, bool signed_svd = true);
    void kroneckerProduct(Matrix3333& dst, const EigenMatrix3& src1, const EigenMatrix3& src2);
    Eigen::Matrix4d poseMatrix(const std::vector<std::pair<std::string, std::vector<double>>>& pose);
    Eigen::Vector3d poseTransform(const std::vector<std::pair<std::string, std::vector<double>>>& pose, const Eigen::Vector3d& position);
    Eigen::Vector3d poseTransform(const Eigen::Matrix4d& pose_matrix, const Eigen::Vector3d& position);

    bool isReal(double n) { return !std::isnan(n) && !std::isinf(n); }
};
}

#endif